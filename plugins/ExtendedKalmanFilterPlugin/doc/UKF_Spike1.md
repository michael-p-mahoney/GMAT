# Unscented Kalman Filter — Spike 1 Design Note

Branch: `feature/unscented-kalman-filter`

## Goal of Spike 1

Prove that GMAT's propagation subsystem can **re-seed an arbitrary state and
re-propagate it from inside the estimator**. This is the single capability a true
UKF time update needs and was the only real unknown in the effort estimate. If
this works cleanly, the rest of the UKF is well-understood linear algebra layered
on the existing (already square-root) EKF infrastructure.

## Key architectural findings

1. **Estimator hierarchy.** `ExtendedKalmanFilter : SeqEstimator : Estimator`.
   `SeqEstimator` (~3800 lines) owns the finite state machine, process noise,
   data-file/JSON/`.mat` output, covariance reporting (incl. VNB), coordinate /
   epsilon conversions, and the smoother forward-pass hooks. `ExtendedKalmanFilter`
   (~1400 lines) contributes only the filter math: `TimeUpdate`, `SetupMeas`,
   `ComputeObs` (builds analytic `H`), `ComputeGain`, `UpdateElements`. A UKF
   inherits everything in `SeqEstimator` unchanged.

2. **The architecture was designed for a UKF.** Both estimator factories already
   carry commented-out `UnscentedKalmanFilter` hooks
   (`ExtendedKalmanFilterFactory.cpp`, `EstimatorFactory.cpp`).

3. **The EKF is already square-root.** Covariance is stored/updated as `sqrtP`
   via Cholesky / `thinQR`. The natural target is the **Square-Root UKF**
   (Van der Merwe & Wan), which reuses this plumbing directly.

4. **The estimator does NOT own the propagation loop — the command does.**
   `RunEstimator::Execute()` runs a finite state machine
   (`PROPAGATING → CALCULATING → ESTIMATING`). In the `PROPAGATING` state,
   `RunEstimator::Propagate()` calls `Step(dt)` where `dt = theEstimator->GetTimeStep()`,
   re-seeding the integrator from the spacecraft via `fm[i]->UpdateFromSpaceObject()`
   when `theEstimator->ResetState()` is true. The estimator only *reports* the step
   size; the command advances the single reference trajectory (carrying the STM via
   variational equations).

   → Consequence: to propagate `2n+1` sigma points, the UKF must drive a
   propagator **directly**, outside the command's one-step-at-a-time loop.

5. **Direct propagation is already done elsewhere in estimation code.**
   `EstimationRootFinder` (light-time / event location) drives
   `propagator->GetPropagator()->Step(dt)` then `ode->UpdateSpaceObject(newEpoch)`
   directly (see `EstimationRootFinder.cpp:185, 275`). Signal/measurement code does
   the same (`MeasureModel.cpp:1540`, `SignalBase.cpp:2414`). This is the proven
   template the UKF reuses.

## Re-seed / re-propagate mechanism (the spike)

Per sigma point, from the base `Estimator`'s held `propagators[0]` (`PropSetup*`):

```
esm.SetEstimationState(seed);     // write the sigma point into the state
esm.MapVectorToObjects();         // push it onto the spacecraft objects
ode->UpdateFromSpaceObject();     // re-seed the integrator from the objects
p->Step(dt);                      // propagate by dt
ode->UpdateSpaceObject(newEpoch); // read the propagated state back to objects
GmatState after = esm.GetEstimationState();
```

`PropSetup::GetPropagator()/GetODEModel()/GetPropStateManager()` provide the
handles; `Propagator::Step(Real)` and `ODEModel::UpdateFromSpaceObject()/
UpdateSpaceObject()` are the re-seed/step/read API.

## What was built in this spike commit

- `UnscentedKalmanFilter` class registered under its own type name (clone of the
  EKF math, so `Create UnscentedKalmanFilter` builds a working, runnable filter
  whose numerical behavior currently equals the EKF). Wired into
  `ExtendedKalmanFilterFactory` and the plugin `CMakeLists.txt`.
- `Alpha` / `Beta` / `Kappa` tuning parameters (scriptable Real fields) plus the
  Van der Merwe sigma-point weights computed in `CompleteInitialization()`.
- `GenerateSigmaPoints()` — the `2n+1` unscented sigma points from `(mean, sqrtP)`.
  Pure linear algebra, low risk.
- `PropagateSigmaPoints(dt)` — the instrumented experiment: saves the reference
  state, generates sigma points, re-seeds + steps + reads each one, recombines the
  weighted mean, logs it against the central point, then restores the reference so
  the live filter is unperturbed.
- A dedicated, fully-initialized sigma-point propagator (`BuildSigmaPropagator()` →
  `sigmaProp`) — see "Why the first run crashed" below.
- Gated behind `DEBUG_UKF_SPIKE`. With the flag off, the class is a clean EKF
  baseline; with it on, `CompleteInitialization()` builds `sigmaProp` and
  `TimeUpdate()` runs the experiment and logs.

## Why the first run crashed (and the fix)

The first instrumented run segfaulted on sigma point #0, immediately at the
re-seed/step sequence. Root cause: the spike stepped `propagators[0]`, the base
`Estimator`'s configured `PropSetup`. **That object is never `Initialize()`d for
standalone stepping.** During a run it is `RunEstimator` that clones the propagator
into its own private list (`p[]`/`fm[]`) and initializes and steps *those*
(`RunEstimator::Execute`, `p[i]->Initialize()`); the estimator's own copy keeps no
built propagation-state vector. Calling `ODEModel::UpdateFromSpaceObject()` /
`Propagator::Step()` on it dereferences an unallocated state → crash. (This is also
why `EstimationRootFinder` can `Step()` with no setup: its `propagator` member *is*
one of the command's already-initialized propagators.)

Fix: `UnscentedKalmanFilter::BuildSigmaPropagator()` (called from
`CompleteInitialization()`, gated by `DEBUG_UKF_SPIKE`) clones `propagators[0]` into
the owned `sigmaProp` member and prepares it for stepping:

```
ObjectArray objs;  esm.GetStateObjects(objs, Gmat::SPACEOBJECT);   // esm's spacecraft
sigmaProp = (PropSetup*)propagators[0]->Clone();
sigmaProp->SetSolarSystem(solarSystem);
psm = sigmaProp->GetPropStateManager();
psm->SetObject(sc);  sigmaProp->SetRefObject(sc, ...);             // for each space obj
psm->SetProperty("CartesianState");  psm->SetProperty("STM");
sigmaProp->PrepareInternals();        // builds state, BuildModelFromMap, inits propagator
```

The critical detail: the clone's `PropagationStateManager` is wired to the **same**
`SpaceObject`s that `esm` maps the estimation state onto. That makes the round-trip
consistent — `esm.MapVectorToObjects()` writes the seed onto those spacecraft,
`ode->UpdateFromSpaceObject()` reads them into `sigmaProp`, we `Step()`,
`ode->UpdateSpaceObject()` writes back, and `esm.GetEstimationState()` reads the
result. `PropagateSigmaPoints()` now steps `sigmaProp` instead of `propagators[0]`.

## Spike 1 result — PROVEN

`Ex_UKF_Spike_GpsPosVec.script`, 3 experiments at consecutive time updates,
dt = 60 s, 7-state (Cartesian + Cd), 15 sigma points, sigma-point prop state dim 55
(6 + 7×7 STM):

- All 3 experiments ran to completion; estimation finished with `EXIT=0`.
- Recombined unscented-transform mean matches the central propagated point to ~8–9
  significant figures.
- Nonlinearity signature (UT mean − central) ≈ 1e-7 km in position, 1e-9 km/s in
  velocity over the 60 s LEO step — correct magnitude for near-linear dynamics, and
  nonzero, confirming each sigma point was independently re-propagated through the
  real force model rather than echoed.

Conclusion: GMAT's propagation subsystem can re-seed and re-propagate an arbitrary
state from inside the estimator. The dominant unknown in the effort estimate is
resolved; the remaining work is the well-understood SR-UKF linear algebra.

## How to run the spike

1. Build the `EKF` plugin target.
2. Define `DEBUG_UKF_SPIKE` at the top of `UnscentedKalmanFilter.cpp` and rebuild.
3. Take an existing EKF estimation script, change `Create ExtendedKalmanFilter` to
   `Create UnscentedKalmanFilter`, and run it.
4. Inspect the `[UKF Spike]` log lines: the recombined sigma-point mean, the
   central propagated point, and their difference (the nonlinearity signature).

## Open items to resolve during the build/run loop

- **Coordinate frame consistency.** The estimation state may be Cartesian
  solve-for, Keplerian, or epsilon-scaled. Sigma points must be generated and
  propagated in a consistent frame (Cartesian MJ2000Eq for the dynamics). The
  draft seeds in the native estimation state; confirm/round-trip first, then refine.
- **STM in the prop vector.** The integrator carries the STM (variational
  equations). The UKF does not need it; once the round-trip is confirmed, disable
  it for sigma-point propagation to cut cost.
- **State restoration fidelity.** Verify the propagator is recovered exactly so the
  (currently EKF) update following the experiment is bit-for-bit unaffected.
- **Runtime.** Expect ~`(2n+1)×` the per-step propagation cost vs. the EKF.

## After the spike

Replace the STM-based covariance map in `TimeUpdate()` and the analytic-`H`
measurement update (`ComputeObs`/`ComputeGain`) with the unscented transform
(SR-UKF), reusing the inherited reporting, data-file, JSON, VNB, process-noise,
and conversion machinery. Unscented smoother (URTSS) is a follow-on; v1 is
filter-only.
