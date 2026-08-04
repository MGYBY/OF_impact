MRF.correctBoundaryVelocity(U);

// Unit vector parallel to the periodic/downslope direction.
const vector eStream(1, 0, 0);

const dimensionedVector gParallel
(
    "gParallel",
    g.dimensions(),
    (g.value() & eStream)*eStream
);

const phase& airPhase = mixture.phases()["air"];

// Full gravity in mud and water.
// Only bed-normal gravity in pure air.
tmp<volVectorField> tGravityForce
(
    rho*g
  - airPhase*airPhase.rho()*gParallel
);

fvVectorMatrix UEqn
(
    fvm::ddt(rho, U)
  + fvm::div(rhoPhi, U)
  + MRF.DDt(rho, U)
 ==
    fvOptions(rho, U)
  + mixture.divTauMF(U)
  + tGravityForce()
);

UEqn.relax();
fvOptions.constrain(UEqn);

spSolverU->solve
(
    UEqn
 ==
    fvc::reconstruct
    (
        (
            mixture.surfaceTensionForce()
          - fvc::snGrad(p)
        )*mesh.magSf()
    )
);

fvOptions.correct(U);
