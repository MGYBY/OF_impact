MRF.correctBoundaryVelocity(U);

fvVectorMatrix UEqn
(
    fvm::ddt(rho, U)
  + fvm::div(rhoPhi, U)
  + MRF.DDt(rho, U)
//+ turbulence->divDevRhoReff(rho, U)
 ==
    fvOptions(rho, U)
  + mixture.divTauMF(U)

    // OpenFOAM stores g as the physical acceleration vector.  From
    // p_rgh = p - rho*(g & x),
    // -grad(p_rgh) - (g & x)grad(rho) = -grad(p) + rho*g.
  + rho*g
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
