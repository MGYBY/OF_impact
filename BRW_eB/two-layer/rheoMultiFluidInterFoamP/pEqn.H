{
    if (correctPhi)
    {
        rAU.ref() = 1.0/UEqn.A();
    }
    else
    {
        rAU = 1.0/UEqn.A();
    }

    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU()));
    volVectorField HbyA(constrainHbyA(rAU()*UEqn.H(), U, p));
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::flux(HbyA)
      + MRF.zeroFilter
        (
            fvc::interpolate(rho*rAU())*fvc::ddtCorr(U, phi, Uf)
        )
    );
    MRF.makeRelative(phiHbyA);

    if (p.needReference())
    {
        fvc::makeRelative(phiHbyA, U);
        adjustCorrPhi(phiHbyA, U, p);
        fvc::makeAbsolute(phiHbyA, U);
    }

    // Gravity is already contained in UEqn.H() through +rho*g.
    // Only the capillary contribution remains as an explicit face
    // flux in the pressure corrector.
    surfaceScalarField phiSt
    (
        "phiSt",
        mixture.surfaceTensionForce()*rAUf*mesh.magSf()
    );

    phiHbyA += phiSt;

    tmp<volScalarField> rAtU(rAU);

    if (simplec)
    {
        rAtU = 1.0/(1.0/rAU() - UEqn.H1());

        phiHbyA +=
            fvc::interpolate(rAtU() - rAU())
           *fvc::snGrad(p)*mesh.magSf();

        HbyA -=
            (rAU() - rAtU())
           *fvc::reconstruct(fvc::snGrad(p)*mesh.magSf());

        constrainPressure(p, U, phiHbyA, rAtU(), MRF);
    }
    else
    {
        constrainPressure(p, U, phiHbyA, rAUf, MRF);
    }

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(rAtU(), p, "laplacian(p|rAtU)")
         == fvc::div(phiHbyA)
        );

        pEqn.setReference(pRefCell, pRefValue);
        spSolverp->solve(pEqn);

        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA - pEqn.flux();
            p.relax();

            if (simplec)
            {
                surfaceScalarField rAtUf
                (
                    "rAtUf",
                    fvc::interpolate(rAtU())
                );

                U =
                    HbyA
                  + rAU()*fvc::reconstruct(phiSt/rAUf)
                  - rAtU()*fvc::reconstruct(pEqn.flux()/rAtUf);
            }
            else
            {
                U =
                    HbyA
                  + rAU()*fvc::reconstruct
                    (
                        (phiSt - pEqn.flux())/rAUf
                    );
            }

            U.correctBoundaryConditions();
            fvOptions.correct(U);
        }
    }

    #include "continuityErrs.H"

    // Correct Uf if the mesh is moving.
    fvc::correctUf(Uf, U, phi);

    // Make the fluxes relative to the mesh motion.
    fvc::makeRelative(phi, U);

    // In a closed/cyclic domain p is defined up to a constant.
    if (p.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
    }

    if (!correctPhi)
    {
        rAU.clear();
    }
}
