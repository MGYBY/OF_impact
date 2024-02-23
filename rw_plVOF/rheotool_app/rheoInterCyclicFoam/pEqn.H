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
    // volVectorField HbyA(constrainHbyA(rAU()*UEqn.H(), U, p_rgh));
    volVectorField HbyA(constrainHbyA(rAU()*UEqn.H(), U, p));
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::flux(HbyA)
      + MRF.zeroFilter(fvc::interpolate(rho*rAU())*fvc::ddtCorr(U, phi, Uf))
    );
    MRF.makeRelative(phiHbyA);

    // if (p_rgh.needReference())
    if (p.needReference())
    {
        fvc::makeRelative(phiHbyA, U);
        // adjustCorrPhi(phiHbyA, U, p_rgh);  
        adjustCorrPhi(phiHbyA, U, p);
        fvc::makeAbsolute(phiHbyA, U); 
    }

    surfaceScalarField phig
    (
        (
            mixture.surfaceTensionForce()
          - ghf*fvc::snGrad(rho)
        )*rAUf*mesh.magSf()
    );

    phiHbyA += phig;

    tmp<volScalarField> rAtU(rAU);
     
    if (simplec)
    {
      rAtU = 1.0/(1.0/rAU() - UEqn.H1());
    //   phiHbyA += fvc::interpolate(rAtU() - rAU())*fvc::snGrad(p_rgh)*mesh.magSf();
    //   HbyA -= (rAU() - rAtU()) * fvc::reconstruct( fvc::snGrad(p_rgh) * mesh.magSf() );

      phiHbyA += fvc::interpolate(rAtU() - rAU())*fvc::snGrad(p)*mesh.magSf();
      HbyA -= (rAU() - rAtU()) * fvc::reconstruct( fvc::snGrad(p) * mesh.magSf() );
      
      // Update the pressure BCs to ensure flux consistency
    //   constrainPressure(p_rgh, U, phiHbyA, rAtU(), MRF);
    constrainPressure(p, U, phiHbyA, rAtU(), MRF);
    }
    else
    {
      // Update the pressure BCs to ensure flux consistency
    //   constrainPressure(p_rgh, U, phiHbyA, rAUf, MRF);  
        constrainPressure(p, U, phiHbyA, rAUf, MRF);
    }

    // Cache the phase change pressure source
    // fvScalarMatrix Sp_rgh(phaseChange.Sp_rgh(rho, gh, p_rgh));
    fvScalarMatrix Sp(phaseChange.Sp(rho, gh, p));

    while (pimple.correctNonOrthogonal())
    {
        // fvScalarMatrix p_rghEqn
        fvScalarMatrix pEqn
        (
            // fvc::div(phiHbyA) - fvm::laplacian(rAtU(), p_rgh)
            fvc::div(phiHbyA) - fvm::laplacian(rAtU(), p)
         == Sp
        );

        // p_rghEqn.setReference
        pEqn.setReference
        (
            pressureReference.refCell(),
            // getRefCellValue(p_rgh, pressureReference.refCell())
            getRefCellValue(p, pressureReference.refCell())
        );

        // spSolverprgh->solve(p_rghEqn); 
        spSolverp->solve(pEqn); 

        if (pimple.finalNonOrthogonalIter())
        {
            // phi = phiHbyA + p_rghEqn.flux();
            phi = phiHbyA + pEqn.flux();

            // p_rgh.relax();
            p.relax();
            
            if (simplec)
            {
               surfaceScalarField rAtUf("rAtUf", fvc::interpolate(rAtU()));
            //    U = HbyA + rAU()*fvc::reconstruct(phig/rAUf) + rAtU()*fvc::reconstruct(p_rghEqn.flux()/rAtUf);
                U = HbyA + rAU()*fvc::reconstruct(phig/rAUf) + rAtU()*fvc::reconstruct(pEqn.flux()/rAtUf);
            }
            else
            {
            //    U = HbyA + rAU()*fvc::reconstruct((phig + p_rghEqn.flux())/rAUf);
                U = HbyA + rAU()*fvc::reconstruct((phig + pEqn.flux())/rAUf);
            }
            
            U.correctBoundaryConditions();
            fvConstraints.constrain(U);
        }
    }

    #include "continuityErrs.H"

    // Correct Uf if the mesh is moving
    fvc::correctUf(Uf, U, phi);

    // Make the fluxes relative to the mesh motion
    fvc::makeRelative(phi, U);

    // p == p_rgh + rho*gh;

    // if (p_rgh.needReference())
    if (p.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pressureReference.refValue()
          - getRefCellValue(p, pressureReference.refCell())
        );
        // p_rgh = p - rho*gh;
    }

    if (!correctPhi)
    {
        rAU.clear();
    }
}
