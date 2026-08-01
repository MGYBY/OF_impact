/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

Application
    rheoMultiFluidInterFoamP

Description
    Physical-pressure variant of rheoMultiFluidInterFoam for n
    incompressible Newtonian or non-Newtonian phases.  It advances the
    mechanical/static pressure p and retains gravity as +rho*g.  This
    permits ordinary translational cyclic conditions in a bed-aligned
    periodic inclined-flow domain when gravity has a streamwise
    component.

    VOF transport, phase-specific surface tension, optional curvature
    smoothing, SIMPLEC coupling, dynamic-mesh data structures and
    RheoTool constitutive equations are inherited from the parent.

\*---------------------------------------------------------------------------*/

#include "sparseMatrixSolvers.H" // Avoid namespace Foam

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "multiphaseMixture.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "CorrectPhi.H"

#include "adjustCorrPhi.H"
#include "ppUtilInterface.H"

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createDyMControls.H"
    #include "createFields.H"
    #include "createSolver.H"
    #include "initCorrectPhi.H"
    #include "createUfIfPresent.H"
    #include "createPPutil.H"

    turbulence->validate();

    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    const surfaceScalarField& rhoPhi(mixture.rhoPhi());

    Info<< "\nUsing physical pressure p with gravity source +rho*g."
        << "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readDyMControls.H"
        #include "CourantNo.H"
        #include "alphaCourantNo.H"
        #include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        while (pimple.loop())
        {
            if (pimple.firstPimpleIter() || moveMeshOuterCorrectors)
            {
                scalar timeBeforeMeshUpdate = runTime.elapsedCpuTime();

                mesh.update();

                if (mesh.changing())
                {
                    Info<< "Execution time for mesh.update() = "
                        << runTime.elapsedCpuTime() - timeBeforeMeshUpdate
                        << " s" << endl;

                    MRF.update();

                    if (correctPhi)
                    {
                        phi = mesh.Sf() & Uf();
                        #include "correctPhi.H"
                        fvc::makeRelative(phi, U);
                    }

                    mixture.correct();

                    if (checkMeshCourantNo)
                    {
                        #include "meshCourantNo.H"
                    }
                }
            }

            mixture.solve();
            rho = mixture.rho();

            #include "UEqn.H"

            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;
    return 0;
}

// ************************************************************************* //
