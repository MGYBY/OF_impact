if
(
    !runTime.restart()
 || isType<twoPhaseChangeModels::noPhaseChange>(phaseChange)
)
{
    if (correctPhi)
    {
        rAU = new volScalarField
        (
            IOobject
            (
                "rAU",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar(dimTime/dimDensity, 1)
        );

        correctUphiBCs(U, phi, true);

        CorrectPhi
        (
            phi,
            U,
            // p_rgh,
            p,
            surfaceScalarField("rAUf", fvc::interpolate(rAU())),
            geometricZeroField(),
            pressureReference,
            pimple
        );
    }
    else
    {
        correctUphiBCs(U, phi, true);

        CorrectPhi
        (
            phi,
            U,
            // p_rgh,
            p,
            dimensionedScalar(dimTime/rho.dimensions(), 1),
            geometricZeroField(),
            pressureReference,
            pimple
        );
    }
}

#include "continuityErrs.H"
