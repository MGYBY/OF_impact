Info<< "Reading physical pressure field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"

multiphaseMixture mixture(U, phi);

// Need to store rho for ddt(rho, U).
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT
    ),
    mixture.rho()
);
rho.oldTime();

// The parent solver's turbulence stress is intentionally inactive in
// UEqn.H.  The object is retained for API/algorithm compatibility.
autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, mixture)
);

// Physical-pressure form: retain g, but do not construct hRef, gh,
// ghf, or p_rgh.
#include "readGravitationalAcceleration.H"

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    pimple.dict(),
    pRefCell,
    pRefValue
);

if (p.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
}

mesh.setFluxRequired(p.name());

#include "createMRF.H"
#include "createFvOptions.H"
