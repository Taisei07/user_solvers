//.H形式のファイルはヘッダファイルと呼ばれる
Info<< "Reading transportProperties\n" << endl;

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

dimensionedScalar nu//動粘度nuをスカラーとして定義
(
    "nu",
    dimViscosity,
    transportProperties.lookup("nu")
);

dimensionedScalar q//電荷密度qをスカラーとして定義
(
    "q",
    transportProperties.lookup("q")
);

dimensionedScalar epsilon//誘電率epsilonをスカラーとして定義
(
    "epsilon",
    transportProperties.lookup("epsilon")
);

Info<< "Reading field p\n" << endl;
volScalarField p//圧力pをスカラーとして定義
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

Info<< "Reading field V\n" << endl;//↓オリジナルで追加
volScalarField V//電位Vをスカラーとして定義。
(
    IOobject
    (
        "V",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


Info<< "Reading field U\n" << endl;
volVectorField U//速度Uをベクトルとして定義
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


#include "createPhi.H"//phiは他のファイルから持ってきている


label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, mesh.solutionDict().subDict("PISO"), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());
