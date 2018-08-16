/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    icoFoam

Description
    Transient solver for incompressible, laminar flow of Newtonian fluids.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"  //include：ファイルを取り込んでいる。""で囲まれている部分は今のフォルダのもの。
#include "pisoControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//プログラムの本体はmain関数。プログラムを起動するとmain関数中の文が順次実行される。

int main(int argc, char *argv[])//argcには引数の個数・acgvには引数の実際の値が入る。
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    pisoControl piso(mesh);

    #include "createFields.H"
    #include "initContinuityErrs.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;//Time=〜と表示させる。

        #include "CourantNo.H"//CourantNo.Hを取り込んでいる。

        // Momentum predictor（運動量予測）

        fvScalarMatrix VEqn//fvScalarMatrixクラスからVEqnインスタンスを作成。*ラプラス方程式を用いて電位を求める
        (
            fvm::laplacian(V)
            == -q/epsilon
        );

        fvVectorMatrix UEqn//fvVectorMatrixクラスからUEqnインスタンスを作成。
        (
            fvm::ddt(U)//∂U/∂t
          + fvm::div(phi, U)//div(UU)
          - fvm::laplacian(nu, U)//∇^2(νU)
          + fvm::div(q, V)//クーロン力の項を追加（オリジナル）
        );

        if (piso.momentumPredictor())
        {
            solve(UEqn == -fvc::grad(p));//UEqn=-∇(p)
        }

        // --- PISO loop
        /*piso法の計算手順
        1.運動方程式を解き、仮の速度を求める
        2.圧力方程式を解き、圧力を求める
        3.速度を更新する
        4.圧力の計算および速度の更新を指定回数だけ繰り返す（通常は2回）*/

        while (piso.correct())
        {
            volScalarField rAU(1.0/UEqn.A());
            volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
            surfaceScalarField phiHbyA
            (
                "phiHbyA",
                fvc::flux(HbyA)
              + fvc::interpolate(rAU)*fvc::ddtCorr(U, phi)
            );

            adjustPhi(phiHbyA, U, p);//流速が成り行きで決まる境界の流束を、質量が保存するように調整。

            // Update the pressure BCs to ensure flux consistency：流れの連続性を確実にするために圧力BCを更新する
            constrainPressure(p, U, phiHbyA, rAU);

            // Non-orthogonal pressure corrector loop：非直行圧力補正ループ
            while (piso.correctNonOrthogonal())
            {
                // Pressure corrector：圧力補正

                fvScalarMatrix pEqn
                (
                    fvm::laplacian(rAU, p) == fvc::div(phiHbyA)
                );

                pEqn.setReference(pRefCell, pRefValue);

                pEqn.solve(mesh.solver(p.select(piso.finalInnerIter())));

                if (piso.finalNonOrthogonalIter())
                {
                    phi = phiHbyA - pEqn.flux();
                }
            }

            #include "continuityErrs.H"

            U = HbyA - rAU*fvc::grad(p);
            U.correctBoundaryConditions();
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
