#include "mbconfig.h"
#include "SphericalQR.h"

#include "JacSubMatrix.h"
#include "fullmh.h"

#include <algorithm>
#include <cmath>
#include <random>


int main(void) {

    std::mt19937 generator;
    std::uniform_real_distribution<double> rand01(0., 10.);
    std::uniform_real_distribution<double> delta_rand01(-1.E-3, 1.E-3);
    std::random_device rd;
    generator.seed(rd());

    for (int i = 0; i < 10; i++) {
        Vec3 F0(rand01(generator), rand01(generator), rand01(generator));
        doublereal F0Norm = F0.Norm();
        Vec3 n0 = F0 / F0Norm;
        Mat3x3 Q;
        SphericalQR(F0, Q);
        Mat3x3 Qold = Q;

        Vec3 deltaF(delta_rand01(generator), delta_rand01(generator), delta_rand01(generator));
        deltaF = deltaF * F0Norm;
        // deltaF = deltaF - Q.GetCol(1) * Q.GetCol(1).Dot(deltaF);
        Vec3 delta_dir = deltaF / F0Norm;
        doublereal delta_dir1 = delta_dir.Dot(Q.GetCol(2));
        doublereal delta_dir2 = delta_dir.Dot(Q.GetCol(3));
        Vec3 F = F0 + deltaF;
        SphericalQR(F, Q, true, Qold);
        Mat3x3 Q1 = Q;

        // VERIFIED
        Vec3 deltaF1 = deltaF - Q.GetCol(1) * Q.GetCol(1).Dot(deltaF);
        Vec3 finite_diff_Q1 = (Q.GetCol(1) - Qold.GetCol(1));// / F0Norm;
        Vec3 finite_diff_Q2 = (Q.GetCol(2) - Qold.GetCol(2));// / F0Norm;
        Vec3 finite_diff_Q3 = (Q.GetCol(3) - Qold.GetCol(3));// / F0Norm;

        std::cout << "F0: " << F0 << std::endl<< std::endl;
        std::cout << "Qold:\n" << Qold.GetRow(1) << "\n"
            << Qold.GetRow(2) << "\n" << Qold.GetRow(3) << std::endl << std::endl;
        std::cout << "deltaF: " << deltaF << std::endl<< std::endl;
        std::cout << "F: " << F << std::endl<< std::endl;
        std::cout << "Q:\n" << Q.GetRow(1) << "\n"
            << Q.GetRow(2) << "\n" << Q.GetRow(3) << std::endl << std::endl;

        std::cout << "Z1 " << deltaF1 / F0Norm << " " << Q.GetCol(1) - Qold.GetCol(1) << std::endl;

        Vec3 der_Q1 = Qold.GetCol(2) * delta_dir1 + Qold.GetCol(3) * delta_dir2;
        Vec3 der_Q2 = -Qold.GetCol(1) * delta_dir1;
        Vec3 der_Q3 = -Qold.GetCol(1) * delta_dir2;
        std::cout << "Z2 " << Qold.GetCol(1) << std::endl;
        std::cout << "Z3 " << der_Q1 << " " << finite_diff_Q1 << std::endl;
        std::cout << "Z4 " << der_Q2 << " " << finite_diff_Q2 << std::endl;
        std::cout << "Z5 " << der_Q3 << " " << finite_diff_Q3 << std::endl;
        std::cout << std::endl;

        std::cout << "***************" << std::endl;
        ExpandableRowVector drealmodF;
        drealmodF.ReDim(3);
        doublereal real_modF = F0.Norm();
        drealmodF.Set(F0/real_modF, 1, 1);
        Vec3 dFNorm;
        drealmodF.Add(dFNorm, 1.);
        // std::cout << F.Norm() - F0Norm << " " << dFNorm.Dot(deltaF) << std::endl;

        // drealmodF.Set(Vec3(2., 3., 4.), 1, 1);
        // ExpandableMatrix dQ1;
        // dQ1.ReDim(3, 2);
        // dQ1.SetBlockDim(1, 3);
        // dQ1.SetBlockIdx(1, 1);
        // dQ1.SetBlockDim(2, 1);
        // dQ1.Set(Eye3 * 0., 1, 1, 1);
        // dQ1.SetCol(Vec3(0.0, 1.1, 0.), 1, 2, 1);
        // dQ1.Link(2, &drealmodF);
        // FullSubMatrixHandler X1(3, 3);
        // for (int i = 1; i <= 3; i++) {
        //     X1.PutRowIndex(i, i);
        //     X1.PutColIndex(i, i);
        // }
        // X1.Reset();
        // FullMatrixHandler X1M(3, 3);
        // X1M.Reset();
        //
        // dQ1.Add(X1, 1, 1.);
        // X1.AddTo(X1M);
        // std::cout << X1 << std::endl;
        // std::cout << ".....\n";
        // std::cout << X1M << std::endl;


        ExpandableMatrix dQ1, dQ2, dQ3;
        // dQ1.ReDim(3, 2);
        // dQ1.SetBlockDim(1, 3);
        // dQ1.SetBlockIdx(1, 1);
        // dQ1.SetBlockDim(2, 1);
        // dQ1.Set(Eye3/real_modF, 1, 1, 1);
        // dQ1.SetCol(-F/(real_modF*real_modF), 1, 2, 1);
        // dQ1.Link(2, &drealmodF);

        dQ1.ReDim(3, 1);
        dQ1.SetBlockDim(1,3);
        dQ1.SetBlockIdx(1, 1);
        dQ1.Set(Eye3/real_modF - F0.Tens(F0)/std::pow(real_modF, 3), 1, 1, 1);

        dQ2.ReDim(3, 1);
        dQ2.SetBlockDim(1, 3);
        dQ2.Set(-Q.GetCol(1).Tens(Q.GetCol(2)), 1, 1, 1);
        dQ2.Link(1, &dQ1);

        dQ3.ReDim(3, 1);
        dQ3.SetBlockDim(1, 3);
        dQ3.Set(-Q.GetCol(1).Tens(Q.GetCol(3)), 1, 1, 1);
        dQ3.Link(1, &dQ1);

        FullSubMatrixHandler X1(3, 3), X2(3, 3), X3(3,3);
        for (int ii = 1; ii <= 3; ii++) {
            X1.PutRowIndex(ii, ii);
            X1.PutColIndex(ii, ii);
            X2.PutRowIndex(ii, ii);
            X2.PutColIndex(ii, ii);
            X3.PutRowIndex(ii, ii);
            X3.PutColIndex(ii, ii);
        }
        X1.Reset();
        X2.Reset();
        X3.Reset();
        FullMatrixHandler X1M(3, 3);
        FullMatrixHandler X2M(3, 3);
        FullMatrixHandler X3M(3, 3);
        X1M.Reset();
        X2M.Reset();
        X3M.Reset();


        dQ1.AddTo(X1, 1, 1.);
        X1.AddTo(X1M);
        dQ2.AddTo(X2, 1, 1.);
        X2.AddTo(X2M);
        dQ3.AddTo(X3, 1, 1.);
        X3.AddTo(X3M);

        // std::cout << X1M << std::endl;
        Mat3x3 X1MM, X2MM, X3MM;
        for (int ii=1; ii<=3; ii++) {
            for (int j=1; j<=3; j++) {
                X1MM(ii, j) = X1M(ii, j);
                X2MM(ii, j) = X2M(ii, j);
                X3MM(ii, j) = X3M(ii, j);
            }

        }

        std::cout << X1MM * deltaF << " " << finite_diff_Q1 << std::endl;
        std::cout << X2MM * deltaF << " " << finite_diff_Q2 << std::endl;
        std::cout << X3MM * deltaF << " " << finite_diff_Q3 << std::endl;


        std::cout << "================" << std::endl;
    }
}


        // Vec3 F0(rand01(generator), rand01(generator), rand01(generator));
        // Vec3 n0 = F0 / F0.Norm();
        // Mat3x3 Q;
        // SphericalQR(F0, Q);
        // Mat3x3 Qold = Q;
        // Vec3 deltaF0(delta_rand01(generator), delta_rand01(generator), delta_rand01(generator));
        // Vec3 deltaF = deltaF0 - Q.GetCol(1) * Q.GetCol(1).Dot(deltaF0);
        // doublereal F0Norm = F0.Norm();
        // Vec3 delta_dir = deltaF / F0Norm;
        // doublereal delta_dir1 = delta_dir.Dot(Q.GetCol(2));
        // doublereal delta_dir2 = delta_dir.Dot(Q.GetCol(3));
        // Vec3 F = F0 + deltaF;
        // SphericalQR(F, Q, true, Qold);
        // // VERIFIED
        // std::cout << deltaF / F0.Norm() << " " << Q.GetCol(1) - Qold.GetCol(1) << std::endl;
        // std::cout << Qold.GetCol(2) * delta_dir1 + Qold.GetCol(3) * delta_dir2 << " " << Q.GetCol(1) - Qold.GetCol(1) << std::endl;
        // Vec3 finite_diff_Q2 = (Q.GetCol(2) - Qold.GetCol(2));// / F0Norm;
        // Vec3 finite_diff_Q3 = (Q.GetCol(3) - Qold.GetCol(3));// / F0Norm;
        // Vec3 der_Q2 = -Qold.GetCol(1) * delta_dir1;
        // Vec3 der_Q3 = -Qold.GetCol(1) * delta_dir2;
        // std::cout << Qold.GetCol(1) << std::endl;
        // std::cout << finite_diff_Q2 << " " << der_Q2 << std::endl;
        // std::cout << finite_diff_Q3 << " " << der_Q3 << std::endl;
        // std::cout << "---" << std::endl;
        // ExpandableRowVector drealmodF;
        // ExpandableMatrix dQ1, dQ2, dQ3;
        // drealmodF.ReDim(3);
        // dQ1.ReDim(3, 2);
        // dQ1.SetBlockDim(1, 3);
        // dQ1.SetBlockIdx(1, 1);
        // dQ1.SetBlockDim(2, 1);
        // doublereal real_modF = F0.Norm();
        // drealmodF.Set(F0/real_modF, 1, 1);
        // dQ1.Set(Eye3/real_modF, 1, 1, 1);
        // dQ1.SetCol(-F/(real_modF*real_modF), 1, 2, 1);
        // dQ1.Link(2, &drealmodF);
        //
        // dQ2.ReDim(3, 1);
        // dQ2.SetBlockDim(1, 3);
        // dQ2.Set(-Q.GetCol(1).Tens(Q.GetCol(2)), 1, 1, 1);
        // dQ2.Link(1, &dQ1);
        //
        // dQ3.ReDim(3, 1);
        // dQ3.SetBlockDim(1, 3);
        // dQ3.Set(-Q.GetCol(1).Tens(Q.GetCol(3)), 1, 1, 1);
        // dQ3.Link(1, &dQ1);
        //
        // Mat3x3 X1, X2, X3;
        // dQ2.Add(X2, 1, 1.);
        // dQ3.Add(X3, 1, 1.);
        // std::cout << finite_diff_Q2 << " " << X2 * deltaF0 << std::endl;
        // std::cout << finite_diff_Q3 << " " << X3 * deltaF0 << std::endl;
        // std::cout << std::endl << std::endl;
