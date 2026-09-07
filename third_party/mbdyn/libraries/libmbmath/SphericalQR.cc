/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 2003-2023
 *
 * This code is a partial merge of HmFe and MBDyn.
 *
 * Pierangelo Masarati  <pierangelo.masarati@polimi.it>
 * Paolo Mantegazza     <paolo.mantegazza@polimi.it>
 * Marco Morandini  <morandini@aero.polimi.it>
 *
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 * via La Masa, 34 - 20156 Milano, Italy
 * http://www.aero.polimi.it
 *
 * Changing this copyright notice is forbidden.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation (version 2 of the License).
 *
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include "mbconfig.h"
#include "matvec3.h"
#include "Rot.hh"

#include <algorithm>
#include <cmath>
#include <numbers> 

// #ifdef HAVE_FENV_H
// #include <fenv.h>
// static void __attribute__ ((constructor))
// trapfpe ()
// {
//         /* Enable some exceptions.  At startup all exceptions are masked.  */
//
//         feenableexcept(FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);
// }
// #endif

/*
 * Compute the SVD of the symmetric matrix A, U Sigma V^T = svd(A)
 *
 * a[0] = A(1,1)
 * a[1] = A(1,2)
 * a[2] = A(2,2)
 *
 * u[0][0] = U(1,1)
 * u[0][1] = U(1,2)
 * u[1][0] = U(2,1)
 * u[1][1] = U(2,2)
 *
 * s[0] = Sigma(1, 2)
 * s[1] = Sigma(2, 2)
 *
 * u[0][0] = U(1,1)
 * u[0][1] = U(1,2)
 * u[1][0] = U(2,1)
 * u[1][1] = U(2,2)
 *
 * Formulae taken from https://lucidar.me/en/mathematics/singular-value-decomposition-of-a-2x2-matrix/
 *
 */

// THE SVD NEEDS TO BE VALIDATED
//
// void svd2x2S(const doublereal a[3], doublereal u[1][1], doublereal s[2], doublereal v[1][1]) {
//     doublereal theta = 0.5 * std::atan2(2. * a[1] * (a[0] + a[2]), a[0] * a[0] - a[2] * a[2]);
//     u[0][0] = std::cos(theta);
//     u[0][1] = -std::sin(theta);
//     u[1][0] = -u[0][1];
//     u[1][1] = u[0][0];
//
//     doublereal S1 = a[0] * a[0] + 2. * a[1] * a[1] + a[2] * a[2];
//     doublereal S2 = std::sqrt(std::pow(a[0] * a[0] - a[2] * a[2], 2) + 4. * (a[1] * (a[0] + a[2])));
//     s[0] = std::sqrt((S1 + S2) / 2.);
//     s[1] = std::sqrt((S1 - S2) / 2.);
//
//     doublereal phi = theta;
//     doublereal cphitheta = std::cos(phi);
//     doublereal sphitheta = std::sin(phi);
//     doublereal s11 = (a[0] * cphitheta + a[2] * sphitheta) * cphitheta + ( a[1] * cphitheta + a[2] * sphitheta) * sphitheta;
//     doublereal s22 = (a[0] * sphitheta - a[2] * cphitheta) * sphitheta + (-a[1] * sphitheta + a[2] * cphitheta) * cphitheta;
//     if (s11) s11 = std::copysign(1., s11);
//     if (s22) s22 = std::copysign(1., s22);
//
//     v[0][0] =  s11 * cphitheta;
//     v[0][1] = -s22 * sphitheta;
//     v[1][0] =  s11 * sphitheta;
//     v[1][1] =  s22 * cphitheta;
//
//     return;
// }

/*
 * Compute the determinant of the symmetric matrix A, d = det(A)
 *
 * a[0] = A(1,1)
 * a[1] = A(1,2)
 * a[2] = A(2,2)
*/
doublereal det2x2S(const doublereal a[3]) {
    return a[0] * a[2] - a[1] * a[1];
}

/*
 * Compute the square root of the symmetric matrix A, B = sqrtm(A)
 *
 * a[0] = A(1,1)
 * a[1] = A(1,2)
 * a[2] = A(2,2)
 *
 * b[0] = B(1,1)
 * b[1] = B(1,2)
 * b[2] = B(2,2)
 *
 * Formulae taken from https://en.wikipedia.org/wiki/Square_root_of_a_2_by_2_matrix
 *
 */
void sqrtm2x2S(const doublereal a[3], doublereal b[3]) {

    doublereal det = det2x2S(a);
    doublereal trace = a[0] + a[2];
    doublereal s = std::sqrt(det);
    doublereal t = std::sqrt(trace + 2. * s);
    b[0] = (a[0] + s) / t;
    b[1] = a[1] / t;
    b[2] = (a[2] + s) / t;


    // std::cout << "(1,1): " << b[0]*b[0] + b[1]*b[1] - a[0] << std::endl;
    // std::cout << "(1,2): " << b[0]*b[1] + b[1]*b[2] - a[1] << std::endl;
    // std::cout << "(2,2): " << b[1]*b[1] + b[2]*b[2] - a[2] << std::endl;
    return;
}

/*
 * Compute the inverse of the symmetric matrix A, B = inv(A)
 *
 * a[0] = A(1,1)
 * a[1] = A(1,2)
 * a[2] = A(2,2)
 *
 * b[0] = B(1,1)
 * b[1] = B(1,2)
 * b[2] = B(2,2)
 *
 */
void inv2x2S(const doublereal a[3], doublereal b[3]) {
    doublereal det = a[0]*a[2] - a[1]*a[1]; //t2+t4;
    doublereal inv_det = 1.0 / det;

    b[0] = a[2] * inv_det;
    b[1] = -a[1] * inv_det; // t8;
    b[2] = a[0] * inv_det;

    // std::cout << "(1,1): " << b[0]*a[0] + b[1]*a[1] << std::endl;
    // std::cout << "(1,2): " << b[0]*a[1] + b[1]*a[2] << std::endl;
    // std::cout << "(1,2): " << b[1]*a[1] + b[2]*a[2] << std::endl;
    return;
}

/*
 * Sort in place v in ascending order wrt |v|, and return the sorting indices
 *
 * Equivalent to Matlab's [v, idx] = sort(v)
 */
std::vector<integer> sort_vector(Vec3& v) {
    std::vector<integer> idx({1, 2, 3});
    std::stable_sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2) {return std::abs(v(i1)) < std::abs(v(i2));});
    Vec3 r;
    r(1) = v(idx[0]);
    r(2) = v(idx[1]);
    r(3) = v(idx[2]);
    v = r;
    return idx;
}


/*
 * Given the vector r, with r.Norm() > 0, compute an orthogonal matrix Q with
 * Q(1,:) aligned with r
 * Q(2,:) and Q(3,:) tangent to the unit sphere
 *
 * if(update) then minimize the rotation around r in ordert to match as mush a possible Qold
 */
void SphericalQR(const Vec3 & r, Mat3x3 &Q, const bool update = false, const Mat3x3& Qold = Eye3) {
    // std::cout << "QDOld: " << Qold << std::endl;
    // std::cout << "QD: " << Q << std::endl;
    // std::cout << "r: " << r << std::endl;

    Vec3 q1 = r / r.Norm();
    Vec3 q2 = q1;
    // std::cout << "q2: " << q2 << std::endl;
    std::vector<integer> sort_idx = sort_vector(q2);
    // std::cout << "sort_idx: " << sort_idx[0] << " " << sort_idx[1] << " " << sort_idx[2] << std::endl;

    // std::cout << q1 << "\n\n" << q2 << "\n\n" << sort_idx[0] <<
    //     " " << sort_idx[1] << " " << sort_idx[2] << std::endl;

    Vec3 v2;
    v2(sort_idx[2]) = q2(2);
    v2(sort_idx[1]) = -q2(3);
    v2(sort_idx[0]) = 0.;
    // std::cout << "v2: " << v2 << std::endl;
    v2 = v2 / v2.Norm();
    // std::cout << q1.Dot(v2) << std::endl;
    Vec3 v3 = q1.Cross(v2);
    Q = Mat3x3(q1, v2, v3);
    // std::cout << "\nortocheck: " << Q.MulMT(Q) << std::endl;

    // std::cout << "update: " << update << std::endl;
    if (update ) {
        Vec3 n1 = Qold.GetCol(1); n1 = n1 / n1.Norm();
        Vec3 n2 = q1;
        doublereal costheta = n1.Dot(n2);
        if (costheta >= 1. - std::numeric_limits<doublereal>::epsilon() * 10.) {
            Q = Qold;
        } else if (costheta <= -1. + std::numeric_limits<doublereal>::epsilon() * 10.) {
            Mat3x3 R = RotManip::Rot(Qold.GetCol(2) * std::numbers::pi);
            Q = R * Qold;
        } else {
            Vec3 k = n1.Cross(n2);
            doublereal theta = std::acos(costheta);
            k = k / k.Norm();
            Mat3x3 R = RotManip::Rot(k * theta);
            Q = R * Qold;
        }
        // std::cerr << "q1: " << Q.GetCol(1) << std::endl;
        // std::cerr << "Q1: " << Q.GetCol(1) << std::endl;
        // std::cerr << "Q2: " << Q.GetCol(2) << std::endl;
        // std::cerr << "Q3: " << Q.GetCol(3) << std::endl;
    }
}

// void SphericalQR(const Vec3 & r, Mat3x3 &Q, const bool update = false, const Mat3x3& Qold = Eye3) {
//     // std::cout << "QDOld: " << Qold << std::endl;
//     // std::cout << "QD: " << Q << std::endl;
//     // std::cout << "r: " << r << std::endl;
//
//     Vec3 q1 = r / r.Norm();
//     Vec3 q2 = q1;
//     // std::cout << "q2: " << q2 << std::endl;
//     std::vector<integer> sort_idx = sort_vector(q2);
//     // std::cout << "sort_idx: " << sort_idx[0] << " " << sort_idx[1] << " " << sort_idx[2] << std::endl;
//
//     // std::cout << q1 << "\n\n" << q2 << "\n\n" << sort_idx[0] <<
//     //     " " << sort_idx[1] << " " << sort_idx[2] << std::endl;
//
//     Vec3 v2;
//     v2(sort_idx[2]) = q2(2);
//     v2(sort_idx[1]) = -q2(3);
//     v2(sort_idx[0]) = 0.;
//     // std::cout << "v2: " << v2 << std::endl;
//     v2 = v2 / v2.Norm();
//     // std::cout << q1.Dot(v2) << std::endl;
//     Vec3 v3 = q1.Cross(v2);
//     Q = Mat3x3(q1, v2, v3);
//     // std::cout << "\nortocheck: " << Q.MulMT(Q) << std::endl;
//
//     // std::cout << "update: " << update << std::endl;
//     if (update ) {
//         // doublereal phi = RotManip::VecRot(Qold.MulTM(Q)).Norm();
//         // if (phi < 0.785398163397448) {
//
//             // std::cout << "old: " << Qold << std::endl;
//             doublereal c[2][2], cct[3], sqrtc[3], isqrtc[3], u[2][2];
//             // std::cout << "c" << std::endl;
//             for (int i = 0; i < 2; i++) {
//                 for (int j = 0; j < 2; j++) {
//                     c[i][j] = Qold.GetCol(i+2).Dot(Q.GetCol(j+2));
//                     // std::cout << c[i][j] << " ";
//                 }
//                 // std::cout << std::endl;
//             }
//             cct[0] = c[0][0] * c[0][0] + c[0][1] * c[0][1];
//             cct[1] = c[0][0] * c[1][0] + c[0][1] * c[1][1];
//             cct[2] = c[1][0] * c[1][0] + c[1][1] * c[1][1];
//             if (std::abs(det2x2S(cct)) > 0.) {
//                 sqrtm2x2S(cct, sqrtc);
//                 inv2x2S(sqrtc, isqrtc);
//                 u[0][0] = c[0][0] * isqrtc[0] + c[1][0] * isqrtc[1];
//                 u[1][0] = c[0][0] * isqrtc[1] + c[1][0] * isqrtc[2];
//                 u[0][1] = c[0][1] * isqrtc[0] + c[1][1] * isqrtc[1];
//                 u[1][1] = c[0][1] * isqrtc[1] + c[1][1] * isqrtc[2];
//                 // std::cout << "old: " << Qold << std::endl;
//                 // std::cout << "old1: " << Qold.GetVec(1) << std::endl;
//                 // std::cout << "old2: " << Qold.GetVec(2) << std::endl;
//                 // std::cout << "old3: " << Qold.GetVec(3) << std::endl;
//                 // std::cout << "bef: " << Q << std::endl;
//                 // std::cout << "bef1: " << Q.GetVec(1) << std::endl;
//                 // std::cout << "bef2: " << Q.GetVec(2) << std::endl;
//                 // std::cout << "bef3: " << Q.GetVec(3) << std::endl;
//                 Vec3 q1 = Q.GetCol(2) * u[0][0] + Q.GetCol(3) * u[1][0];
//                 Vec3 q2 = Q.GetCol(2) * u[0][1] + Q.GetCol(3) * u[1][1];
//                 Q.PutVec(2, q1);
//                 Q.PutVec(3, q2);
//                 // std::cout << "aft: " << Q << std::endl;
//                 // std::cout << "aft1: " << Q.GetVec(1) << std::endl;
//                 // std::cout << "aft2: " << Q.GetVec(2) << std::endl;
//                 // std::cout << "aft3: " << Q.GetVec(3) << std::endl;
//                 // std::cout << "\nortocheck: " << Q.MulMT(Q) << std::endl;
//                 // doublereal phi = std::abs(RotManip::VecRot(Qold.MulTM(Q)).Dot(Q.GetVec(1)));
//                 doublereal phi = RotManip::VecRot(Qold.MulTM(Q)).Norm();
//                 // std::cerr << "phi: " << RotManip::VecRot(Qold.MulTM(Q)) << std::endl
//                 //     << "\t" << RotManip::VecRot(Qold.MulTM(Q)).Dot(Q.GetVec(1)) << std::endl
//                 //     << "\t" << std::abs(RotManip::VecRot(Qold.MulTM(Q)).Dot(Q.GetVec(1))) << std::endl
//                 //     << "\t" << RotManip::VecRot(Qold.MulTM(Q)).Norm() << std::endl;
//                 if (phi > std::numbers::pi / 2.) {
//                     // std::cout << "phi before " << phi << std::endl;
//                     // Vec3 p = Q.GetVec(1) * std::numbers::pi;
//                     Vec3 p = Vec3(1., 0., 0.) * std::numbers::pi;
//                     Mat3x3 R = RotManip::Rot(p);
//                     Q = Q * R;
//                     // std::cout << "aftaft: " << Q << std::endl;
//                     // std::cout << "aftaft1: " << Q.GetVec(1) << std::endl;
//                     // std::cout << "aftaft2: " << Q.GetVec(2) << std::endl;
//                     // std::cout << "aftaft3: " << Q.GetVec(3) << std::endl;
//                     // std::cout << "\nortocheck: " << Q.MulMT(Q) << std::endl;
//                     phi = std::abs(RotManip::VecRot(Qold.MulTM(Q)).Dot(Q.GetVec(1)));
//                     // std::cout << "phi after " << phi << std::endl;
//                 }
//             } else {
//                 // std::cout << "________________________________" << std::endl;
//                 // std::cout << det2x2S(cct) << std::endl;
//                 // std::cout << "________________________________" << std::endl;
//             }
//         // }
//     }
//     // std::cout << "QD after: " << Q << std::endl;
// }

// int main(void) {
//     Vec3 r(878., 234., 1123.); r = r / r.Norm();
//     Mat3x3 Q;
//     SphericalQR(r, Q);
//     Mat3x3 Qold = Q;
//     r = r + Vec3(0.08, 0.0158, 0.0038);
//     SphericalQR(r, Q, true, Qold);
//
//     return 0;
// }
