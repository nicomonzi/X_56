// MIT License
// 
// Copyright (c) 2017 Martin Bisson
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#include "mbconfig.h"

// Inclusion of the actual implementation of the algorithm.
#include "polar_decomposition_3x3.h"
#include "polar_decomposition_3x3_impl.h"


namespace polar
{


    // Compute the polar decomposition of the input matrix.
    //
    // This method decomposes the input matrix into one orthonormal (rotation) matrix
    // and a remaining positive semi-definite (scale) matrix.
    //
    // Matrices are represented as 3x3 matrices in column-major representation.
    //
    // A [in]  : Matrix to decompose.
    // Q [out] : Rotation part of the polar decomposition.
    // H [out] : Symmetric matrix represention the scale part of the polar decomposition.
    void polar_decomposition(Mat3x3& Q, Mat3x3& H, const Mat3x3& A)
    {
        // Convert parameters for the algorithm implementation.
        typedef detail::matrix<doublereal, 3, 3> TMatrix;

        // Depending on the implementation, this part might involve a copy,
        // but here the memory layout requirement is for matrices to be 3x3,
        // column-major, so we can just cast directly.
        TMatrix* matrixQ = reinterpret_cast<TMatrix*>(Q.pGetMat());
        TMatrix* matrixH = reinterpret_cast<TMatrix*>(H.pGetMat());
        const TMatrix* matrixA = reinterpret_cast<const TMatrix*>(A.pGetMat());

        detail::run_algorithm_3_5(*matrixQ, *matrixH, *matrixA);
    }


}; // End of namespace polar.

