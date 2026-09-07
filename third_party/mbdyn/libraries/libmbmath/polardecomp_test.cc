// This is a example test of the polar decomposition algorithm implementation.

#include "mbconfig.h"
#include "polar_decomposition_3x3.h"
#include "polar_decomposition_3x3_impl.h"

#include "matvec3.h"

#include <iostream>


namespace polar
{


    namespace detail
    {


        template < typename doublereal >
        void multiply(
            matrix<doublereal, 3, 3>& result,
            const matrix<doublereal, 3, 3>& a,
            const matrix<doublereal, 3, 3>& b
            )
        {
            result(0,0) = a(0,0) * b(0,0) + a(1,0) * b(0,1) + a(2,0) * b(0,2);
            result(0,1) = a(0,1) * b(0,0) + a(1,1) * b(0,1) + a(2,1) * b(0,2);
            result(0,2) = a(0,2) * b(0,0) + a(1,2) * b(0,1) + a(2,2) * b(0,2);
            result(1,0) = a(0,0) * b(1,0) + a(1,0) * b(1,1) + a(2,0) * b(1,2);
            result(1,1) = a(0,1) * b(1,0) + a(1,1) * b(1,1) + a(2,1) * b(1,2);
            result(1,2) = a(0,2) * b(1,0) + a(1,2) * b(1,1) + a(2,2) * b(1,2);
            result(2,0) = a(0,0) * b(2,0) + a(1,0) * b(2,1) + a(2,0) * b(2,2);
            result(2,1) = a(0,1) * b(2,0) + a(1,1) * b(2,1) + a(2,1) * b(2,2);
            result(2,2) = a(0,2) * b(2,0) + a(1,2) * b(2,1) + a(2,2) * b(2,2);
        }


        template < typename doublereal >
        void subtract(
            matrix<doublereal, 3, 3>& result,
            const matrix<doublereal, 3, 3>& a,
            const matrix<doublereal, 3, 3>& b
            )
        {
            result(0) = a(0) - b(0);
            result(1) = a(1) - b(1);
            result(2) = a(2) - b(2);
            result(3) = a(3) - b(3);
            result(4) = a(4) - b(4);
            result(5) = a(5) - b(5);
            result(6) = a(6) - b(6);
            result(7) = a(7) - b(7);
            result(8) = a(8) - b(8);
        }


        template < typename doublereal >
        inline doublereal norm(const matrix<doublereal, 3, 3>& m)
        {
            doublereal length = m(0) * m(0);
            length += m(1) * m(1);
            length += m(2) * m(2);
            length += m(3) * m(3);
            length += m(4) * m(4);
            length += m(5) * m(5);
            length += m(6) * m(6);
            length += m(7) * m(7);
            length += m(8) * m(8);

            return math_utils<doublereal>::sqrt(length);
        }


        template < typename doublereal >
        inline matrix<doublereal, 3, 3> identity()
        {
            matrix<doublereal, 3, 3> m;
            m(0) = 1;
            m(1) = 0;
            m(2) = 0;
            m(3) = 0;
            m(4) = 1;
            m(5) = 0;
            m(6) = 0;
            m(7) = 0;
            m(8) = 1;
            return m;
        }


    }; // End of namespace detail.


}; // End of namespace polar.




namespace
{


    // These methods validate the output of the algorithm.
    void Check(
        const Mat3x3 valuesA,
        const Mat3x3 valuesQ,
        const Mat3x3 valuesH
        )
    {
        const polar::detail::matrix<doublereal, 3, 3>& A = *reinterpret_cast<const polar::detail::matrix<doublereal, 3, 3>*>(valuesA.pGetMat());
        const polar::detail::matrix<doublereal, 3, 3>& Q = *reinterpret_cast<const polar::detail::matrix<doublereal, 3, 3>*>(valuesQ.pGetMat());
        const polar::detail::matrix<doublereal, 3, 3>& H = *reinterpret_cast<const polar::detail::matrix<doublereal, 3, 3>*>(valuesH.pGetMat());

        polar::detail::matrix<doublereal, 3, 3> temp;
        polar::detail::multiply(temp, Q, H);
        polar::detail::subtract(temp, A, temp);
        const doublereal residual = polar::detail::norm(temp) / (polar::detail::norm(A) + std::numeric_limits<doublereal>::min());

        polar::detail::transpose_multiply(temp, Q, Q);
        polar::detail::subtract(temp, temp, polar::detail::identity<doublereal>());
        const doublereal orthogonality = norm(temp);

        std::cout << "Relative residual = " << residual << " ,  orthogonality = " << orthogonality << std::endl;
        if (residual > 1E-14 || orthogonality > 1E-14) {
            std::cout << "==> failed!" << std::endl;
            MBDYN_TESTSUITE_ASSERT(0);
        }
    }


    template <typename doublereal>
    void CheckExactQ(
        const doublereal y,
        const Mat3x3 valuesQ
        )
    {
        const polar::detail::matrix<doublereal, 3, 3>& Q = *reinterpret_cast<const polar::detail::matrix<doublereal, 3, 3>*>(valuesQ.pGetMat());

        static const Mat3x3 valuesQ1(
            139 / static_cast<doublereal>(255), 466 / static_cast<doublereal>(1275), 962 / static_cast<doublereal>(1275),
            -14 / static_cast<doublereal>(51), -197 / static_cast<doublereal>(255), 146 / static_cast<doublereal>(255),
            202 / static_cast<doublereal>(255), -662 / static_cast<doublereal>(1275), -409 / static_cast<doublereal>(1275)
        );
        const polar::detail::matrix<doublereal, 3, 3>& Q1 = *reinterpret_cast<const polar::detail::matrix<doublereal, 3, 3>*>(valuesQ1.pGetMat());
        // Condition number of U.
        const doublereal kappa = polar::detail::math_utils<doublereal>::sqrt((1 + 2 * y * y) / (3 * y * y));

        polar::detail::matrix<doublereal, 3, 3> temp;
        polar::detail::subtract(temp, Q, Q1);
        const doublereal error = polar::detail::norm(temp) / (polar::detail::norm(Q1) * kappa);
        std::cout << "Condition number of U = " << kappa << ",  scaled relative error in Q = " << error << std::endl;
        if (error > 1E-14) {
            std::cout << "==> failed!" << std::endl;
            MBDYN_TESTSUITE_ASSERT(0);
        }

    }


    MBDYN_TESTSUITE_TEST(itertest, itertest1)
    {
        // This implementation is column-major.
        std::cout << "Test (5.1) from paper:" << std::endl;
        {
            const Mat3x3 A(0.1, 0.1, 0.3, 0.2, 0.1, 0.2, 0.3, 0.0, 0.1);
            Mat3x3 Q;
            Mat3x3 H;
            polar::polar_decomposition(Q, H, A);
            Check(A, Q, H);
        }
        std::cout << "***************************" << std::endl;

        std::cout << "Test (5.2) from paper:" << std::endl;
        static const doublereal yValues[] = {
            static_cast<doublereal>(1),
            static_cast<doublereal>(1.0e-4),
            static_cast<doublereal>(1.0e-8),
            static_cast<doublereal>(1.0e-12),
            static_cast<doublereal>(1.0e-16),
        };
        static const int yValuesCount = sizeof(yValues) / sizeof(yValues[0]);
        for (int i = 0; i < yValuesCount; ++i)
        {
            const doublereal y = yValues[i];
            std::cout << i << " : y = " << y << std::endl;

            const Mat3x3 A(
                (720*y - 25) / 1275, (396*y + 70) / 1275, (972*y - 10) / 1275,
                (-650*y + 300) / 1275, (-145*y - 840) / 1275, (610*y + 120) / 1275,
                (710*y + 300) / 1275, (178*y - 840) / 1275, (-529*y + 120) / 1275
            );
            Mat3x3 Q;
            Mat3x3 H;
            polar::polar_decomposition(Q, H, A);
            Check(A, Q, H);
            CheckExactQ(y, Q);
        }
    }


}; // End of anonymous namespace.




MBDYN_DEFINE_OPERATOR_NEW_DELETE

int main(int argc, char* argv[])
{
     MBDYN_TESTSUITE_INIT(&argc, argv);

     return MBDYN_RUN_ALL_TESTS();
}

