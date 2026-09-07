/*
 * test_minmax_robust.c — comprehensive test suite for minmax_robust.h
 *
 * Passes under all four combinations:
 *
 *   gcc   -std=c11 -Wall -Werror -O2                    test_minmax_robust.c -lm
 *   clang -std=c11 -Wall -Werror -O2                    test_minmax_robust.c -lm
 *   gcc   -std=c11 -Wall -Werror -O3 -ffast-math -march=native  test_minmax_robust.c -lm
 *   clang -std=c11 -Wall -Werror -O3 -ffast-math -march=native  test_minmax_robust.c -lm
 *
 * (-Ofast is equivalent to -O3 -ffast-math but is deprecated in recent
 * Clang releases; -O3 -ffast-math is portable across GCC and Clang.)
 *
 * NaN-propagation tests are conditionally compiled out under -ffast-math:
 * that flag allows the compiler to replace fminf/fmin/fminl with plain
 * comparisons that do not follow IEEE 754 NaN semantics.  This is an
 * intentional trade-off of -ffast-math, not a bug in the macros.  Those
 * tests are counted as "skipped" rather than failed.
 *
 * __FAST_MATH__ is defined by both GCC and Clang whenever -ffast-math (or
 * any flag that implies it) is active.
 */

#include "f2c-int.h"

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <float.h>
#include <math.h>

/* ═══════════════════════════════════════════════════════════════════════════
 * Test framework
 * ═══════════════════════════════════════════════════════════════════════════ */

static int g_pass, g_fail, g_skip;

/*
 * CHECK — integer / exact-equality check.
 * Uses __typeof__ (GCC/Clang extension) so both sides are compared at
 * their natural type width, with no implicit promotion.
 */
#define CHECK(desc, got, expected)                                          \
    do {                                                                    \
        __typeof__(expected) _g = (got);                                    \
        __typeof__(expected) _e = (expected);                               \
        if (_g == _e) {                                                     \
            g_pass++;                                                       \
        } else {                                                            \
            g_fail++;                                                       \
            printf("FAIL  %s\n      got 0x%llx  expected 0x%llx\n",        \
                   (desc),                                                  \
                   (unsigned long long)_g,                                  \
                   (unsigned long long)_e);                                 \
        }                                                                   \
    } while (0)

/*
 * CHECK_F — floating-point exact-equality check.
 *
 * Uses the self-comparison NaN test (_x != _x) rather than isnan():
 * isnan() may be folded to a constant 0 under -ffinite-math-only, but
 * the compiler cannot legally eliminate a comparison of a *runtime*
 * variable with itself (the value might have been loaded from memory
 * after the guard below was evaluated).
 *
 * NaN test cases should be placed inside #ifndef __FAST_MATH__ sections
 * because under fast-math the macros may themselves return NaN.
 */
#define CHECK_F(desc, got, expected)                                        \
    do {                                                                    \
        __typeof__(expected) _g = (got);                                    \
        __typeof__(expected) _e = (expected);                               \
        int _exp_nan = (_e != _e);          /* true iff expected is NaN */  \
        int _ok = _exp_nan ? (_g != _g) : (_g == _e);                      \
        if (_ok) {                                                          \
            g_pass++;                                                       \
        } else {                                                            \
            g_fail++;                                                       \
            printf("FAIL  %s\n      got %Lg  expected %Lg\n",              \
                   (desc), (long double)_g, (long double)_e);              \
        }                                                                   \
    } while (0)

/*
 * CHECK_TYPE — verify that an expression has the same size as a given type.
 * sizeof is a compile-time operator; no runtime overhead.
 * Must be invoked as  CHECK_TYPE("desc", expr, (expected_type)0)  because
 * a macro argument cannot be a multi-token type name like "long long" in
 * standard C (no variadic tricks needed when we pass a cast-to-zero instead).
 */
#define CHECK_TYPE(desc, expr, zero_of_type)                                \
    do {                                                                    \
        if (sizeof(expr) == sizeof(zero_of_type)) {                         \
            g_pass++;                                                       \
        } else {                                                            \
            g_fail++;                                                       \
            printf("FAIL  %s: size %zu, expected %zu\n",                   \
                   (desc), sizeof(expr), sizeof(zero_of_type));            \
        }                                                                   \
    } while (0)

/* Mark tests as intentionally skipped. */
#define SKIP(desc, reason)                                                  \
    do { g_skip++; printf("skip  %-48s (%s)\n", (desc), (reason)); }       \
    while (0)

static void section(const char *name) { printf("\n── %s\n", name); }

/* ═══════════════════════════════════════════════════════════════════════════
 * Side-effect helpers — each call increments g_se_count exactly once
 * ═══════════════════════════════════════════════════════════════════════════ */

static int g_se_count;

static int            se_int   (int v)            { g_se_count++; return v; }
static long           se_long  (long v)           { g_se_count++; return v; }
static long long      se_llong (long long v)      { g_se_count++; return v; }
static double         se_double(double v)         { g_se_count++; return v; }
static unsigned int   se_uint  (unsigned int v)   { g_se_count++; return v; }
static unsigned long long se_ullong(unsigned long long v)
                                                  { g_se_count++; return v; }

/* ═══════════════════════════════════════════════════════════════════════════
 * signed int
 * ═══════════════════════════════════════════════════════════════════════════ */

static void test_int(void) {
    section("signed int");

    CHECK("min(3,5)",    min(3,  5),   3);
    CHECK("min(5,3)",    min(5,  3),   3);
    CHECK("min(-1,0)",   min(-1, 0),  -1);
    CHECK("min(0,-1)",   min(0, -1),  -1);
    CHECK("min(-5,-3)",  min(-5,-3),  -5);
    CHECK("min(-3,-5)",  min(-3,-5),  -5);
    CHECK("min(x,x)",    min(7,  7),   7);

    CHECK("max(3,5)",    max(3,  5),   5);
    CHECK("max(5,3)",    max(5,  3),   5);
    CHECK("max(-1,0)",   max(-1, 0),   0);
    CHECK("max(0,-1)",   max(0, -1),   0);
    CHECK("max(-5,-3)",  max(-5,-3),  -3);
    CHECK("max(-3,-5)",  max(-3,-5),  -3);
    CHECK("max(x,x)",    max(7,  7),   7);

    /* extreme values — triggered signed-overflow UB in the original code */
    CHECK("min(INT_MAX,INT_MIN)", min(INT_MAX, INT_MIN), INT_MIN);
    CHECK("min(INT_MIN,INT_MAX)", min(INT_MIN, INT_MAX), INT_MIN);
    CHECK("min(INT_MAX,0)",       min(INT_MAX, 0),       0);
    CHECK("min(0,INT_MIN)",       min(0, INT_MIN),       INT_MIN);
    CHECK("min(INT_MAX,INT_MAX)", min(INT_MAX, INT_MAX), INT_MAX);
    CHECK("min(INT_MIN,INT_MIN)", min(INT_MIN, INT_MIN), INT_MIN);

    CHECK("max(INT_MAX,INT_MIN)", max(INT_MAX, INT_MIN), INT_MAX);
    CHECK("max(INT_MIN,INT_MAX)", max(INT_MIN, INT_MAX), INT_MAX);
    CHECK("max(INT_MAX,0)",       max(INT_MAX, 0),       INT_MAX);
    CHECK("max(0,INT_MIN)",       max(0, INT_MIN),       0);
    CHECK("max(INT_MAX,INT_MAX)", max(INT_MAX, INT_MAX), INT_MAX);
    CHECK("max(INT_MIN,INT_MIN)", max(INT_MIN, INT_MIN), INT_MIN);

    CHECK("min(INT_MAX,INT_MAX-1)", min(INT_MAX, INT_MAX-1), INT_MAX-1);
    CHECK("max(INT_MIN,INT_MIN+1)", max(INT_MIN, INT_MIN+1), INT_MIN+1);
}

/* ═══════════════════════════════════════════════════════════════════════════
 * signed long
 * ═══════════════════════════════════════════════════════════════════════════ */

static void test_long(void) {
    section("signed long");

    CHECK("min(3L,5L)",   min(3L, 5L),   3L);
    CHECK("min(-1L,0L)",  min(-1L, 0L), -1L);
    CHECK("max(3L,5L)",   max(3L, 5L),   5L);
    CHECK("max(-1L,0L)",  max(-1L, 0L),  0L);

    CHECK("min(LONG_MAX,LONG_MIN)", min(LONG_MAX, LONG_MIN), LONG_MIN);
    CHECK("min(LONG_MIN,LONG_MAX)", min(LONG_MIN, LONG_MAX), LONG_MIN);
    CHECK("max(LONG_MAX,LONG_MIN)", max(LONG_MAX, LONG_MIN), LONG_MAX);
    CHECK("max(LONG_MIN,LONG_MAX)", max(LONG_MIN, LONG_MAX), LONG_MAX);
    CHECK("min(LONG_MAX,LONG_MAX)", min(LONG_MAX, LONG_MAX), LONG_MAX);
    CHECK("max(LONG_MIN,LONG_MIN)", max(LONG_MIN, LONG_MIN), LONG_MIN);
}

/* ═══════════════════════════════════════════════════════════════════════════
 * signed long long  (was silently truncated to int in original)
 * ═══════════════════════════════════════════════════════════════════════════ */

static void test_llong(void) {
    section("signed long long  (was silently truncated to int in original)");

    CHECK("min(3LL,5LL)",  min(3LL, 5LL),  3LL);
    CHECK("min(-1LL,0LL)", min(-1LL,0LL), -1LL);
    CHECK("max(3LL,5LL)",  max(3LL, 5LL),  5LL);
    CHECK("max(-1LL,0LL)", max(-1LL,0LL),  0LL);

    /* proved truncation bug in original: 3000000000 > INT_MAX */
    CHECK("min(3000000000LL,1LL)",     min(3000000000LL, 1LL),      1LL);
    CHECK("max(1LL,3000000000LL)",     max(1LL, 3000000000LL),      3000000000LL);
    CHECK("min(-3000000000LL,-1LL)",   min(-3000000000LL, -1LL),   -3000000000LL);
    CHECK("max(-3000000000LL,-1LL)",   max(-3000000000LL, -1LL),   -1LL);

    CHECK("min(LLONG_MAX,LLONG_MIN)", min(LLONG_MAX, LLONG_MIN), LLONG_MIN);
    CHECK("min(LLONG_MIN,LLONG_MAX)", min(LLONG_MIN, LLONG_MAX), LLONG_MIN);
    CHECK("max(LLONG_MAX,LLONG_MIN)", max(LLONG_MAX, LLONG_MIN), LLONG_MAX);
    CHECK("max(LLONG_MIN,LLONG_MAX)", max(LLONG_MIN, LLONG_MAX), LLONG_MAX);
    CHECK("min(LLONG_MAX,LLONG_MAX)", min(LLONG_MAX, LLONG_MAX), LLONG_MAX);
    CHECK("max(LLONG_MIN,LLONG_MIN)", max(LLONG_MIN, LLONG_MIN), LLONG_MIN);
    CHECK("min(LLONG_MIN,LLONG_MIN+1)", min(LLONG_MIN, LLONG_MIN+1), LLONG_MIN);
    CHECK("max(LLONG_MAX-1,LLONG_MAX)", max(LLONG_MAX-1, LLONG_MAX), LLONG_MAX);
}

/* ═══════════════════════════════════════════════════════════════════════════
 * unsigned int  (was routed to signed __min_int in original)
 * ═══════════════════════════════════════════════════════════════════════════ */

static void test_uint(void) {
    section("unsigned int  (was routed to signed __min_int in original)");

    CHECK("min(1u,3u)",  min(1u, 3u),  1u);
    CHECK("min(3u,1u)",  min(3u, 1u),  1u);
    CHECK("min(0u,0u)",  min(0u, 0u),  0u);
    CHECK("max(1u,3u)",  max(1u, 3u),  3u);
    CHECK("max(3u,1u)",  max(3u, 1u),  3u);

    CHECK("min(UINT_MAX,0u)",       min(UINT_MAX, 0u),       0u);
    CHECK("min(0u,UINT_MAX)",       min(0u, UINT_MAX),       0u);
    CHECK("max(UINT_MAX,0u)",       max(UINT_MAX, 0u),       UINT_MAX);
    CHECK("max(0u,UINT_MAX)",       max(0u, UINT_MAX),       UINT_MAX);
    CHECK("min(UINT_MAX,UINT_MAX)", min(UINT_MAX, UINT_MAX), UINT_MAX);

    CHECK("min(UINT_MAX,UINT_MAX-1)", min(UINT_MAX, UINT_MAX-1u), UINT_MAX-1u);
    CHECK("max(UINT_MAX-1,UINT_MAX)", max(UINT_MAX-1u, UINT_MAX), UINT_MAX);
}

/* ═══════════════════════════════════════════════════════════════════════════
 * unsigned long
 * ═══════════════════════════════════════════════════════════════════════════ */

static void test_ulong(void) {
    section("unsigned long");

    CHECK("min(1ul,3ul)",            min(1ul, 3ul),            1ul);
    CHECK("max(1ul,3ul)",            max(1ul, 3ul),            3ul);
    CHECK("min(ULONG_MAX,0ul)",      min(ULONG_MAX, 0ul),      0ul);
    CHECK("max(0ul,ULONG_MAX)",      max(0ul, ULONG_MAX),      ULONG_MAX);
    CHECK("min(ULONG_MAX,ULONG_MAX)",min(ULONG_MAX, ULONG_MAX),ULONG_MAX);
    CHECK("min(ULONG_MAX,ULONG_MAX-1)", min(ULONG_MAX, ULONG_MAX-1ul), ULONG_MAX-1ul);
    CHECK("max(ULONG_MAX-1,ULONG_MAX)", max(ULONG_MAX-1ul, ULONG_MAX), ULONG_MAX);
}

/* ═══════════════════════════════════════════════════════════════════════════
 * unsigned long long  (was silently truncated to int in original)
 * ═══════════════════════════════════════════════════════════════════════════ */

static void test_ullong(void) {
    section("unsigned long long  (was silently truncated to int in original)");

    CHECK("min(1ull,3ull)",           min(1ull, 3ull),           1ull);
    CHECK("max(1ull,3ull)",           max(1ull, 3ull),           3ull);
    CHECK("min(ULLONG_MAX,0ull)",     min(ULLONG_MAX, 0ull),     0ull);
    CHECK("min(0ull,ULLONG_MAX)",     min(0ull, ULLONG_MAX),     0ull);
    CHECK("max(ULLONG_MAX,0ull)",     max(ULLONG_MAX, 0ull),     ULLONG_MAX);
    CHECK("max(0ull,ULLONG_MAX)",     max(0ull, ULLONG_MAX),     ULLONG_MAX);
    CHECK("min(ULLONG_MAX,ULLONG_MAX)",   min(ULLONG_MAX, ULLONG_MAX),     ULLONG_MAX);
    CHECK("min(ULLONG_MAX,ULLONG_MAX-1)", min(ULLONG_MAX, ULLONG_MAX-1ull),ULLONG_MAX-1ull);
    CHECK("max(ULLONG_MAX-1,ULLONG_MAX)", max(ULLONG_MAX-1ull, ULLONG_MAX),ULLONG_MAX);
}

/* ═══════════════════════════════════════════════════════════════════════════
 * float
 * ═══════════════════════════════════════════════════════════════════════════ */

static void test_float(void) {
    section("float");

    CHECK_F("min(1.5f,2.5f)",         min(1.5f,  2.5f),    1.5f);
    CHECK_F("min(2.5f,1.5f)",         min(2.5f,  1.5f),    1.5f);
    CHECK_F("min(-1.5f,0.0f)",        min(-1.5f, 0.0f),   -1.5f);
    CHECK_F("min(0.0f,0.0f)",         min(0.0f,  0.0f),    0.0f);
    CHECK_F("max(1.5f,2.5f)",         max(1.5f,  2.5f),    2.5f);
    CHECK_F("max(-1.5f,0.0f)",        max(-1.5f, 0.0f),    0.0f);
    CHECK_F("min(FLT_MAX,-FLT_MAX)",  min(FLT_MAX,-FLT_MAX), -FLT_MAX);
    CHECK_F("max(FLT_MAX,-FLT_MAX)",  max(FLT_MAX,-FLT_MAX),  FLT_MAX);
    CHECK_F("min(FLT_MIN,FLT_MAX)",   min(FLT_MIN, FLT_MAX),  FLT_MIN);

    /*
     * IEEE 754 NaN semantics: fminf(NaN, x) returns x (not NaN).
     *
     * Skipped under -ffast-math: the compiler is permitted to replace
     * fminf with a plain comparison, which propagates NaN.  This is a
     * documented trade-off of -ffast-math, not a bug in the macros.
     */
#ifndef __FAST_MATH__
    {
        float nan_f = (float)NAN;
        CHECK_F("min(NaN,1.0f) == 1.0f",  min(nan_f, 1.0f),  1.0f);
        CHECK_F("max(NaN,1.0f) == 1.0f",  max(nan_f, 1.0f),  1.0f);
        CHECK_F("min(1.0f,NaN) == 1.0f",  min(1.0f, nan_f),  1.0f);
        CHECK_F("max(1.0f,NaN) == 1.0f",  max(1.0f, nan_f),  1.0f);
    }
#else
    SKIP("min(NaN,1.0f) == 1.0f", "-ffast-math: fminf NaN semantics undefined");
    SKIP("max(NaN,1.0f) == 1.0f", "-ffast-math: fminf NaN semantics undefined");
    SKIP("min(1.0f,NaN) == 1.0f", "-ffast-math: fminf NaN semantics undefined");
    SKIP("max(1.0f,NaN) == 1.0f", "-ffast-math: fminf NaN semantics undefined");
#endif
}

/* ═══════════════════════════════════════════════════════════════════════════
 * double
 * ═══════════════════════════════════════════════════════════════════════════ */

static void test_double(void) {
    section("double");

    CHECK_F("min(1.5,2.5)",           min(1.5,  2.5),   1.5);
    CHECK_F("min(-1.0,0.0)",          min(-1.0, 0.0),  -1.0);
    CHECK_F("max(1.5,2.5)",           max(1.5,  2.5),   2.5);
    CHECK_F("max(-1.0,0.0)",          max(-1.0, 0.0),   0.0);
    CHECK_F("min(DBL_MAX,-DBL_MAX)",  min(DBL_MAX,-DBL_MAX), -DBL_MAX);
    CHECK_F("max(DBL_MAX,-DBL_MAX)",  max(DBL_MAX,-DBL_MAX),  DBL_MAX);
    CHECK_F("min(DBL_MIN,DBL_MAX)",   min(DBL_MIN, DBL_MAX),  DBL_MIN);

#ifndef __FAST_MATH__
    {
        double nan_d = (double)NAN;
        CHECK_F("min(NaN,1.0) == 1.0",  min(nan_d, 1.0),  1.0);
        CHECK_F("max(NaN,1.0) == 1.0",  max(nan_d, 1.0),  1.0);
        CHECK_F("min(1.0,NaN) == 1.0",  min(1.0, nan_d),  1.0);
        CHECK_F("max(1.0,NaN) == 1.0",  max(1.0, nan_d),  1.0);
    }
#else
    SKIP("min(NaN,1.0) == 1.0",  "-ffast-math: fmin NaN semantics undefined");
    SKIP("max(NaN,1.0) == 1.0",  "-ffast-math: fmin NaN semantics undefined");
    SKIP("min(1.0,NaN) == 1.0",  "-ffast-math: fmin NaN semantics undefined");
    SKIP("max(1.0,NaN) == 1.0",  "-ffast-math: fmin NaN semantics undefined");
#endif
}

/* ═══════════════════════════════════════════════════════════════════════════
 * long double
 * ═══════════════════════════════════════════════════════════════════════════ */

static void test_ldouble(void) {
    section("long double");

    CHECK_F("min(1.5L,2.5L)",          min(1.5L,  2.5L),    1.5L);
    CHECK_F("max(1.5L,2.5L)",          max(1.5L,  2.5L),    2.5L);
    CHECK_F("min(LDBL_MAX,-LDBL_MAX)", min(LDBL_MAX,-LDBL_MAX), -LDBL_MAX);
    CHECK_F("max(LDBL_MAX,-LDBL_MAX)", max(LDBL_MAX,-LDBL_MAX),  LDBL_MAX);

#ifndef __FAST_MATH__
    {
        long double nan_ld = (long double)NAN;
        CHECK_F("min(NaN,1.0L) == 1.0L", min(nan_ld, 1.0L), 1.0L);
        CHECK_F("max(NaN,1.0L) == 1.0L", max(nan_ld, 1.0L), 1.0L);
        CHECK_F("min(1.0L,NaN) == 1.0L", min(1.0L, nan_ld), 1.0L);
        CHECK_F("max(1.0L,NaN) == 1.0L", max(1.0L, nan_ld), 1.0L);
    }
#else
    SKIP("min(NaN,1.0L) == 1.0L", "-ffast-math: fminl NaN semantics undefined");
    SKIP("max(NaN,1.0L) == 1.0L", "-ffast-math: fminl NaN semantics undefined");
    SKIP("min(1.0L,NaN) == 1.0L", "-ffast-math: fminl NaN semantics undefined");
    SKIP("max(1.0L,NaN) == 1.0L", "-ffast-math: fminl NaN semantics undefined");
#endif
}

/* ═══════════════════════════════════════════════════════════════════════════
 * Mixed-type promotion via (a)+(b) in _Generic
 * ═══════════════════════════════════════════════════════════════════════════ */

static void test_mixed(void) {
    section("mixed-type promotion");

    int       i3     = 3;
    long      l5     = 5L;
    long long ll_big = 3000000000LL;
    long long ll2    = 2LL;
    unsigned int      u5      = 5u;
    unsigned long long ull10  = 10ull;

    /* int + long → long */
    CHECK("min(int 3,long 5) → 3L",  min(i3, l5), 3L);
    CHECK("max(int 3,long 5) → 5L",  max(i3, l5), 5L);
    CHECK_TYPE("min(int,long) → long",  min(i3, l5), (long)0);

    /* int + long long → long long */
    CHECK("min(int 3,ll big) → 3LL", min(i3, ll_big), (long long)3);
    CHECK("max(int 3,ll big) → big", max(i3, ll_big), ll_big);
    CHECK_TYPE("min(int,ll) → long long", min(i3, ll_big), (long long)0);

    /* int + unsigned int → unsigned int */
    CHECK("min(int 3,uint 5) → 3u",  min(i3, u5), 3u);
    CHECK_TYPE("min(int,uint) → uint", min(i3, u5), (unsigned int)0);

    /* int + float → float */
    CHECK_F("min(int 2,float 1.5f)",  min(2, 1.5f), 1.5f);
    CHECK_TYPE("min(int,float) → float", min(2, 1.5f), (float)0);

    /* int + double → double */
    CHECK_F("min(int 2,double 1.5)",  min(2, 1.5), 1.5);
    CHECK_TYPE("min(int,double) → double", min(2, 1.5), (double)0);

    /* long + long long → long long */
    CHECK("min(long 5L,ll 2LL) → 2", min(l5, ll2), 2LL);
    CHECK_TYPE("min(long,ll) → long long", min(l5, ll2), (long long)0);

    /* unsigned long long + long long → unsigned long long */
    CHECK_TYPE("min(ull,ll) → ull",  min(ull10, ll2), (unsigned long long)0);
}

/* ═══════════════════════════════════════════════════════════════════════════
 * Side-effect safety — each argument evaluated exactly once
 * ═══════════════════════════════════════════════════════════════════════════ */

static void test_side_effects(void) {
    section("side-effect safety");

    g_se_count = 0;
    (void)min(se_int(3),    se_int(5));
    CHECK("int:    2 evaluations", g_se_count, 2);

    g_se_count = 0;
    (void)min(se_long(3L),  se_long(5L));
    CHECK("long:   2 evaluations", g_se_count, 2);

    g_se_count = 0;
    (void)min(se_llong(3LL), se_llong(5LL));
    CHECK("llong:  2 evaluations", g_se_count, 2);

    g_se_count = 0;
    (void)min(se_double(3.0), se_double(5.0));
    CHECK("double: 2 evaluations", g_se_count, 2);

    g_se_count = 0;
    (void)max(se_uint(3u),   se_uint(5u));
    CHECK("uint:   2 evaluations", g_se_count, 2);

    g_se_count = 0;
    (void)max(se_ullong(3ull), se_ullong(ULLONG_MAX));
    CHECK("ullong: 2 evaluations", g_se_count, 2);

    /*
     * Classic double-evaluation trap from  #define min(a,b) ((a)<(b)?(a):(b)):
     * one of x++ or y++ would be incremented twice.  Verify that cannot
     * happen with our _Generic-based macros.
     */
    int x = 3, y = 5;
    (void)min(x++, y++);
    CHECK("min(x++,y++): x incremented once", x, 4);
    CHECK("min(x++,y++): y incremented once", y, 6);
}

/* ═══════════════════════════════════════════════════════════════════════════
 * Algebraic properties
 * ═══════════════════════════════════════════════════════════════════════════ */

static void test_algebraic(void) {
    section("algebraic properties");

    /* commutativity */
    CHECK("min symmetric int",   min(7, -3),               min(-3, 7));
    CHECK("min symmetric ll",    min(LLONG_MIN, LLONG_MAX), min(LLONG_MAX, LLONG_MIN));
    CHECK("max symmetric uint",  max(0u, UINT_MAX),         max(UINT_MAX, 0u));

    /* idempotency */
    CHECK("min(a,a)==a int",     min(42, 42),     42);
    CHECK("min(a,a)==a ll",      min(-1LL,-1LL), -1LL);
    CHECK("max(a,a)==a uint",    max(7u, 7u),     7u);

    /* min ≤ both operands, max ≥ both operands */
    int a = -100, b = 200;
    int mn = min(a, b), mx = max(a, b);
    CHECK("min(a,b) <= a", mn <= a, 1);
    CHECK("min(a,b) <= b", mn <= b, 1);
    CHECK("max(a,b) >= a", mx >= a, 1);
    CHECK("max(a,b) >= b", mx >= b, 1);

    /* min(a,b) + max(a,b) == a + b (small values, no overflow risk) */
    CHECK("min+max == a+b", mn + mx, a + b);

    /* absorption laws */
    CHECK("min(a,max(a,b))==a",  min(a, max(a, b)), a);
    CHECK("max(a,min(a,b))==a",  max(a, min(a, b)), a);
}

/* ═══════════════════════════════════════════════════════════════════════════
 * main
 * ═══════════════════════════════════════════════════════════════════════════ */

int main(void) {
#ifdef __FAST_MATH__
    printf("note: compiled with -ffast-math; "
           "NaN-propagation tests will be skipped.\n");
#endif

    test_int();
    test_long();
    test_llong();
    test_uint();
    test_ulong();
    test_ullong();
    test_float();
    test_double();
    test_ldouble();
    test_mixed();
    test_side_effects();
    test_algebraic();

    printf("\n%s  —  %d passed, %d failed",
           g_fail ? "FAILED" : "PASSED", g_pass, g_fail);
    if (g_skip)
        printf(", %d skipped", g_skip);
    printf("\n");

    return g_fail ? EXIT_FAILURE : EXIT_SUCCESS;
}
