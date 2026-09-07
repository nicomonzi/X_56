/* $Header$ */
/*
 * This library comes with MBDyn (C), a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2023
 *
 * Pierangelo Masarati  <pierangelo.masarati@polimi.it>
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

#ifndef AC_F2C_INT_H
#define AC_F2C_INT_H

#ifdef __cplusplus
#warning "<ac/f2c-int.h> is not meant to be included in C++ code"
extern "C" {
#endif /* __cplusplus */

#include <ac/f2c.h>

#include <limits.h>   /* CHAR_BIT                           */
#include <math.h>     /* fminf/fmin/fminl, fmaxf/fmax/fmaxl */

/* Number of value bits in type T. */
#define __BITS(T)  (sizeof(T) * (unsigned)CHAR_BIT)

/* Number of value bits in type T. */
#define __BITS_FOR_MINMAX(T)  (sizeof(T) * (unsigned)CHAR_BIT)

/*
 * BRANCHLESS INTEGER MIN / MAX — design notes
 * ────────────────────────────────────────────
 * Naive bit-hack  b + ((a-b) & ((a-b) >> 31))  has two defects:
 *   (a) signed subtraction (a-b) is undefined behaviour on overflow;
 *   (b) right-shifting a negative signed value is implementation-defined.
 *
 * Fix: do the subtraction in the *unsigned* domain (well-defined wrap),
 * then recover the true signed comparison bit via the two's-complement
 * overflow detector:
 *
 *   d        = (UT)a − (UT)b               [unsigned wrap, no UB]
 *   overflow = MSB( (ua ^ ub) & (ua ^ d) ) [1 iff subtraction overflowed]
 *   true_lt  = MSB(d) ^ overflow           [1 iff a < b, signed]
 *   mask     = −(T)true_lt                 [0 or all-ones]
 *
 *   min(a,b) = (a & mask) | (b & ~mask)
 *   max(a,b) = (a & ~mask) | (b & mask)
 *
 * All right-shifts are on *unsigned* types → always logical, never
 * implementation-defined.
 *
 * For unsigned types: XOR both operands with SIGN_BIT first, which maps
 * the unsigned total order onto the signed total order; the same formula
 * then applies unchanged.
 *
 * _Generic dispatch: explicit casts in every branch suppress
 * implicit-conversion warnings that arise because all branches are
 * syntactically type-checked even when not selected (C11 §6.5.1.1 — the
 * controlling expression is not evaluated, but each association expression
 * is still parsed and type-checked by the compiler).  The casts in
 * non-selected branches are never evaluated and carry zero runtime cost.
 *
 * ── Floating-point and -ffast-math ──────────────────────────────────────
 * The float/double/long double cases delegate to fminf/fmin/fminl (and
 * their fmax counterparts), which follow IEEE 754 under standard
 * optimisation levels.
 *
 * Under -ffast-math (implied by -Ofast), the compiler is permitted to
 * replace these library calls with plain comparisons, which do NOT preserve
 * IEEE 754 NaN semantics: min(NaN, x) may return NaN instead of x.
 * This is an intentional trade-off of -ffast-math; it cannot be worked
 * around portably without disabling fast-math for the affected translation
 * unit.  All integer helpers are pure unsigned arithmetic and are
 * completely unaffected by -ffast-math.
 */

/* ── signed int ─────────────────────────────────────────────────────────── */
static inline int __min_int(int a, int b) {
    unsigned int ua = (unsigned int)a, ub = (unsigned int)b;
    unsigned int d  = ua - ub;
    unsigned int ov = ((ua ^ ub) & (ua ^ d)) >> (__BITS_FOR_MINMAX(int) - 1);
    int mask = -(int)((d >> (__BITS_FOR_MINMAX(int) - 1)) ^ ov);
    return (a & mask) | (b & ~mask);
}
static inline int __max_int(int a, int b) {
    unsigned int ua = (unsigned int)a, ub = (unsigned int)b;
    unsigned int d  = ua - ub;
    unsigned int ov = ((ua ^ ub) & (ua ^ d)) >> (__BITS_FOR_MINMAX(int) - 1);
    int mask = -(int)((d >> (__BITS_FOR_MINMAX(int) - 1)) ^ ov);
    return (a & ~mask) | (b & mask);
}

/* ── signed long ─────────────────────────────────────────────────────────── */
static inline long __min_long(long a, long b) {
    unsigned long ua = (unsigned long)a, ub = (unsigned long)b;
    unsigned long d  = ua - ub;
    unsigned long ov = ((ua ^ ub) & (ua ^ d)) >> (__BITS_FOR_MINMAX(long) - 1);
    long mask = -(long)((d >> (__BITS_FOR_MINMAX(long) - 1)) ^ ov);
    return (a & mask) | (b & ~mask);
}
static inline long __max_long(long a, long b) {
    unsigned long ua = (unsigned long)a, ub = (unsigned long)b;
    unsigned long d  = ua - ub;
    unsigned long ov = ((ua ^ ub) & (ua ^ d)) >> (__BITS_FOR_MINMAX(long) - 1);
    long mask = -(long)((d >> (__BITS_FOR_MINMAX(long) - 1)) ^ ov);
    return (a & ~mask) | (b & mask);
}

/* ── signed long long ───────────────────────────────────────────────────── */
static inline long long __min_llong(long long a, long long b) {
    unsigned long long ua = (unsigned long long)a, ub = (unsigned long long)b;
    unsigned long long d  = ua - ub;
    unsigned long long ov = ((ua ^ ub) & (ua ^ d)) >> (__BITS_FOR_MINMAX(long long) - 1);
    long long mask = -(long long)((d >> (__BITS_FOR_MINMAX(long long) - 1)) ^ ov);
    return (a & mask) | (b & ~mask);
}
static inline long long __max_llong(long long a, long long b) {
    unsigned long long ua = (unsigned long long)a, ub = (unsigned long long)b;
    unsigned long long d  = ua - ub;
    unsigned long long ov = ((ua ^ ub) & (ua ^ d)) >> (__BITS_FOR_MINMAX(long long) - 1);
    long long mask = -(long long)((d >> (__BITS_FOR_MINMAX(long long) - 1)) ^ ov);
    return (a & ~mask) | (b & mask);
}

/* ── unsigned helpers (flip-MSB trick) ──────────────────────────────────── */
/*
 * XOR with SIGN_BIT maps the unsigned total order to the signed total order:
 *   a <(u) b   iff   (a ^ SIGN_BIT) <(s) (b ^ SIGN_BIT)
 * After the flip, the signed formula above yields the correct mask.
 */
static inline unsigned int __min_uint(unsigned int a, unsigned int b) {
    unsigned int F  = 1u << (__BITS_FOR_MINMAX(unsigned int) - 1);
    unsigned int sa = a ^ F, sb = b ^ F, d = sa - sb;
    unsigned int ov = ((sa ^ sb) & (sa ^ d)) >> (__BITS_FOR_MINMAX(unsigned int) - 1);
    unsigned int mask = -((d >> (__BITS_FOR_MINMAX(unsigned int) - 1)) ^ ov);
    return (a & mask) | (b & ~mask);
}
static inline unsigned int __max_uint(unsigned int a, unsigned int b) {
    unsigned int F  = 1u << (__BITS_FOR_MINMAX(unsigned int) - 1);
    unsigned int sa = a ^ F, sb = b ^ F, d = sa - sb;
    unsigned int ov = ((sa ^ sb) & (sa ^ d)) >> (__BITS_FOR_MINMAX(unsigned int) - 1);
    unsigned int mask = -((d >> (__BITS_FOR_MINMAX(unsigned int) - 1)) ^ ov);
    return (a & ~mask) | (b & mask);
}

static inline unsigned long __min_ulong(unsigned long a, unsigned long b) {
    unsigned long F  = 1ul << (__BITS_FOR_MINMAX(unsigned long) - 1);
    unsigned long sa = a ^ F, sb = b ^ F, d = sa - sb;
    unsigned long ov = ((sa ^ sb) & (sa ^ d)) >> (__BITS_FOR_MINMAX(unsigned long) - 1);
    unsigned long mask = -((d >> (__BITS_FOR_MINMAX(unsigned long) - 1)) ^ ov);
    return (a & mask) | (b & ~mask);
}
static inline unsigned long __max_ulong(unsigned long a, unsigned long b) {
    unsigned long F  = 1ul << (__BITS_FOR_MINMAX(unsigned long) - 1);
    unsigned long sa = a ^ F, sb = b ^ F, d = sa - sb;
    unsigned long ov = ((sa ^ sb) & (sa ^ d)) >> (__BITS_FOR_MINMAX(unsigned long) - 1);
    unsigned long mask = -((d >> (__BITS_FOR_MINMAX(unsigned long) - 1)) ^ ov);
    return (a & ~mask) | (b & mask);
}

static inline unsigned long long __min_ullong(unsigned long long a, unsigned long long b) {
    unsigned long long F  = 1ull << (__BITS_FOR_MINMAX(unsigned long long) - 1);
    unsigned long long sa = a ^ F, sb = b ^ F, d = sa - sb;
    unsigned long long ov = ((sa ^ sb) & (sa ^ d)) >> (__BITS_FOR_MINMAX(unsigned long long) - 1);
    unsigned long long mask = -((d >> (__BITS_FOR_MINMAX(unsigned long long) - 1)) ^ ov);
    return (a & mask) | (b & ~mask);
}
static inline unsigned long long __max_ullong(unsigned long long a, unsigned long long b) {
    unsigned long long F  = 1ull << (__BITS_FOR_MINMAX(unsigned long long) - 1);
    unsigned long long sa = a ^ F, sb = b ^ F, d = sa - sb;
    unsigned long long ov = ((sa ^ sb) & (sa ^ d)) >> (__BITS_FOR_MINMAX(unsigned long long) - 1);
    unsigned long long mask = -((d >> (__BITS_FOR_MINMAX(unsigned long long) - 1)) ^ ov);
    return (a & ~mask) | (b & mask);
}

/* ── Master C11 macros ──────────────────────────────────────────────────── */
/*
 * Dispatch on typeof((a)+(b)) after usual arithmetic conversions, so
 * mixed-type pairs naturally route to the wider/promoted type:
 *   int + long long  → long long  → __min_llong   ✓
 *   int + float      → float      → fminf          ✓
 *   int + uint       → uint       → __min_uint     ✓
 *
 * The controlling expression (a)+(b) is NOT evaluated (C11 §6.5.1.1),
 * so arguments with side-effects are evaluated exactly once (in the
 * selected branch only).
 */
#define min(a, b) _Generic((a) + (b),                                               \
    long double:        fminl((a), (b)),                                             \
    double:             fmin((a), (b)),                                              \
    float:              fminf((a), (b)),                                             \
    unsigned long long: __min_ullong((unsigned long long)(a), (unsigned long long)(b)), \
    long long:          __min_llong((long long)(a), (long long)(b)),                 \
    unsigned long:      __min_ulong((unsigned long)(a), (unsigned long)(b)),         \
    long:               __min_long((long)(a), (long)(b)),                            \
    unsigned int:       __min_uint((unsigned int)(a), (unsigned int)(b)),            \
    default:            __min_int((int)(a), (int)(b))                                \
)

#define max(a, b) _Generic((a) + (b),                                               \
    long double:        fmaxl((a), (b)),                                             \
    double:             fmax((a), (b)),                                              \
    float:              fmaxf((a), (b)),                                             \
    unsigned long long: __max_ullong((unsigned long long)(a), (unsigned long long)(b)), \
    long long:          __max_llong((long long)(a), (long long)(b)),                 \
    unsigned long:      __max_ulong((unsigned long)(a), (unsigned long)(b)),         \
    long:               __max_long((long)(a), (long)(b)),                            \
    unsigned int:       __max_uint((unsigned int)(a), (unsigned int)(b)),            \
    default:            __max_int((int)(a), (int)(b))                                \
)

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* AC_F2C_INT_H */

