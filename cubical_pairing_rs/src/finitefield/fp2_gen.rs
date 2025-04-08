#![allow(non_snake_case)]

// NOTE:
// The majority of this code was written by Thomas Pornin, which is part of a
// larger work in progress implementing another project which is still ongoing.
// It is similar in form to the macros in the cryptographic library for rust:
// crrl. https://github.com/pornin/crrl
// This code was also used in the (2,2)-isogeny implementation
// https://github.com/ThetaIsogenies/two-isogenies

// A macro to define the finite field Fp^2 with modulus x^2 + 1 = 0 (uses that p
// = 3 mod 4) given a finite field of type Fp. All functions are designed to run
// in constant time. To construct the finite field, see the examples in
// fields.rs which include three distinct cases. To generate your own constants
// for the macro, the sage file `gen_fp.sage` will compute everything needed
// given a prime modulus.

// Macro expectations:
// A finite field Fp = GF(p) with p = 3 mod 4
// NQR_RE a Fp type such that (i + NQR_RE) is a NQR
macro_rules! define_fp2_core {
    () => {
        use core::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};
        use num_bigint::{BigInt, Sign};
        use rand_core::{CryptoRng, RngCore};
        use std::fmt;

        /// GF(p^2) implementation.
        #[derive(Clone, Copy, Debug)]
        pub struct Fp2 {
            x0: Fp,
            x1: Fp,
        }

        impl Fp2 {
            pub const ZERO: Self = Self {
                x0: Fp::ZERO,
                x1: Fp::ZERO,
            };
            pub const ONE: Self = Self {
                x0: Fp::ONE,
                x1: Fp::ZERO,
            };
            pub const TWO: Self = Self {
                x0: Fp::TWO,
                x1: Fp::ZERO,
            };
            pub const THREE: Self = Self {
                x0: Fp::THREE,
                x1: Fp::ZERO,
            };
            pub const FOUR: Self = Self {
                x0: Fp::FOUR,
                x1: Fp::ZERO,
            };
            pub const MINUS_ONE: Self = Self {
                x0: Fp::MINUS_ONE,
                x1: Fp::ZERO,
            };
            pub const ZETA: Self = Self {
                x0: Fp::ZERO,
                x1: Fp::ONE,
            };
            pub const MINUS_ZETA: Self = Self {
                x0: Fp::ZERO,
                x1: Fp::MINUS_ONE,
            };

            pub const ENCODED_LENGTH: usize = 2 * Fp::ENCODED_LENGTH;
            pub const BIT_LENGTH: usize = Fp::BIT_LENGTH;

            /// Non-quadratic residue.
            pub const NQR: Self = Self {
                x0: NQR_RE,
                x1: Fp::ONE,
            };

            pub const fn new(re: &Fp, im: &Fp) -> Self {
                Self { x0: *re, x1: *im }
            }

            #[inline]
            pub fn iszero(self) -> u32 {
                self.x0.iszero() & self.x1.iszero()
            }

            #[inline]
            pub fn equals(self, rhs: &Self) -> u32 {
                self.x0.equals(&rhs.x0) & self.x1.equals(&rhs.x1)
            }

            #[inline]
            fn set_add(&mut self, rhs: &Self) {
                self.x0 += &rhs.x0;
                self.x1 += &rhs.x1;
            }

            #[inline]
            fn set_sub(&mut self, rhs: &Self) {
                self.x0 -= &rhs.x0;
                self.x1 -= &rhs.x1;
            }

            #[inline]
            pub fn set_neg(&mut self) {
                self.x0.set_neg();
                self.x1.set_neg();
            }

            #[inline]
            pub fn set_conj(&mut self) {
                self.x1.set_neg();
            }

            #[inline]
            pub fn conj(self) -> Self {
                Self {
                    x0: self.x0,
                    x1: -&self.x1,
                }
            }

            #[inline]
            // NOTE: old multiplication has been replaced in favour of Longa's
            // algorithm which efficiently computes sums and differences of
            // products. For more information, see the `sum_of_products()`
            // function in `fp_gen.rs`.
            fn set_mul_old(&mut self, rhs: &Self) {
                // a <- x0*y0
                // b <- x1*y1
                // c <- (x0 + x1)*(y0 + y1)
                // (x0 + i*x1)*(y0 + i*y1) = (x0*y0 - x1*y1) + i*(x0*y1 + y0*x1)
                //                         = (a - b) + i*(c - a - b)
                let a = &self.x0 * &rhs.x0;
                let b = &self.x1 * &rhs.x1;
                let c = &(&self.x0 + &self.x1) * &(&rhs.x0 + &rhs.x1);
                self.x0 = a;
                self.x0 -= &b;
                self.x1 = c;
                self.x1 -= &a;
                self.x1 -= &b;
            }

            #[inline]
            pub fn mul_old(self, rhs: Self) -> Self {
                let mut r = self;
                r.set_mul_old(&rhs);
                r
            }

            #[inline(always)]
            pub fn set_mul(&mut self, rhs: &Self) {
                // Computes x*y from:
                // x = (x0 + i*x1)
                // y = (y0 + i*y1)
                // x*y = (x0 + i*x1)*(y0 + i*y1)
                //     = (x0*y0 - x1*y1) + i*(x0*y1 + y0*x1)
                // Computes (x0*y0 - x1*y1)
                let x0 = Fp::difference_of_products(&self.x0, &rhs.x0, &self.x1, &rhs.x1);
                // Computes (x0*y1 + y0*x1)
                let x1 = Fp::sum_of_products(&self.x0, &rhs.x1, &self.x1, &rhs.x0);

                self.x0 = x0;
                self.x1 = x1;
            }

            #[inline]
            pub fn mul_new(self, rhs: Self) -> Self {
                let mut r = self;
                r.set_mul(&rhs);
                r
            }

            #[inline]
            pub fn set_square(&mut self) {
                // (x0 + i*x1)^2 = (x0^2 - x1^2) + 2*i*(x0*x1)
                //               = (x0 + x1)*(x0 - x1) + i*(2*x0*x1)
                let a = &self.x0 + &self.x1;
                let b = &self.x0 - &self.x1;
                self.x1 *= &self.x0;
                self.x1.set_mul2();
                self.x0 = a;
                self.x0 *= &b;
            }

            #[inline]
            pub fn square(self) -> Self {
                let mut r = self;
                r.set_square();
                r
            }

            #[inline]
            pub fn set_half(&mut self) {
                self.x0.set_half();
                self.x1.set_half();
            }

            #[inline]
            pub fn half(self) -> Self {
                let mut r = self;
                r.set_half();
                r
            }

            #[inline]
            pub fn set_mul2(&mut self) {
                self.x0.set_mul2();
                self.x1.set_mul2();
            }

            #[inline]
            pub fn mul2(self) -> Self {
                let mut r = self;
                r.set_mul2();
                r
            }

            #[inline]
            pub fn set_mul3(&mut self) {
                self.x0.set_mul3();
                self.x1.set_mul3();
            }

            #[inline]
            pub fn mul3(self) -> Self {
                let mut r = self;
                r.set_mul3();
                r
            }

            #[inline]
            pub fn set_mul4(&mut self) {
                self.x0.set_mul4();
                self.x1.set_mul4();
            }

            #[inline]
            pub fn mul4(self) -> Self {
                let mut r = self;
                r.set_mul4();
                r
            }

            #[inline]
            pub fn set_mul8(&mut self) {
                self.x0.set_mul8();
                self.x1.set_mul8();
            }

            #[inline]
            pub fn mul8(self) -> Self {
                let mut r = self;
                r.set_mul8();
                r
            }

            #[inline]
            pub fn conjugate(self) -> Self {
                let mut r = self;
                r.x1.set_neg();
                r
            }

            #[inline]
            pub fn set_mul_small(&mut self, k: i32) {
                self.x0.set_mul_small(k);
                self.x1.set_mul_small(k);
            }

            #[inline]
            pub fn mul_small(self, k: i32) -> Self {
                let mut r = self;
                r.set_mul_small(k);
                r
            }

            #[inline]
            pub fn set_select(&mut self, a: &Self, b: &Self, ctl: u32) {
                self.x0.set_select(&a.x0, &b.x0, ctl);
                self.x1.set_select(&a.x1, &b.x1, ctl);
            }

            #[inline]
            pub fn select(a: &Self, b: &Self, ctl: u32) -> Self {
                Self {
                    x0: Fp::select(&a.x0, &b.x0, ctl),
                    x1: Fp::select(&a.x1, &b.x1, ctl),
                }
            }

            #[inline]
            pub fn set_cond(&mut self, rhs: &Self, ctl: u32) {
                self.x0.set_cond(&rhs.x0, ctl);
                self.x1.set_cond(&rhs.x1, ctl);
            }

            #[inline]
            pub fn set_condneg(&mut self, ctl: u32) {
                let y0 = -(&self.x0);
                let y1 = -(&self.x1);
                self.x0.set_cond(&y0, ctl);
                self.x1.set_cond(&y1, ctl);
            }

            #[inline]
            pub fn condswap(a: &mut Self, b: &mut Self, ctl: u32) {
                Fp::condswap(&mut a.x0, &mut b.x0, ctl);
                Fp::condswap(&mut a.x1, &mut b.x1, ctl);
            }

            #[inline]
            fn set_div(&mut self, rhs: &Self) {
                // 1/(x0 + i*x1) = (x0 - i*x1)/(x0^2 + x1^2)
                let mut z = rhs.x0.square();
                z += &rhs.x1.square();
                z.set_invert();
                let mut r = *rhs;
                r.x1.set_neg();
                r.x0 *= &z;
                r.x1 *= &z;
                self.set_mul(&r);
            }

            #[inline]
            pub fn set_invert(&mut self) {
                // 1/(x0 + i*x1) = (x0 - i*x1)/(x0^2 + x1^2)
                let mut z = self.x0.square();
                z += &self.x1.square();
                z.set_invert();
                self.x0 *= &z;
                self.x1 *= &z;
                self.x1.set_neg();
            }

            #[inline]
            pub fn invert(self) -> Self {
                let mut r = self;
                r.set_invert();
                r
            }

            /// Legendre symbol on this value. Return value is:
            ///   0   if this value is zero
            ///  +1   if this value is a non-zero quadratic residue
            ///  -1   if this value is not a quadratic residue
            #[inline]
            pub fn legendre(self) -> i32 {
                // x = x0 + i*x1 is a square in GF(p^2) if and only if
                // x0^2 + x1^2 is a square in GF(p). Moreover, x0^2 + x1^2 is
                // zero if and only if x is zero.
                (self.x0.square() + self.x1.square()).legendre()
            }

            /// Set this value to its square root. Returned value is 0xFFFFFFFF if
            /// the operation succeeded (value was indeed a quadratic residue), or
            /// 0x00000000 otherwise. On success, the chosen root is the one whose
            /// sign is 0 (i.e. if the "real part" is non-zero, then it is an even
            /// integer; if the "real part" is zero, then the "imaginary part" is
            /// an even integer). On failure, this value is set to 0.
            pub fn set_sqrt(&mut self) -> u32 {
                // x^p = (x0 + i*x1)^p = x0 - i*x1  (Frobenius automorphism)
                // Thus: x^(p+1) = (x0 + i*x1)*(x0 - i*x1) = x0^2 + x1^2, which
                // is an element of GF(p). All elements of GF(p) are squares in
                // GF(p^2), but x0^2 + x1^2 is not necessarily a square in GF(p).
                //
                // Let conj(p) = x^p = x0 - i*x1. Note that conj() is analogous to
                // the conjugate in complex numbers. In particular:
                //    conj(a + b) = conj(a) + conj(b)
                //    conj(a * b) = conj(a) * conj(b)
                // This implies that conj(x) is a square if and only if x is a
                // square, and conj(sqrt(x)) = sqrt(conj(x)). Thus, if x is a
                // square, then:
                //    (sqrt(x)*conj(sqrt(x)))^2 = x*conj(x) = x0^2 + x1^2
                // But sqrt(x)*conj(sqrt(x)) is in GF(p); therefore, if x is a
                // square, then x0^2 + x1^2 must be a square in GF(p).
                //
                // Suppose that y = y0 + i*y1 such that y^2 = x. Then:
                //   y0^2 - y1^2 = x0
                //   2*y0*y1 = x1
                // If x1 = 0 then:
                //    if x0.legendre() >= 0 then y = sqrt(x0)
                //                          else y = i*sqrt(-x0)
                // else:
                //    y0 != 0 (necessarily) and y1 = x1 / (2*y0)
                //    Thus:
                //       y0^4 - x0*y0^2 - (x1^2)/4 = 0
                //    Discriminant is delta = x0^2 + x1^2, which is always a square
                //    (see above). Therefore:
                //       y0^2 = (x0 +/- sqrt(delta))/2
                //    We can thus compute (x0 + sqrt(delta))/2 and check its
                //    Legendre symbol; we subtract sqrt(delta) from it if it is
                //    not a square. We then extract y0 as a square root of the
                //    result, and compute y1 from it.
                //
                // Main cost is the two square roots in GF(p) (for delta and
                // for y0); Legendre symbols and inversions are vastly faster.

                // sqrt_delta <- sqrt(x0^2 + x1^2)
                let (sqrt_delta, r1) = (self.x0.square() + self.x1.square()).sqrt();
                // y0sq <- (x0 + sqrt(delta)) / 2
                let mut y0sq = (self.x0 + sqrt_delta).half();
                // If x1 = 0, then replace y0sq with x0
                let x1z = self.x1.iszero();
                y0sq.set_cond(&self.x0, x1z);
                // Get the Legendre symbol and set nqr to 0xFFFFFFFF when y0sq
                // is not a square
                let ls = y0sq.legendre();
                let nqr = (ls >> 1) as u32;
                // If not a square:
                //    if x1 = 0, then y0sq contains x0 and we want -x0
                //    if x1 != 0, then y0sq <- y0sq - sqrt(delta)
                y0sq.set_condneg(nqr & x1z);
                y0sq.set_cond(&(y0sq - sqrt_delta), nqr & !x1z);
                // Get the square root.
                let (mut y0, r2) = y0sq.sqrt();
                let r = r1 & r2;
                // Compute y1 = x1 / (2*y0).
                let mut y1 = self.x1 / y0.mul2();
                // If x1 = 0, then the square root worked, and y1 = 0 at this point;
                // we must still exchange y0 and y1 if x0 was not a square.
                Fp::condswap(&mut y0, &mut y1, nqr & x1z);
                // Result goes into this object. If there was a failure (r == 0),
                // then we must clear both x0 and x1.
                self.x0.set_select(&Fp::ZERO, &y0, r);
                self.x1.set_select(&Fp::ZERO, &y1, r);
                // Sign mangement: negate the result if needed.
                let x0odd = ((self.x0.encode()[0] as u32) & 1).wrapping_neg();
                let x1odd = ((self.x1.encode()[0] as u32) & 1).wrapping_neg();
                let x0z = self.x0.iszero();
                self.set_condneg(x0odd | (x0z & x1odd));
                r
            }

            pub fn sqrt(self) -> (Self, u32) {
                let mut y = self;
                let r = y.set_sqrt();
                (y, r)
            }

            /// Raise this value to the power e. Exponent e is encoded in
            /// unsigned little-endian convention over exactly ebitlen bits.
            pub fn set_pow(&mut self, e: &[u8], ebitlen: usize) {
                self.set_pow_ext(e, 0, ebitlen);
            }

            /// Raise this value to the power e. Exponent e is encoded in
            /// unsigned little-endian convention, over exactly ebitlen bits,
            /// and starting at the bit offset eoff.
            pub fn set_pow_ext(&mut self, e: &[u8], eoff: usize, ebitlen: usize) {
                // TODO: implement a window optimization to make fewer
                // multiplications.
                let x = *self;
                *self = Self::ONE;
                for i in (eoff..(eoff + ebitlen)).rev() {
                    let y = &*self * &x;
                    let ctl = (((e[i >> 3] >> (i & 7)) as u32) & 1).wrapping_neg();
                    self.set_cond(&y, ctl);
                    if i == eoff {
                        break;
                    }
                    self.set_square();
                }
            }

            /// Return this value to the power e (as a new element). Exponent e
            /// is encoded in unsigned little-endian convention over exactly
            /// ebitlen bits.
            pub fn pow(self, e: &[u8], ebitlen: usize) -> Self {
                let mut x = self;
                x.set_pow(e, ebitlen);
                x
            }

            /// Return this value to the power e (as a new element). Exponent e
            /// is encoded in unsigned little-endian convention over exactly
            /// ebitlen bits, and starting at the bit offset eoff.
            pub fn pow_ext(self, e: &[u8], eoff: usize, ebitlen: usize) -> Self {
                let mut x = self;
                x.set_pow_ext(e, eoff, ebitlen);
                x
            }

            pub fn encode(self) -> [u8; Self::ENCODED_LENGTH] {
                let mut r = [0u8; Self::ENCODED_LENGTH];
                r[..Fp::ENCODED_LENGTH].copy_from_slice(&self.x0.encode());
                r[Fp::ENCODED_LENGTH..].copy_from_slice(&self.x1.encode());
                r
            }

            pub fn decode(buf: &[u8]) -> (Self, u32) {
                if buf.len() != Self::ENCODED_LENGTH {
                    return (Self::ZERO, 0);
                }
                let (mut x0, c0) = Fp::decode(&buf[..Fp::ENCODED_LENGTH]);
                let (mut x1, c1) = Fp::decode(&buf[Fp::ENCODED_LENGTH..]);
                let cx = c0 & c1;
                x0.set_cond(&Fp::ZERO, !cx);
                x1.set_cond(&Fp::ZERO, !cx);
                (Self { x0, x1 }, cx)
            }

            /// Set this structure to a random field element (indistinguishable
            /// from uniform generation).
            pub fn set_rand<T: CryptoRng + RngCore>(&mut self, rng: &mut T) {
                self.x0.set_rand(rng);
                self.x1.set_rand(rng);
            }

            /// Return a new random field element (indistinguishable from
            /// uniform generation).
            pub fn rand<T: CryptoRng + RngCore>(rng: &mut T) -> Self {
                let mut x = Self::ZERO;
                x.set_rand(rng);
                x
            }

            /// Set this structure to a random non-square field element
            /// (indistinguishable from uniform generation).
            pub fn set_rand_nonsquare<T: CryptoRng + RngCore>(&mut self, rng: &mut T) {
                // We get a random non-square by getting a random non-zero
                // value, squaring it, and multiplying by a known non-square.
                // The loop handles the case of getting zero randomly, which
                // is very unlikely in practice.
                loop {
                    self.set_rand(rng);
                    if self.iszero() != 0 {
                        continue;
                    }
                    break;
                }
                self.set_square();
                self.set_mul(&Self::NQR);
            }

            /// Return a new random non-square field element
            /// (indistinguishable from uniform generation).
            pub fn rand_nonsquare<T: CryptoRng + RngCore>(self, rng: &mut T) -> Self {
                let mut x = Self::ZERO;
                x.set_rand_nonsquare(rng);
                x
            }

            /// Raise this value to the power e. The exponent length (in bits)
            /// MUST be at most ebitlen. This is constant-time for both the
            /// base value (self) and the exponent (e); the exponent maximum
            /// size (ebitlen) is considered non-secret.
            pub fn set_pow_u64(&mut self, e: u64, ebitlen: usize) {
                match ebitlen {
                    0 => {
                        *self = Self::ONE;
                    }
                    1 => {
                        self.set_cond(&Self::ONE, ((e as u32) & 1).wrapping_sub(1));
                    }
                    _ => {
                        let x = *self;
                        self.set_cond(
                            &Self::ONE,
                            (((e >> (ebitlen - 1)) as u32) & 1).wrapping_sub(1),
                        );
                        for i in (0..(ebitlen - 1)).rev() {
                            self.set_square();
                            let y = &*self * &x;
                            self.set_cond(&y, (((e >> i) as u32) & 1).wrapping_neg());
                        }
                    }
                }
            }

            /// Return this value to the power e. The exponent length (in bits)
            /// MUST be at most ebitlen. This is constant-time for both the
            /// base value (self) and the exponent (e); the exponent maximum
            /// size (ebitlen) is considered non-secret.
            pub fn pow_u64(self, e: u64, ebitlen: usize) -> Self {
                let mut x = self;
                x.set_pow_u64(e, ebitlen);
                x
            }

            /// Raise this value to the power e. The exponent is considered
            /// non-secret.
            pub fn set_pow_u64_vartime(&mut self, e: u64) {
                match e {
                    0 => {
                        *self = Self::ONE;
                    }
                    1 => {
                        return;
                    }
                    2 => {
                        self.set_square();
                    }
                    3 => {
                        *self *= self.square();
                    }
                    4 => {
                        self.set_square();
                        self.set_square();
                    }
                    _ => {
                        let xx = self.square();
                        let xw = [*self, xx, xx * &*self];
                        let mut j = 63 - e.leading_zeros();
                        j &= !1u32;
                        *self = xw[((e >> j) as usize) - 1];
                        while j > 0 {
                            j -= 2;
                            self.set_square();
                            self.set_square();
                            let k = ((e >> j) as usize) & 3;
                            if k > 0 {
                                self.set_mul(&xw[k - 1]);
                            }
                        }
                    }
                }
            }

            /// Return this value to the power e. The exponent is considered
            /// non-secret.
            pub fn pow_u64_vartime(self, e: u64) -> Self {
                let mut x = self;
                x.set_pow_u64_vartime(e);
                x
            }

            /// Get the "hash" of the value. For x = x0 + i*x1, this is:
            ///    (hashcode(x0) << 1) | (hashcode(x1) & 1)
            /// i.e. bit 0 is bit 0 of x1, and bits 1..63 are bits 0..62 of x0
            /// (both in Montgomery representation).
            pub fn hashcode(self) -> u64 {
                (self.x0.hashcode() << 1) | (self.x1.hashcode() & 1)
            }

            pub fn batch_invert(xx: &mut [Self]) {
                // We use Montgomery's trick:
                //   1/u = v*(1/(u*v))
                //   1/v = u*(1/(u*v))
                // Applied recursively on n elements, this computes an inversion
                // with a single inversion in the field, and 3*(n-1) multiplications.
                // We use batches of 200 elements; larger batches only yield
                // moderate improvements, while sticking to a fixed moderate batch
                // size allows stack-based allocation.
                let n = xx.len();
                let mut i = 0;
                while i < n {
                    let blen = if (n - i) > 200 { 200 } else { n - i };
                    let mut tt = [Self::ZERO; 200];
                    tt[0] = xx[i];
                    let zz0 = tt[0].iszero();
                    tt[0].set_cond(&Self::ONE, zz0);
                    for j in 1..blen {
                        tt[j] = xx[i + j];
                        tt[j].set_cond(&Self::ONE, tt[j].iszero());
                        tt[j] *= tt[j - 1];
                    }
                    let mut k = Self::ONE / tt[blen - 1];
                    for j in (1..blen).rev() {
                        let mut x = xx[i + j];
                        let zz = x.iszero();
                        x.set_cond(&Self::ONE, zz);
                        xx[i + j].set_cond(&(k * tt[j - 1]), !zz);
                        k *= x;
                    }
                    xx[i].set_cond(&k, !zz0);
                    i += blen;
                }
            }

            /// Inner function for solving DLP with order 2^e.
            /// If gk = -1, then base is self; otherwise, it is gpp[gk]
            /// (equal to self^(2^dlog_table[gk])). Order of the base is 2^lg.
            /// Output (of size lg bits) is written in v[], starting at offset
            /// voff (counted in bits). The target bit values MUST be all zero
            /// initially in v[] (non-target bits are not modified).
            /// Returned value is 0xFFFFFFFF on success, 0x00000000 on error.
            fn solve_dlp_n_inner(
                self,
                gpp: &[Self],
                gk: usize,
                x: &Self,
                v: &mut [u8],
                voff: usize,
                e: usize,
                dlog_table: &[usize],
            ) -> u32 {
                let lg = e - dlog_table[gk];

                // At the deepest recursion level, lg = 1, g = -1,
                // and x = 1 or -1.
                if lg == 1 {
                    let hz = x.x1.iszero();
                    let lp = x.x0.equals(&Fp::ONE);
                    let ln = x.x0.equals(&Fp::MINUS_ONE);
                    v[voff >> 3] |= ((ln & 1) << (voff & 7)) as u8;
                    return hz & (lp | ln);
                }

                // Split lg = lg0 + lg1.
                // Precomputed indices (in dlog_table) assume that the split is
                // done such that lg0 = floor(lg/2).
                let lg0 = lg >> 1;
                let lg1 = lg - lg0;

                // Solve for v0.
                //   g' = g^(2^lg1)
                //   x' = x^(2^lg1)
                let mut gk0 = gk + 1;
                while dlog_table[gk0] != e - lg0 {
                    gk0 += 1;
                }
                let mut x0 = *x;
                for _ in 0..lg1 {
                    x0.set_square();
                }
                let ok0 = self.solve_dlp_n_inner(gpp, gk0, &x0, v, voff, e, dlog_table);

                // Solve for v1.
                //   g' = g^(2^lg0)
                //   x' = x/g^v0
                let mut gk1 = gk + 1;
                while dlog_table[gk1] != e - lg1 {
                    gk1 += 1;
                }
                let mut x1 = gpp[gk].conj();
                x1.set_pow_ext(v, voff, lg0);
                x1 *= x;
                let ok1 = self.solve_dlp_n_inner(gpp, gk1, &x1, v, voff + lg0, e, dlog_table);

                ok0 & ok1
            }

            /// Precompute a look-up table for computing discrete logs
            /// Use the function precompute_dlp_table to precompute the
            /// values g^(2^lg), keeping the relevant values
            /// in the gpp[] array. We have:
            ///    gpp[j] = g^(2^dlog_table[j])
            /// Note that the first value (gpp[0]) is g itself, and the
            /// last one must be -1 (otherwise, g does not have order
            /// exactly 2^e).
            pub fn precompute_dlp_table(self, dlog_table: &[usize]) -> (Vec<Self>, u32) {
                let mut gpp = vec![Self::ZERO; dlog_table.len()];
                gpp[0] = self;
                let mut j = 1;
                let mut g = self;
                let mut lg = 0;
                while j < gpp.len() {
                    g.set_square();
                    lg += 1;
                    if lg == dlog_table[j] {
                        gpp[j] = g;
                        j += 1;
                    }
                }

                // Also check that g has indeed order exactly n.
                let ok1 = g.equals(&Self::MINUS_ONE);

                (gpp, ok1)
            }

            /// Find integer v (modulo 2^e) such that x = self^v. If self
            /// has order exactly 2^e, and there is a solution v, then this
            /// function returns (v, 0xFFFFFFFF). If self does not have order
            /// exactly 2^e (including if self^(2^(e-1)) = 1, i.e. the order of
            /// self is a strict divisor or 2^e), or if there is no solution,
            /// then this function returns (0, 0).
            pub fn solve_dlp_2e(
                self,
                x: &Self,
                e: usize,
                gpp: &[Self],
                dlog_table: &[usize],
            ) -> (Vec<u8>, u32) {
                // Method: consider g, x and lg such that:
                //   g has multiplicative order 2^lg
                //   x = g^v for some v (in the 0 to 2^lg-1 range)
                // If lg = 1 then g = -1, and x = 0 or -1.
                //   -> if g != -1, or x is not 0 or -1, then the input is
                //      erroneous and we can report it
                // If lg > 1:
                //   Let lg0 = floor(lg / 2) and lg1 = lg - lg0.
                //   Let v = v0 + (2^lg0)*v1
                //   Then:
                //      x^(2^lg1) = (g^(2^lg1))^v0
                //   We get v0 with a recursive call on base g^(2^lg1) and
                //   value x^(2^lg1). Once we have v0:
                //      x/g^v0 = (g^(2^lg0))^v1
                //   Another recursive call on base g^(2^lg0) and value
                //   x/g^v0 yields v1, from which we easily obtain v.
                //   Note that 1/g = conj(g), since g is a 2n-th root of 1.
                //
                // We avoid recomputing the same g^(2^lg) values by keeping
                // the relevant values in a local array; the important indices
                // are the ones specified in the dlog_table array.

                // Use the function precompute_dlp_table to precompute the
                // values g^(2^lg), keeping the relevant values
                // in the gpp[] array. We have:
                //    gpp[j] = g^(2^dlog_table[j])
                // Note that the first value (gpp[0]) is g itself, and the
                // last one must be -1 (otherwise, g does not have order
                // exactly n).

                // Apply the recursion.
                let mut v = vec![0u8; (e + 7) >> 3];
                let ok = self.solve_dlp_n_inner(gpp, 0, x, &mut v, 0, e, dlog_table);
                (v, ok)
            }
        }

        // ========================================================================

        impl fmt::Display for Fp2 {
            fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
                let r = self.encode();

                let x0_bytes = &r[..Fp::ENCODED_LENGTH];
                let x1_bytes = &r[Fp::ENCODED_LENGTH..];

                let x0_big = BigInt::from_bytes_le(Sign::Plus, x0_bytes);
                let x1_big = BigInt::from_bytes_le(Sign::Plus, x1_bytes);

                write!(f, "i*{} + {}", x1_big, x0_big)
            }
        }

        impl Add<Fp2> for Fp2 {
            type Output = Fp2;

            #[inline(always)]
            fn add(self, other: Fp2) -> Fp2 {
                let mut r = self;
                r.set_add(&other);
                r
            }
        }

        impl Add<&Fp2> for Fp2 {
            type Output = Fp2;

            #[inline(always)]
            fn add(self, other: &Fp2) -> Fp2 {
                let mut r = self;
                r.set_add(other);
                r
            }
        }

        impl Add<Fp2> for &Fp2 {
            type Output = Fp2;

            #[inline(always)]
            fn add(self, other: Fp2) -> Fp2 {
                let mut r = *self;
                r.set_add(&other);
                r
            }
        }

        impl Add<&Fp2> for &Fp2 {
            type Output = Fp2;

            #[inline(always)]
            fn add(self, other: &Fp2) -> Fp2 {
                let mut r = *self;
                r.set_add(other);
                r
            }
        }

        impl AddAssign<Fp2> for Fp2 {
            #[inline(always)]
            fn add_assign(&mut self, other: Fp2) {
                self.set_add(&other);
            }
        }

        impl AddAssign<&Fp2> for Fp2 {
            #[inline(always)]
            fn add_assign(&mut self, other: &Fp2) {
                self.set_add(other);
            }
        }

        impl Div<Fp2> for Fp2 {
            type Output = Fp2;

            #[inline(always)]
            fn div(self, other: Fp2) -> Fp2 {
                let mut r = self;
                r.set_div(&other);
                r
            }
        }

        impl Div<&Fp2> for Fp2 {
            type Output = Fp2;

            #[inline(always)]
            fn div(self, other: &Fp2) -> Fp2 {
                let mut r = self;
                r.set_div(other);
                r
            }
        }

        impl Div<Fp2> for &Fp2 {
            type Output = Fp2;

            #[inline(always)]
            fn div(self, other: Fp2) -> Fp2 {
                let mut r = *self;
                r.set_div(&other);
                r
            }
        }

        impl Div<&Fp2> for &Fp2 {
            type Output = Fp2;

            #[inline(always)]
            fn div(self, other: &Fp2) -> Fp2 {
                let mut r = *self;
                r.set_div(other);
                r
            }
        }

        impl DivAssign<Fp2> for Fp2 {
            #[inline(always)]
            fn div_assign(&mut self, other: Fp2) {
                self.set_div(&other);
            }
        }

        impl DivAssign<&Fp2> for Fp2 {
            #[inline(always)]
            fn div_assign(&mut self, other: &Fp2) {
                self.set_div(other);
            }
        }

        impl Mul<Fp2> for Fp2 {
            type Output = Fp2;

            #[inline(always)]
            fn mul(self, other: Fp2) -> Fp2 {
                let mut r = self;
                r.set_mul(&other);
                r
            }
        }

        impl Mul<&Fp2> for Fp2 {
            type Output = Fp2;

            #[inline(always)]
            fn mul(self, other: &Fp2) -> Fp2 {
                let mut r = self;
                r.set_mul(other);
                r
            }
        }

        impl Mul<Fp2> for &Fp2 {
            type Output = Fp2;

            #[inline(always)]
            fn mul(self, other: Fp2) -> Fp2 {
                let mut r = *self;
                r.set_mul(&other);
                r
            }
        }

        impl Mul<&Fp2> for &Fp2 {
            type Output = Fp2;

            #[inline(always)]
            fn mul(self, other: &Fp2) -> Fp2 {
                let mut r = *self;
                r.set_mul(other);
                r
            }
        }

        impl MulAssign<Fp2> for Fp2 {
            #[inline(always)]
            fn mul_assign(&mut self, other: Fp2) {
                self.set_mul(&other);
            }
        }

        impl MulAssign<&Fp2> for Fp2 {
            #[inline(always)]
            fn mul_assign(&mut self, other: &Fp2) {
                self.set_mul(other);
            }
        }

        impl Neg for Fp2 {
            type Output = Fp2;

            #[inline(always)]
            fn neg(self) -> Fp2 {
                let mut r = self;
                r.set_neg();
                r
            }
        }

        impl Neg for &Fp2 {
            type Output = Fp2;

            #[inline(always)]
            fn neg(self) -> Fp2 {
                let mut r = *self;
                r.set_neg();
                r
            }
        }

        impl Sub<Fp2> for Fp2 {
            type Output = Fp2;

            #[inline(always)]
            fn sub(self, other: Fp2) -> Fp2 {
                let mut r = self;
                r.set_sub(&other);
                r
            }
        }

        impl Sub<&Fp2> for Fp2 {
            type Output = Fp2;

            #[inline(always)]
            fn sub(self, other: &Fp2) -> Fp2 {
                let mut r = self;
                r.set_sub(other);
                r
            }
        }

        impl Sub<Fp2> for &Fp2 {
            type Output = Fp2;

            #[inline(always)]
            fn sub(self, other: Fp2) -> Fp2 {
                let mut r = *self;
                r.set_sub(&other);
                r
            }
        }

        impl Sub<&Fp2> for &Fp2 {
            type Output = Fp2;

            #[inline(always)]
            fn sub(self, other: &Fp2) -> Fp2 {
                let mut r = *self;
                r.set_sub(other);
                r
            }
        }

        impl SubAssign<Fp2> for Fp2 {
            #[inline(always)]
            fn sub_assign(&mut self, other: Fp2) {
                self.set_sub(&other);
            }
        }

        impl SubAssign<&Fp2> for Fp2 {
            #[inline(always)]
            fn sub_assign(&mut self, other: &Fp2) {
                self.set_sub(other);
            }
        }
    };
} // End of macro: define_fp2_core

pub(crate) use define_fp2_core;

// ========================================================================

// Macro expectations:
#[cfg(test)]
macro_rules! define_fp2_tests {
    () => {
        use super::{Fp, Fp2};
        use num_bigint::{BigInt, Sign};
        use sha2::{Digest, Sha256};

        fn check_fp2_ops(va: &[u8], vb: &[u8], with_sqrt: bool) {
            let mut zpww = [0u32; Fp::N * 2];
            for i in 0..Fp::N {
                zpww[2 * i] = Fp::MODULUS[i] as u32;
                zpww[2 * i + 1] = (Fp::MODULUS[i] >> 32) as u32;
            }
            let zp = BigInt::from_slice(Sign::Plus, &zpww);

            let alen = va.len() >> 1;
            let blen = vb.len() >> 1;

            let a0 = Fp::decode_reduce(&va[..alen]);
            let a1 = Fp::decode_reduce(&va[alen..]);
            let b0 = Fp::decode_reduce(&vb[..blen]);
            let b1 = Fp::decode_reduce(&vb[blen..]);
            let za0 = BigInt::from_bytes_le(Sign::Plus, &a0.encode());
            let za1 = BigInt::from_bytes_le(Sign::Plus, &a1.encode());
            let zb0 = BigInt::from_bytes_le(Sign::Plus, &b0.encode());
            let zb1 = BigInt::from_bytes_le(Sign::Plus, &b1.encode());
            let a = Fp2::new(&a0, &a1);
            let b = Fp2::new(&b0, &b1);

            let c = a + b;
            let vc = c.encode();
            let zc0 = BigInt::from_bytes_le(Sign::Plus, &vc[..Fp::ENCODED_LENGTH]);
            let zc1 = BigInt::from_bytes_le(Sign::Plus, &vc[Fp::ENCODED_LENGTH..]);
            let zd0 = (&za0 + &zb0) % &zp;
            let zd1 = (&za1 + &zb1) % &zp;
            assert!(zc0 == zd0 && zc1 == zd1);

            let c = a - b;
            let vc = c.encode();
            let zc0 = BigInt::from_bytes_le(Sign::Plus, &vc[..Fp::ENCODED_LENGTH]);
            let zc1 = BigInt::from_bytes_le(Sign::Plus, &vc[Fp::ENCODED_LENGTH..]);
            let zd0 = (&zp + &za0 - &zb0) % &zp;
            let zd1 = (&zp + &za1 - &zb1) % &zp;
            assert!(zc0 == zd0 && zc1 == zd1);

            let c = a * b;
            let vc = c.encode();
            let zc0 = BigInt::from_bytes_le(Sign::Plus, &vc[..Fp::ENCODED_LENGTH]);
            let zc1 = BigInt::from_bytes_le(Sign::Plus, &vc[Fp::ENCODED_LENGTH..]);
            let zd0 = (&zp + ((&za0 * &zb0) % &zp) - ((&za1 * &zb1) % &zp)) % &zp;
            let zd1 = ((&za0 * &zb1) + (&za1 * &zb0)) % &zp;
            assert!(zc0 == zd0 && zc1 == zd1);

            // let c = a.mul_new(b);
            // let vc = c.encode();
            // let zc0 = BigInt::from_bytes_le(Sign::Plus, &vc[..Fp::ENCODED_LENGTH]);
            // let zc1 = BigInt::from_bytes_le(Sign::Plus, &vc[Fp::ENCODED_LENGTH..]);
            // let zd0 = (&zp + ((&za0 * &zb0) % &zp) - ((&za1 * &zb1) % &zp)) % &zp;
            // let zd1 = ((&za0 * &zb1) + (&za1 * &zb0)) % &zp;
            // println!("4");
            // assert!(zc0 == zd0 && zc1 == zd1);

            let c = a.square();
            let vc = c.encode();
            let zc0 = BigInt::from_bytes_le(Sign::Plus, &vc[..Fp::ENCODED_LENGTH]);
            let zc1 = BigInt::from_bytes_le(Sign::Plus, &vc[Fp::ENCODED_LENGTH..]);
            let zd0 = (&zp + ((&za0 * &za0) % &zp) - ((&za1 * &za1) % &zp)) % &zp;
            let zd1 = ((&za0 * &za1) + (&za1 * &za0)) % &zp;
            assert!(zc0 == zd0 && zc1 == zd1);

            let c = a / b;
            if b.iszero() != 0 {
                assert!(c.iszero() == 0xFFFFFFFF);
            } else {
                let c = c * b;
                let vc = c.encode();
                let zc0 = BigInt::from_bytes_le(Sign::Plus, &vc[..Fp::ENCODED_LENGTH]);
                let zc1 = BigInt::from_bytes_le(Sign::Plus, &vc[Fp::ENCODED_LENGTH..]);
                assert!(zc0 == za0 && zc1 == za1);
            }

            let c = b.invert();
            if b.iszero() != 0 {
                assert!(c.iszero() == 0xFFFFFFFF);
            } else {
                let c = c * b;
                assert!(c.equals(&Fp2::ONE) == 0xFFFFFFFF);
            }

            if with_sqrt {
                let e = a * a;
                let (c, r) = e.sqrt();
                assert!(r == 0xFFFFFFFF);
                assert!((c * c).equals(&e) == 0xFFFFFFFF);
                let vc = c.encode();
                let zc0 = BigInt::from_bytes_le(Sign::Plus, &vc[..Fp::ENCODED_LENGTH]);
                let zc1 = BigInt::from_bytes_le(Sign::Plus, &vc[Fp::ENCODED_LENGTH..]);
                assert!(zc0.bit(0) == false);
                if zc0.sign() == Sign::NoSign {
                    assert!(zc1.bit(0) == false);
                }
                if a.iszero() == 0 {
                    assert!(e.legendre() == 1);
                    let e = a * a * Fp2::NQR;
                    assert!(e.legendre() == -1);
                    let (c, r) = e.sqrt();
                    assert!(r == 0);
                    assert!(c.iszero() == 0xFFFFFFFF);
                } else {
                    assert!(e.legendre() == 0);
                }

                if a0.iszero() == 0 {
                    let f = Fp2::new(&a0, &Fp::ZERO);
                    let (c, r) = f.sqrt();
                    assert!(r == 0xFFFFFFFF);
                    assert!((c * c).equals(&f) == 0xFFFFFFFF);
                    let g = -f;
                    let (c, r) = g.sqrt();
                    assert!(r == 0xFFFFFFFF);
                    assert!((c * c).equals(&g) == 0xFFFFFFFF);
                }
            }
        }

        #[test]
        fn fp2_ops() {
            let mut va = [0u8; (2 * Fp::ENCODED_LENGTH + 64) & !31usize];
            let mut vb = [0u8; (2 * Fp::ENCODED_LENGTH + 64) & !31usize];
            for i in 0..100 {
                let mut sh = Sha256::new();
                for j in 0..(va.len() >> 5) {
                    sh.update(((16 * i + 8 * j + 0) as u64).to_le_bytes());
                    va[(32 * j)..(32 * j + 32)].copy_from_slice(&sh.finalize_reset());
                }
                for j in 0..(vb.len() >> 5) {
                    sh.update(((16 * i + 8 * j + 1) as u64).to_le_bytes());
                    vb[(32 * j)..(32 * j + 32)].copy_from_slice(&sh.finalize_reset());
                }
                if i == 10 || i == 12 {
                    va.fill(0);
                }
                if i == 11 || i == 12 {
                    vb.fill(0);
                }
                check_fp2_ops(&va, &vb, true);
            }
        }
    };
} // End of macro: define_fp2_tests

#[cfg(test)]
pub(crate) use define_fp2_tests;
