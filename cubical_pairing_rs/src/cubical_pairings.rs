// Macro expections:
//    Fq    type of field element
macro_rules! define_cubical_core {
    ($Fq:ty) => {
        use std::fmt;

        /// Projective representation of a point (X : Y : Z)
        #[derive(Clone, Copy, Debug)]
        pub struct Point {
            X: $Fq,
            Y: $Fq,
            Z: $Fq,
        }

        impl Point {
            /// The point-at-infinity (neutral element of the group law).
            pub const INFINITY: Self = Self {
                X: Fq::ZERO,
                Y: Fq::ONE,
                Z: Fq::ZERO,
            };

            /// Create a new point: WARNING no check is made
            fn new(X: &$Fq, Y: &$Fq, Z: &$Fq) -> Self {
                Self {
                    X: *X,
                    Y: *Y,
                    Z: *Z,
                }
            }

            /// Negate the point
            pub fn set_neg(&mut self) {
                self.Y.set_neg()
            }

            /// Copy rhs into self if ctl == 0xFFFFFFFF.
            /// Do nothing is ctl == 0x00000000.
            /// ctl MUST be either 0xFFFFFFFF or 0x00000000.
            pub fn set_cond(&mut self, rhs: &Self, ctl: u32) {
                self.X.set_cond(&rhs.X, ctl);
                self.Y.set_cond(&rhs.Y, ctl);
                self.Z.set_cond(&rhs.Z, ctl);
            }

            // Return 0xFFFFFFFF if self and rhs represent the same point.
            /// Otherwise, return 0x00000000.
            pub fn equals(self, rhs: &Self) -> u32 {
                // P1 == P2 if and only if:
                //    P1 == inf AND P2 == inf
                //  OR:
                //    P1 != inf AND P2 != inf AND X1*Z2 = X2*Z1 AND Y1*Z2 = Y2*Z1
                let lz = self.Z.iszero();
                let rz = rhs.Z.iszero();
                let vx = (&self.X * &rhs.Z).equals(&(&rhs.X * &self.Z));
                let vy = (&self.Y * &rhs.Z).equals(&(&rhs.Y * &self.Z));
                (lz & rz) | (!lz & !rz & vx & vy)
            }
        }

        impl core::ops::Neg for Point {
            type Output = Point;

            #[inline(always)]
            fn neg(self) -> Point {
                let mut r = self;
                r.set_neg();
                r
            }
        }

        impl core::ops::Neg for &Point {
            type Output = Point;

            #[inline(always)]
            fn neg(self) -> Point {
                let mut r = *self;
                r.set_neg();
                r
            }
        }

        /// Special X-only representation of a point (or a pair of points,
        /// since two Y coordinates may match a given X).
        #[derive(Clone, Copy, Debug)]
        pub struct PointX {
            X: $Fq,
            Z: $Fq,
        }

        impl PointX {
            /// The neutral point of the group
            pub const INFINITY: Self = Self {
                X: <$Fq>::ONE,
                Z: <$Fq>::ZERO,
            };

            /// Create a point from coordinates
            pub fn new(X: &$Fq, Z: &Fq) -> Self {
                Self { X: *X, Z: *Z }
            }

            /// Recover the (X : Z) coordinates of a PointX
            fn coords(self) -> ($Fq, Fq) {
                (self.X, self.Z)
            }

            /// Point at infinity is of the form (X : 0)
            /// Returns 0xFFFFFFFF if Z == 0 and 0 otherwise
            pub fn is_zero(self) -> u32 {
                self.Z.iszero()
            }

            /// Returns the affine `x` coordinate of a point
            /// x = X / Z
            pub fn x(self) -> Fq {
                self.X / self.Z
            }

            /// Return 0xFFFFFFFF if self and rhs represent the same point.
            /// Otherwise, return 0x00000000.
            pub fn equals(self, rhs: &PointX) -> u32 {
                let inf1 = self.is_zero();
                let inf2 = rhs.is_zero();
                let e = (&self.X * &rhs.Z).equals(&(&rhs.X * &self.Z));
                (inf1 & inf2) | (!inf1 & !inf2 & e)
            }

            /// Affine translation by a two torsion point needed for even degree
            /// Tate pairings
            fn translate(self, T: Self) -> Self {
                let (A, B) = T.coords();
                // When we translates three things can happen.
                // - If T = (X : 0) then the translation of P is P
                // - If T = (0 : Z) with Z != 0, the translation of P = (X : Z)
                //   is (Z : X)
                // - Otherwise T = (A : B) and the translation of P = (X : Z) is
                // (A X - B Z : B X - A Z)

                // Compute generic values
                let mut X_new = &(&A * &self.X) - &(&B * &self.Z);
                let mut Z_new = &(&B * &self.X) - &(&A * &self.Z);

                // If A is zero, we should return (Z : X) instead
                let A_is_zero = A.iszero();
                X_new.set_cond(&self.Z, A_is_zero);
                Z_new.set_cond(&self.X, A_is_zero);

                // If B is zero, we should return (X : Z) instead
                let B_is_zero = B.iszero();
                X_new.set_cond(&self.X, B_is_zero);
                Z_new.set_cond(&self.Z, B_is_zero);

                Self::new(&X_new, &Z_new)
            }
        }

        /// Pretty print a point, used for debugging
        impl fmt::Display for PointX {
            fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
                write!(f, "Elliptic Curve PointX: ({} : {})", self.X, self.Z)
            }
        }

        /// Curve y^2 = x^3 + A*x^2 + x, for a given constant A
        /// (special case of a Montgomery curve).
        #[derive(Clone, Copy, Debug)]
        pub struct Curve {
            A: $Fq,   // A
            A24: $Fq, // (A+2)/4
        }

        impl Curve {
            /// Create a new curve instance, with the provided constant.
            pub fn new(A: &Fq) -> Self {
                // We check that the curve is not singular, i.e. A^2 != 4.
                assert!(A.equals(&<$Fq>::TWO) == 0);
                assert!((A + &<$Fq>::TWO).iszero() == 0);

                Self {
                    A: *A,
                    A24: (A + &<$Fq>::TWO).half().half(),
                }
            }

            /// x-only differential formula with cubical normalisation
            /// (XP : ZP), (XQ : ZQ), (1 : ixPQ) For PQ = P - Q
            #[inline(always)]
            fn cubical_xadd(xP: &$Fq, zP: &$Fq, xQ: &$Fq, zQ: &Fq, ixPQ: &$Fq) -> (Fq, Fq) {
                let V1 = &(xP - zP) * &(&*xQ + &*zQ);
                let V2 = &(xP + zP) * &(&*xQ - &*zQ);
                let X = ixPQ * &(&V1 + &V2).square();
                let Z = (&V1 - &V2).square();
                (X, Z)
            }

            /// x-only doubling and differential addition formula
            /// Assuming P - Q has been normalised.
            /// Note: order of arguments:
            /// (XP : ZP), (XQ : ZQ), (1: ixPQ) For PQ = P - Q
            /// Sets P = [2]P and Q = P + Q in place
            #[inline(always)]
            fn cubical_xdbladd(
                self,
                XP: &mut $Fq,
                ZP: &mut $Fq,
                XQ: &mut $Fq,
                ZQ: &mut $Fq,
                iXQP: &$Fq,
            ) {
                let mut t0 = &*XP + &*ZP;
                let mut t1 = &*XP - &*ZP;
                let mut X2P = t0.square();
                let mut t2 = &*XQ - &*ZQ;
                let mut XPQ = &*XQ + &*ZQ;
                t0 *= t2;
                let mut Z2P = t1.square();
                t1 *= &XPQ;
                t2 = &X2P - &Z2P;
                X2P *= &Z2P;
                XPQ = &self.A24 * &t2;
                let mut ZPQ = &t0 - &t1;
                Z2P = &XPQ + &Z2P;
                XPQ = &t0 + &t1;
                Z2P *= t2;
                ZPQ = ZPQ.square();
                XPQ = XPQ.square();
                XPQ *= iXQP;

                // Modify in place
                *XP = X2P;
                *ZP = Z2P;
                *XQ = XPQ;
                *ZQ = ZPQ;
            }

            /// Given xP, xQ xPQ with
            /// xP = x(P), xQ = x(Q) and xPQ = x(P - Q) sets:
            /// nP <- [n]P
            /// nPQ <- [n]P + Q
            /// Integer n is encoded as unsigned little-endian, with length
            /// nbitlen bits. Bits beyond that length are ignored.
            fn cubical_ladder_into(
                self,
                nP: &mut PointX,
                nPQ: &mut PointX,
                xP: &$Fq,
                xQ: &$Fq,
                xPQ: &$Fq,
                n: &[u8],
                nbitlen: usize,
                div_by_two: bool,
            ) {
                // Batch invert to obtain ixP, ixQ and ixPQ
                let mut inverses = [*xP, *xQ, *xPQ];
                <$Fq>::batch_invert(&mut inverses);
                let ixP = inverses[0];
                let mut ixQ = inverses[1];
                let mut ixPQ = inverses[2];

                // Set initial points for the ladder computation
                let mut xS0 = <$Fq>::ONE;
                let mut zS0 = <$Fq>::ZERO;

                let mut xS1 = *xP;
                let mut zS1 = <$Fq>::ONE;

                let mut xR;
                let mut zR;

                let mut xT = *xQ;
                let mut zT = <$Fq>::ONE;

                // If we want to divide by two, we skip the bottom bit in the ladder
                // this is used for the even pairings where the last factor of two is
                // computed using a translate.
                let start = div_by_two as usize;

                let mut cc = 0u32;
                for i in (start..nbitlen).rev() {
                    // First compute R = S0 + S1
                    (xR, zR) = Self::cubical_xadd(&xS0, &zS0, &xS1, &zS1, &ixP);

                    // Compute [2]Si and T + Si for i = {0, 1} depending on the bit of n
                    let ctl = (((n[i >> 3] >> (i & 7)) as u32) & 1).wrapping_neg();
                    <$Fq>::condswap(&mut xS0, &mut xS1, ctl ^ cc);
                    <$Fq>::condswap(&mut zS0, &mut zS1, ctl ^ cc);
                    <$Fq>::condswap(&mut ixQ, &mut ixPQ, ctl ^ cc);
                    self.cubical_xdbladd(&mut xS0, &mut zS0, &mut xT, &mut zT, &ixQ);

                    // Update Sj to be R
                    xS1 = xR;
                    zS1 = zR;

                    // Save swap variable
                    cc = ctl;
                }
                // Perform a final swap
                <$Fq>::condswap(&mut xS0, &mut xS1, cc);
                <$Fq>::condswap(&mut zS0, &mut zS1, cc);

                // Set the variables to return S0, T = [n]P, [n]P + Q
                nP.X = xS0;
                nP.Z = zS0;
                nPQ.X = xT;
                nPQ.Z = zT;
            }

            /// Given x(P), x(Q) and x(P - Q) computes [n]P and [n]P + Q for Tate pairing.
            /// Integer n is encoded as unsigned little-endian, with length
            /// nbitlen bits. Bits beyond that length are ignored.
            pub fn cubical_ladder(
                self,
                xP: &$Fq,
                xQ: &$Fq,
                xPQ: &$Fq,
                n: &[u8],
                nbitlen: usize,
                div_by_two: bool,
            ) -> (PointX, PointX) {
                let mut nP = PointX::INFINITY;
                let mut nPQ = PointX::INFINITY;
                self.cubical_ladder_into(&mut nP, &mut nPQ, xP, xQ, xPQ, n, nbitlen, div_by_two);
                (nP, nPQ)
            }

            /// Given xP, xQ xPQ with
            /// xP = x(P), xQ = x(Q) and xPQ = x(P - Q) sets:
            /// nP <- [2^e]P
            /// nPQ <- [2^e]P + Q
            fn cubical_ladder_2exp_into(
                self,
                nP: &mut PointX,
                nPQ: &mut PointX,
                xP: &$Fq,
                xQ: &$Fq,
                xPQ: &$Fq,
                e: usize,
            ) {
                // Batch invert to obtain ixP, ixQ and ixPQ
                let mut inverses = [*xQ, *xPQ];
                <$Fq>::batch_invert(&mut inverses);
                let ixQ = inverses[0];
                let ixPQ = inverses[1];

                let mut X0 = *xP;
                let mut Z0 = <$Fq>::ONE;
                let mut X1 = *xQ;
                let mut Z1 = <$Fq>::ONE;

                // Compute P + Q
                (X1, Z1) = Self::cubical_xadd(&X0, &Z0, &X1, &Z1, &ixPQ);

                // Compute [2^e]P and [2^e]P + Q
                for _ in 0..e {
                    self.cubical_xdbladd(&mut X0, &mut Z0, &mut X1, &mut Z1, &ixQ);
                }

                // Set the variables of [2^e]P, [2^e]P + Q
                nP.X = X0;
                nP.Z = Z0;
                nPQ.X = X1;
                nPQ.Z = Z1;
            }

            /// Given x(P), x(Q) and x(P - Q) computes the points
            /// [2^e]P and [2^e]P + Q for the Tate pairing.
            pub fn cubical_ladder_2exp(
                self,
                xP: &$Fq,
                xQ: &$Fq,
                xPQ: &$Fq,
                e: usize,
            ) -> (PointX, PointX) {
                let mut nP = PointX::INFINITY;
                let mut nPQ = PointX::INFINITY;
                self.cubical_ladder_2exp_into(&mut nP, &mut nPQ, xP, xQ, xPQ, e);
                (nP, nPQ)
            }

            /// Given x(P), x(Q) and x(P - Q) computes the non-reduced tate pairing
            fn tate_pairing_not_reduced(
                self,
                xP: &$Fq,
                xQ: &$Fq,
                xPQ: &$Fq,
                n: &[u8],
                nbitlen: usize,
            ) -> (Fq, Fq) {
                let is_even: bool = n[0] & 1 == 0;
                let (mut nP, mut nPQ) = self.cubical_ladder(xP, xQ, xPQ, n, nbitlen, is_even);
                if is_even {
                    nPQ = nPQ.translate(nP);
                    nP = nP.translate(nP);
                }
                (nPQ.Z, nP.X)
            }

            /// Given x(P), x(Q) and x(P - Q) computes the non-reduced tate
            /// pairing using points of order 2^e
            fn tate_pairing_not_reduced_2exp(
                self,
                xP: &$Fq,
                xQ: &$Fq,
                xPQ: &$Fq,
                e: usize,
            ) -> (Fq, Fq) {
                let (mut nP, mut nPQ) = self.cubical_ladder_2exp(xP, xQ, xPQ, e - 1);
                nPQ = nPQ.translate(nP);
                nP = nP.translate(nP);
                (nPQ.Z, nP.X)
            }

            /// Compute x^(p^2 - 1) = x^(p + 1)(p - 1) to reduce the tate pairing for E / Fp^2
            /// where Fp^2 has modulus x^2 + 1
            fn reduce_tate_pairing(self, num: &$Fq, den: &$Fq, d: &[u8], dbitlen: usize) -> Fq {
                // We can compute the pth power with a conjugate \pi(x) = x^p
                let num_p = num.conjugate();
                let den_p = den.conjugate();

                // First we compute ePQ^(p-1) using Frobenius x^p / x = x^(p-1)
                let ePQ = (num_p * den) / (den_p * num);

                // Now handle the power of (p + 1) / d with supplied exponent
                ePQ.pow(d, dbitlen)
            }

            /// Given x(P), x(Q) and x(P + Q) computes the reduced tate
            /// pairing
            /// d encodes in little endian the value (p + 1) // 2^e
            pub fn tate_pairing(
                self,
                xP: &$Fq,
                xQ: &$Fq,
                xPQ: &$Fq,
                n: &[u8],
                nbitlen: usize,
                d: &[u8],
                dbitlen: usize,
            ) -> Fq {
                // First compute the non-reduced tate pairing which is returned as num / den
                let (num, den) = self.tate_pairing_not_reduced(xP, xQ, xPQ, n, nbitlen);

                // Efficiently compute ePQ^(p^2 - 1)
                self.reduce_tate_pairing(&num, &den, d, dbitlen)
            }

            /// Given x(P), x(Q) and x(P - Q) computes the reduced tate
            /// pairing using points of order 2^e
            /// d encodes in little endian the value (p + 1) // 2^e
            pub fn tate_pairing_2exp(
                self,
                xP: &$Fq,
                xQ: &$Fq,
                xPQ: &$Fq,
                e: usize,
                d: &[u8],
                dbitlen: usize,
            ) -> Fq {
                // First compute the non-reduced tate pairing
                let (num, den) = self.tate_pairing_not_reduced_2exp(xP, xQ, xPQ, e);

                // Efficiently compute ePQ^(p^2 - 1)
                self.reduce_tate_pairing(&num, &den, d, dbitlen)
            }

            /// Given x(P), x(Q) and x(P - Q) computes the Weil pairing
            pub fn weil_pairing(
                self,
                xP: &$Fq,
                xQ: &$Fq,
                xPQ: &$Fq,
                n: &[u8],
                nbitlen: usize,
            ) -> Fq {
                let (e1_num, e1_den) = self.tate_pairing_not_reduced(xP, xQ, xPQ, n, nbitlen);
                let (e2_num, e2_den) = self.tate_pairing_not_reduced(xQ, xP, xPQ, n, nbitlen);
                (e1_num * e2_den) / (e2_num * e1_den)
            }

            /// Given x(P), x(Q) and x(P - Q) computes the weil pairing
            /// using points of order = 2^e
            pub fn weil_pairing_2exp(self, xP: &$Fq, xQ: &$Fq, xPQ: &$Fq, e: usize) -> Fq {
                let (e1_num, e1_den) = self.tate_pairing_not_reduced_2exp(xP, xQ, xPQ, e);
                let (e2_num, e2_den) = self.tate_pairing_not_reduced_2exp(xQ, xP, xPQ, e);
                (e1_num * e2_den) / (e2_num * e1_den)
            }

            //
            // The following code is used to do the full point compression as used in SQIsign. For those
            // interested in the pairing computations only, everything below this comment can be ignored.
            //

            /// P3 <- P1 + P2
            fn add_into(self, P3: &mut Point, P1: &Point, P2: &Point) {
                // Complete routine, to handle all edge cases:
                //   if Z1 == 0:            # P1 == inf
                //       return P2
                //   if Z2 == 0:            # P2 == inf
                //       return P1
                //   L <- Y2*Z1 - Y1*Z2
                //   T <- X2*Z1 - X1*Z2
                //   if T == 0:             # x1 == x2
                //       if L == 0:         # ... and y1 == y2: doubling case
                //           L <- 3*X1^2 + 2*A*X1*Z1 + Z1^2
                //           T <- 2*Y1*Z1
                //       else:              # ... but y1 != y2, thus P2 = -P1
                //           return inf
                //   U <- Z1*Z2*L^2 - (X1*Z2 + X2*Z1 + A*Z1*Z2)*T^2
                //   X3 <- U*T
                //   Y3 <- L*(X1*Z2*T^2 - U) - Y1*Z2*T^3
                //   Z3 <- Z1*Z2*T^3
                //
                // Constant-time processing:
                //   Cases P1 == inf and P2 == inf are handled at the end.
                //   (L,T) are always computed for both normal and doubling cases.
                //   If P1 == -P2 then we can let T == 0 and L != 0, this will
                //   properly lead to Z3 == 0.
                //
                // Formulas from https://eprint.iacr.org/2015/1060 are faster
                // but do not cover the case when P1 - P2 is a point of order 2,
                // which can happen in all generality.
                //
                // : current formulas have cost 16M+5S; this can probably
                // be improved. Main issues to tackle:
                //   - Multiplications by A are expensive (since A can be any value)
                //   - There are three points of order 2; this makes finding
                //     complete formulas challenging.

                // T = X2*Z1 - X1*Z2
                // L = Y2*Z1 - Y1*Z2
                let x1z2 = &P1.X * &P2.Z;
                let x2z1 = &P2.X * &P1.Z;
                let mut T = &x2z1 - &x1z2;
                let y1z2 = &P1.Y * &P2.Z;
                let y2z1 = &P2.Y * &P1.Z;
                let mut L = &y2z1 - &y1z2;

                // Alternate (T,L) for doubling:
                //   Td = 2*Y1*Z1
                //   Ld = 3*X1^2 + 2*A*X1*Z1 + Z1^2
                let dbl = T.iszero() & L.iszero();
                let Td = (&P1.Y * &P1.Z).mul2();
                let x1x1 = P1.X.square();
                let z1z1 = P1.Z.square();
                let dx1z1 = &(&P1.X + &P1.Z).square() - &x1x1 - &z1z1;
                let Ld = &x1x1.mul3() + &z1z1 + &self.A * &dx1z1;
                T.set_cond(&Td, dbl);
                L.set_cond(&Ld, dbl);

                // U = L^2*Z1*Z2 - (X1*Z2 + X2*Z1 + A*Z1*Z2)*T^2
                let T2 = T.square();
                let T3 = &T * &T2;
                let z1z2 = &P1.Z * &P2.Z;
                let U = &(&L.square() * &z1z2) - &(&(&x1z2 + &x2z1 + &(&self.A * &z1z2)) * &T2);

                // X3 = U*T
                // Y3 = L*(X1*Z2*T^2 - U) - Y1*Z2*T^3
                // Z3 = Z1*Z2*T^3
                P3.X = &U * &T;
                P3.Y = &(&L * &(&(&x1z2 * &T2) - &U)) - &(&y1z2 * &T3);
                P3.Z = &z1z2 * &T3;

                // Corrective action in case one of the inputs was the
                // point-at-infinity.
                let inf1 = P1.Z.iszero();
                let inf2 = P2.Z.iszero();
                P3.set_cond(&P2, inf1);
                P3.set_cond(&P1, inf2);
            }

            /// Return P1 + P2 as a new point
            fn add(self, P1: &Point, P2: &Point) -> Point {
                let mut P3 = Point::INFINITY;
                self.add_into(&mut P3, P1, P2);
                P3
            }

            /// Given the x-coordinate of a point, lift it to a projective point
            fn lift_point(self, x: &$Fq) -> Point {
                let mut y = x + self.A; // y = x + A
                y *= x; // y = x^2 + A*x
                y += <$Fq>::ONE; // x^2 + A*x + 1
                y *= x; // y = x^3 + A*x^2 + x
                y.set_sqrt();
                Point::new(&x, &y, &<$Fq>::ONE)
            }

            /// Given the x-coordinates of x(P), x(Q) and x(P - Q) lift the points
            /// onto the curve <P, Q>
            fn lift_basis(self, xP: &$Fq, xQ: &$Fq, xPQ: &$Fq) -> (Point, Point) {
                let P = self.lift_point(xP);

                // Okeya-Sakurai algorithm to recover Q.Y without a sqrt
                let mut v2 = xP + xQ;
                let mut v3 = xP - xQ;
                v3.set_square();
                v3 *= xPQ;
                let mut v1 = self.A.mul2();
                v2 += v1;
                let mut v4 = xP * xQ;
                v4 += <$Fq>::ONE;
                v2 *= v4;
                v2 -= v1;
                let y = v3 - v2;
                v1 = P.Y + P.Y;
                let x = xQ * v1;
                let Q = Point::new(&x, &y, &v1);

                return (P, Q);
            }

            /// Given the x-coordinates of two bases, compute pairs of differences
            /// x(R - P), x(R - Q), x(S - P), x(S - Q)
            fn compute_difference_points(
                self,
                xP: &$Fq,
                xQ: &$Fq,
                xPQ: &$Fq,
                xR: &$Fq,
                xS: &$Fq,
                xRS: &$Fq,
            ) -> ($Fq, $Fq, $Fq, $Fq) {
                // Lift x-coordinates to projective points on curve
                let (P, Q) = self.lift_basis(xP, xQ, xPQ);
                let (R, S) = self.lift_basis(xR, xS, xRS);

                // Compute R - P, R - Q, S - P, S - Q
                let RmP = self.add(&R, &P);
                let RmQ = self.add(&R, &Q);
                let SmP = self.add(&S, &P);
                let SmQ = self.add(&S, &Q);

                // Invert the four Z coordinates of the differences
                let mut zs = [RmP.Z, RmQ.Z, SmP.Z, SmQ.Z];
                <$Fq>::batch_invert(&mut zs);

                // Compute the normalized x-coordinate of the differences
                (RmP.X * zs[0], RmQ.X * zs[1], SmP.X * zs[2], SmQ.X * zs[3])
            }

            /// Given x(P), x(Q) and x(P - Q) for a basis <P, Q> of E[2^f] and
            /// x(R), x(S) and x(R - S) for a basis <P, Q> of E[2^e] for e <= f
            /// compute integers a, b, c, d such that:
            /// R = 2^(e - f) ([a]P + [b]B)
            /// S = 2^(e - f) ([c]P + [d]Q)
            /// represented as little endian values.
            pub fn point_compression(
                self,
                xP: &$Fq,
                xQ: &$Fq,
                xPQ: &$Fq,
                xR: &$Fq,
                xS: &$Fq,
                xRS: &$Fq,
                f: usize,
                e: usize,
                d: &[u8],
                dbitlen: usize,
                dlog_table: &[usize],
            ) -> (Vec<u8>, Vec<u8>, Vec<u8>, Vec<u8>, u32) {
                let e_diff = f - e;

                // We are going to compute five non-reduced tate pairings:
                // e(P, Q), e(R, P), e(R, Q), e(S, P) and e(S, Q)
                // To do this, we first need to precompute the points we need for the five
                // cubical ladders, x(R - P), x(R - Q), x(S - P), x(S - Q)
                let (xRP, xRQ, xSP, xSQ) = self.compute_difference_points(xP, xQ, xPQ, xR, xS, xRS);

                // Now we need to compute the inverses for the ladders below
                let mut inv = [*xP, *xQ];
                <$Fq>::batch_invert(&mut inv);
                let ixP = inv[0];
                let ixQ = inv[1];

                // We compute the pairing e(P, Q) by computing [2^f]P and [2^f]P - Q which for our case
                // is the same as e(P, -Q) = 1/e(P, Q) = e(Q, P)
                let mut nP = PointX::new(xP, &<$Fq>::ONE);
                let mut nPQ = PointX::new(xPQ, &<$Fq>::ONE);
                for _ in 0..(f - 1) {
                    self.cubical_xdbladd(&mut nP.X, &mut nP.Z, &mut nPQ.X, &mut nPQ.Z, &ixQ);
                }
                nPQ = nPQ.translate(nP);
                nP = nP.translate(nP);

                // For the remaining four pairings we want to compute pairings of order e
                // which requires us to compute the values
                // [2^e]R, [2^e]S, [2^e]R - P, [2^e]R - Q, [2^e]S - P, [2^e]S - Q
                let mut nR = PointX::new(xR, &<$Fq>::ONE);
                let mut nS = PointX::new(xS, &<$Fq>::ONE);
                let mut nRP = PointX::new(&xRP, &<$Fq>::ONE);
                let mut nRQ = PointX::new(&xRQ, &<$Fq>::ONE);
                let mut nSP = PointX::new(&xSP, &<$Fq>::ONE);
                let mut nSQ = PointX::new(&xSQ, &<$Fq>::ONE);
                for _ in 0..(e - 1) {
                    (nRP.X, nRP.Z) = Self::cubical_xadd(&nRP.X, &nRP.Z, &nR.X, &nR.Z, &ixP);
                    self.cubical_xdbladd(&mut nR.X, &mut nR.Z, &mut nRQ.X, &mut nRQ.Z, &ixQ);

                    (nSP.X, nSP.Z) = Self::cubical_xadd(&nSP.X, &nSP.Z, &nS.X, &nS.Z, &ixP);
                    self.cubical_xdbladd(&mut nS.X, &mut nS.Z, &mut nSQ.X, &mut nSQ.Z, &ixQ);
                }
                nRP = nRP.translate(nR);
                nRQ = nRQ.translate(nR);
                nSP = nSP.translate(nS);
                nSQ = nSQ.translate(nS);
                nR = nR.translate(nR);
                nS = nS.translate(nS);

                // Represent the pairings projectively
                let num = [nPQ.Z, nRP.Z, nRQ.Z, nSP.Z, nSQ.Z];
                let den = [nP.X, nR.X, nR.X, nS.X, nS.X];

                // We now raise each of these to the power p^2 - 1
                let mut pairings: [$Fq; 5] = [<$Fq>::ONE; 5];
                let mut pairings_den: [$Fq; 5] = [<$Fq>::ONE; 5];

                // First raise num and den to power of p-1
                for i in 0..5 {
                    pairings[i] = num[i].conjugate() * den[i];
                    pairings_den[i] = den[i].conjugate() * num[i];
                }
                // Invert the denominators
                <$Fq>::batch_invert(&mut pairings_den);

                // Now compute power of (p + 1 // 2^f) * 2^(f - e)
                for i in 0..5 {
                    pairings[i] = (pairings[i] * pairings_den[i]).pow(d, dbitlen);
                    for _ in 0..e_diff {
                        pairings[i] = pairings[i].square();
                    }
                }

                let (gpp, ok0) = pairings[0].precompute_dlp_table(dlog_table);
                let (r2, ok1) = pairings[0].solve_dlp_2e(&pairings[1], e, &gpp, dlog_table);
                let (r1, ok2) = pairings[0].solve_dlp_2e(&pairings[2], e, &gpp, dlog_table);
                let (s2, ok3) = pairings[0].solve_dlp_2e(&pairings[3], e, &gpp, dlog_table);
                let (s1, ok4) = pairings[0].solve_dlp_2e(&pairings[4], e, &gpp, dlog_table);

                (r1, r2, s1, s2, ok0 & ok1 & ok2 & ok3 & ok4)
            }
        }
    };
} // End of macro: define_cubical_core

pub(crate) use define_cubical_core;
