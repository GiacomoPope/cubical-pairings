#![allow(non_snake_case)]
#![allow(non_upper_case_globals)]

#[cfg(test)]
macro_rules! define_cubical_test {
    ($Curve:ty, $Fq:ty, $PointX:ty) => {
        #[test]
        fn test_cubical_ladder() {
            // Create E0
            // A supersingular elliptic curve of order 2^ea * 3^eb
            let (A0, _) = <$Fq>::decode(&hex::decode(A0_str).unwrap());
            let E0 = <$Curve>::new(&A0);

            // Generator x(P)
            let (xP, _) = <$Fq>::decode(&hex::decode(PX_str).unwrap());

            // Generator x(Q)
            let (xQ, _) = <$Fq>::decode(&hex::decode(QX_str).unwrap());

            // Difference point x(P - Q)
            let (xPQ, _) = <$Fq>::decode(&hex::decode(PmQX_str).unwrap());

            // P2 = [B] P (Even torsion part)
            let (xP2, _) = <$Fq>::decode(&hex::decode(P2X_str).unwrap());
            let P2 = <$PointX>::new(&xP2, &<$Fq>::ONE);

            // [B]P + Q
            let (BP_Q_X, _) = <$Fq>::decode(&hex::decode(BP_QX_str).unwrap());
            let BP_Q = <$PointX>::new(&BP_Q_X, &<$Fq>::ONE);

            // Compute [n]P and [n]P + Q
            let (nP_test, nPQ_test) = E0.cubical_ladder(&xP, &xQ, &xPQ, &B, B_BITLEN, false);

            assert!(P2.equals(&nP_test) == u32::MAX);
            assert!(BP_Q.equals(&nPQ_test) == u32::MAX);
        }

        #[test]
        fn test_cubical_ladder_2exp() {
            // Create E0
            // A supersingular elliptic curve of order 2^ea * 3^eb
            let (A0, _) = <$Fq>::decode(&hex::decode(A0_str).unwrap());
            let E0 = <$Curve>::new(&A0);

            // Generator x(P)
            let (xP, _) = <$Fq>::decode(&hex::decode(PX_str).unwrap());

            // Generator x(Q)
            let (xQ, _) = <$Fq>::decode(&hex::decode(QX_str).unwrap());

            // Difference point x(P - Q)
            let (xPQ, _) = <$Fq>::decode(&hex::decode(PmQX_str).unwrap());

            // P3 = [2^ea] P (Odd torsion part)
            let (xP3, _) = <$Fq>::decode(&hex::decode(P3X_str).unwrap());
            let P3 = <$PointX>::new(&xP3, &<$Fq>::ONE);

            // [2^ea] P + Q
            let (AP_Q_X, _) = <$Fq>::decode(&hex::decode(AP_QX_str).unwrap());
            let AP_Q = <$PointX>::new(&AP_Q_X, &<$Fq>::ONE);

            let (nP_test, nPQ_test) = E0.cubical_ladder_2exp(&xP, &xQ, &xPQ, ea);

            assert!(P3.equals(&nP_test) == u32::MAX);
            assert!(AP_Q.equals(&nPQ_test) == u32::MAX);
        }

        #[test]
        fn test_tate_pairing_2exp() {
            // Create E0
            // A supersingular elliptic curve of order 2^ea * 3^eb
            let (A0, _) = <$Fq>::decode(&hex::decode(A0_str).unwrap());
            let E0 = <$Curve>::new(&A0);

            // x-coordinate x(P) of order 2^ea
            let (xP, _) = <$Fq>::decode(&hex::decode(P2X_str).unwrap());

            // x-coordinate x(Q) of order 2^ea
            let (xQ, _) = <$Fq>::decode(&hex::decode(Q2X_str).unwrap());

            // x-coordinate x(P - Q) of order 2^ea
            let (xPQ, _) = <$Fq>::decode(&hex::decode(P2mQ2X_str).unwrap());

            // Tate Pairing computed via SageMath
            let (eQ2P2, _) = <$Fq>::decode(&hex::decode(tate_pairing_A_str).unwrap());

            // Compute Tate pairing
            let eQ2P2_test = E0.tate_pairing_2exp(&xP, &xQ, &xPQ, ea, &dA, dA_BITLEN);

            assert!(eQ2P2.equals(&eQ2P2_test) == u32::MAX);
        }

        #[test]
        fn test_weil_pairing_2exp() {
            // Create E0
            // A supersingular elliptic curve of order 2^ea * 3^eb
            let (A0, _) = <$Fq>::decode(&hex::decode(A0_str).unwrap());
            let E0 = <$Curve>::new(&A0);

            // x-coordinate x(P) of order 2^ea
            let (xP, _) = <$Fq>::decode(&hex::decode(P2X_str).unwrap());

            // x-coordinate x(Q) of order 2^ea
            let (xQ, _) = <$Fq>::decode(&hex::decode(Q2X_str).unwrap());

            // x-coordinate x(P - Q) of order 2^ea
            let (xPQ, _) = <$Fq>::decode(&hex::decode(P2mQ2X_str).unwrap());

            // Weil Pairing computed via SageMath
            let (eQ2P2, _) = <$Fq>::decode(&hex::decode(weil_pairing_A_str).unwrap());

            // Compute Weil pairing
            let eQ2P2_test = E0.weil_pairing_2exp(&xP, &xQ, &xPQ, ea);

            assert!(eQ2P2.equals(&eQ2P2_test) == u32::MAX);
        }

        #[test]
        fn test_tate_pairing_odd() {
            // Create E0
            // A supersingular elliptic curve of order 2^ea * 3^eb
            let (A0, _) = <$Fq>::decode(&hex::decode(A0_str).unwrap());
            let E0 = <$Curve>::new(&A0);

            // x-coordinate x(P) of order 3^eb
            let (xP, _) = <$Fq>::decode(&hex::decode(P3X_str).unwrap());

            // x-coordinate x(Q) of order 3^eb
            let (xQ, _) = <$Fq>::decode(&hex::decode(Q3X_str).unwrap());

            // x-coordinate x(P-Q) of order 3^eb
            let (xPQ, _) = <$Fq>::decode(&hex::decode(P3mQ3X_str).unwrap());

            // Tate Pairing computed via SageMath
            let (eQ3P3, _) = <$Fq>::decode(&hex::decode(tate_pairing_B_str).unwrap());

            // Compute Tate pairing
            let eQ3P3_test = E0.tate_pairing(&xP, &xQ, &xPQ, &B, B_BITLEN, &dB, dB_BITLEN);

            assert!(eQ3P3.equals(&eQ3P3_test) == u32::MAX);
        }

        #[test]
        fn test_weil_pairing_odd() {
            // Create E0
            // A supersingular elliptic curve of order 2^ea * 3^eb
            let (A0, _) = <$Fq>::decode(&hex::decode(A0_str).unwrap());
            let E0 = <$Curve>::new(&A0);

            // x-coordinate x(P) of order 3^eb
            let (xP, _) = <$Fq>::decode(&hex::decode(P3X_str).unwrap());

            // x-coordinate x(Q) of order 3^eb
            let (xQ, _) = <$Fq>::decode(&hex::decode(Q3X_str).unwrap());

            // x-coordinate x(P-Q) of order 3^eb
            let (xPQ, _) = <$Fq>::decode(&hex::decode(P3mQ3X_str).unwrap());

            // Weil Pairing computed via SageMath
            let (eQ3P3, _) = <$Fq>::decode(&hex::decode(weil_pairing_B_str).unwrap());

            // Compute Weil pairing
            let eQ3P3_test = E0.weil_pairing(&xP, &xQ, &xPQ, &B, B_BITLEN);

            assert!(eQ3P3.equals(&eQ3P3_test) == u32::MAX);
        }

        #[test]
        fn test_tate_pairing_even() {
            // Create E0
            // A supersingular elliptic curve of order 2^ea * 3^eb
            let (A0, _) = <$Fq>::decode(&hex::decode(A0_str).unwrap());
            let E0 = <$Curve>::new(&A0);

            // Generator x(P)
            let (xP, _) = <$Fq>::decode(&hex::decode(PX_str).unwrap());

            // Generator x(Q)
            let (xQ, _) = <$Fq>::decode(&hex::decode(QX_str).unwrap());

            // Difference point x(P - Q)
            let (xPQ, _) = <$Fq>::decode(&hex::decode(PmQX_str).unwrap());

            // Tate Pairing computed via SageMath
            let (eQP, _) = <$Fq>::decode(&hex::decode(tate_pairing_str).unwrap());

            // Compute Tate pairing
            let eQP_test = E0.tate_pairing(&xP, &xQ, &xPQ, &N, N_BITLEN, &dN, dN_BITLEN);

            assert!(eQP.equals(&eQP_test) == u32::MAX);
        }

        #[test]
        fn test_weil_pairing_even() {
            // Create E0
            // A supersingular elliptic curve of order 2^ea * 3^eb
            let (A0, _) = <$Fq>::decode(&hex::decode(A0_str).unwrap());
            let E0 = <$Curve>::new(&A0);

            // Generator x(P)
            let (xP, _) = <$Fq>::decode(&hex::decode(PX_str).unwrap());

            // Generator x(Q)
            let (xQ, _) = <$Fq>::decode(&hex::decode(QX_str).unwrap());

            // Difference point x(P - Q)
            let (xPQ, _) = <$Fq>::decode(&hex::decode(PmQX_str).unwrap());

            // Weil Pairing computed via SageMath
            let (eQP, _) = <$Fq>::decode(&hex::decode(weil_pairing_str).unwrap());

            // Compute Tate pairing
            let eQP_test = E0.weil_pairing(&xP, &xQ, &xPQ, &N, N_BITLEN);

            assert!(eQP.equals(&eQP_test) == u32::MAX);
        }
    };
}

#[cfg(test)]
mod test_level_one {
    // Run tests with sml values
    use crate::kummer434::{Curve, Fq, PointX};
    use crate::test_data::sidh_one::*;

    define_cubical_test!(Curve, Fq, PointX);
}

#[cfg(test)]
mod test_level_three {
    // Run tests with sml values
    use crate::kummer610::{Curve, Fq, PointX};
    use crate::test_data::sidh_three::*;

    define_cubical_test!(Curve, Fq, PointX);
}

#[cfg(test)]
mod test_level_five {
    // Run tests with sml values
    use crate::kummer751::{Curve, Fq, PointX};
    use crate::test_data::sidh_five::*;

    define_cubical_test!(Curve, Fq, PointX);
}
