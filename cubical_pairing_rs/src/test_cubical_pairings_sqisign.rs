#![allow(non_snake_case)]
#![allow(non_upper_case_globals)]

macro_rules! define_cubical_test_even_only {
    ($Curve:ty, $Fq:ty, $PointX:ty) => {
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
            let Q = <$PointX>::new(&xQ, &<$Fq>::ONE);

            // Difference point x(P - Q)
            let (xPQ, _) = <$Fq>::decode(&hex::decode(PmQX_str).unwrap());

            let (nP_test, nPQ_test) = E0.cubical_ladder_2exp(&xP, &xQ, &xPQ, ea);

            assert!(nP_test.equals(&<$PointX>::INFINITY) == u32::MAX);
            assert!(nPQ_test.equals(&Q) == u32::MAX);
        }

        #[test]
        fn test_tate_pairing_2exp() {
            // Create E0
            // A supersingular elliptic curve of order f * 2^ea
            let (A0, _) = <$Fq>::decode(&hex::decode(A0_str).unwrap());
            let E0 = <$Curve>::new(&A0);

            // Point of order 2^ea x(P)
            let (xP, _) = <$Fq>::decode(&hex::decode(PX_str).unwrap());

            // Point of order 2^ea x(Q)
            let (xQ, _) = <$Fq>::decode(&hex::decode(QX_str).unwrap());

            // Difference point x(P - Q)
            let (xPQ, _) = <$Fq>::decode(&hex::decode(PmQX_str).unwrap());

            // Tate Pairing computed via SageMath
            let (eQ2P2, _) = <$Fq>::decode(&hex::decode(tate_pairing_str).unwrap());

            // Compute Tate pairing
            let eQ2P2_test = E0.tate_pairing_2exp(&xP, &xQ, &xPQ, ea, &dA, dA_BITLEN);

            assert!(eQ2P2.equals(&eQ2P2_test) == u32::MAX);
        }

        #[test]
        fn test_weil_pairing_2exp() {
            // Create E0
            // A supersingular elliptic curve of order 2^74 * 3^41
            let (A0, _) = <$Fq>::decode(&hex::decode(A0_str).unwrap());
            let E0 = <$Curve>::new(&A0);

            // Point of order 2^ea x(P)
            let (xP, _) = <$Fq>::decode(&hex::decode(PX_str).unwrap());

            // Point of order 2^ea x(Q)
            let (xQ, _) = <$Fq>::decode(&hex::decode(QX_str).unwrap());

            // Difference point x(P - Q)
            let (xPQ, _) = <$Fq>::decode(&hex::decode(PmQX_str).unwrap());

            // Weil Pairing computed via SageMath
            let (eQ2P2, _) = <$Fq>::decode(&hex::decode(weil_pairing_str).unwrap());

            // Compute Weil pairing
            let eQ2P2_test = E0.weil_pairing_2exp(&xP, &xQ, &xPQ, ea);

            assert!(eQ2P2.equals(&eQ2P2_test) == u32::MAX);
        }
    };
}

#[cfg(test)]
mod test_level_one {
    // Common imports for both test modules
    use crate::test_data::sqisign_one::*;

    // Import for the first set of tests
    use crate::kummer248::{Curve, Fq, PointX};

    // Define the tests for the first set
    define_cubical_test_even_only!(Curve, Fq, PointX);
}

#[cfg(test)]
mod test_level_three {
    // Common imports for both test modules
    use crate::test_data::sqisign_three::*;

    // Import for the first set of tests
    use crate::kummer383::{Curve, Fq, PointX};

    // Define the tests for the first set
    define_cubical_test_even_only!(Curve, Fq, PointX);
}

#[cfg(test)]
mod test_level_five {
    // Common imports for both test modules
    use crate::test_data::sqisign_five::*;

    // Import for the first set of tests
    use crate::kummer505::{Curve, Fq, PointX};

    // Define the tests for the first set
    define_cubical_test_even_only!(Curve, Fq, PointX);
}
