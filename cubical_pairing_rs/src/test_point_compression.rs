#![allow(non_snake_case)]
#![allow(non_upper_case_globals)]

macro_rules! define_compression_test {
    ($Curve:ty, $Fq:ty) => {
        #[test]
        fn test_point_compression() {
            // Create E0
            // A supersingular elliptic curve of order 2^ea * 3^eb
            let (A0, _) = <$Fq>::decode(&hex::decode(A0_str).unwrap());
            let E0 = <$Curve>::new(&A0);

            // Generator x(P), x(Q), x(P - Q)
            let (xP, _) = <$Fq>::decode(&hex::decode(PX_str).unwrap());
            let (xQ, _) = <$Fq>::decode(&hex::decode(QX_str).unwrap());
            let (xPQ, _) = <$Fq>::decode(&hex::decode(PmQX_str).unwrap());

            // Generator x(R), x(S), x(R - S)
            let (xR, _) = <$Fq>::decode(&hex::decode(RX_str).unwrap());
            let (xS, _) = <$Fq>::decode(&hex::decode(SX_str).unwrap());
            let (xRS, _) = <$Fq>::decode(&hex::decode(RmSX_str).unwrap());

            let (r1, r2, s1, s2, check) = E0.point_compression(
                &xP,
                &xQ,
                &xPQ,
                &xR,
                &xS,
                &xRS,
                f,
                e,
                &dA,
                dA_BITLEN,
                &dlog_table,
            );

            assert!(r1 == a);
            assert!(r2 == b);
            assert!(s1 == c);
            assert!(s2 == d);
            assert!(check == u32::MAX);
        }
    };
}

#[cfg(test)]
mod test_level_one {
    // Common imports for both test modules
    use crate::test_data::compression_one::*;

    // Import for the first set of tests
    use crate::kummer248::{Curve, Fq};

    // Define the tests for the first set
    define_compression_test!(Curve, Fq);
}

#[cfg(test)]
mod test_level_three {
    // Common imports for both test modules
    use crate::test_data::compression_three::*;

    // Import for the first set of tests
    use crate::kummer383::{Curve, Fq};

    // Define the tests for the first set
    define_compression_test!(Curve, Fq);
}

#[cfg(test)]
mod test_level_five {
    // Common imports for both test modules
    use crate::test_data::compression_five::*;

    // Import for the first set of tests
    use crate::kummer505::{Curve, Fq};

    // Define the tests for the first set
    define_compression_test!(Curve, Fq);
}
