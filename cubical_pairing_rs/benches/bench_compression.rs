#![allow(non_snake_case)]
#![allow(non_upper_case_globals)]
#![allow(non_snake_case)]
#![allow(non_upper_case_globals)]
#![allow(unused_macros)] // macro used for benchmarking upsets the linter
#![allow(unused_imports)]

macro_rules! define_compression_bench {
    ($Curve:ty, $Fq:ty) => {
        fn bench_compression(cc: &mut Criterion) {
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

            // Assert we're benchmarking something which works...
            assert!(r1 == a);
            assert!(r2 == b);
            assert!(s1 == c);
            assert!(s2 == d);
            assert!(check == u32::MAX);

            // Benchmark the Tate pairing
            let bench_id = format!(
                "Compression of points of order 2^{} over field of {} bits",
                e,
                <$Fq>::BIT_LENGTH,
            );
            cc.bench_function(&bench_id, |bb| {
                bb.iter(|| {
                    black_box(E0).point_compression(
                        &black_box(xP),
                        &black_box(xQ),
                        &black_box(xPQ),
                        &black_box(xR),
                        &black_box(xS),
                        &black_box(xRS),
                        black_box(f),
                        black_box(e),
                        &black_box(dA),
                        black_box(dA_BITLEN),
                        black_box(&dlog_table),
                    )
                })
            });
        }

        criterion_group! {
            name = compression_benches;
            config = Criterion::default().measurement_time(Duration::from_secs(10));
            targets = bench_compression,
        }
    };
}

mod bench_compression_one {
    use cubical_pairing::kummer248::{Curve, Fq};
    use cubical_pairing::test_data::compression_one::*;

    use criterion::{black_box, criterion_group, criterion_main, Criterion};
    use std::time::Duration;

    define_compression_bench!(Curve, Fq);

    criterion_main!(compression_benches);
}

mod bench_compression_three {
    use cubical_pairing::kummer383::{Curve, Fq};
    use cubical_pairing::test_data::compression_three::*;

    use criterion::{black_box, criterion_group, criterion_main, Criterion};
    use std::time::Duration;

    define_compression_bench!(Curve, Fq);

    criterion_main!(compression_benches);
}

mod bench_compression_five {
    use cubical_pairing::kummer505::{Curve, Fq};
    use cubical_pairing::test_data::compression_five::*;

    use criterion::{black_box, criterion_group, criterion_main, Criterion};
    use std::time::Duration;

    define_compression_bench!(Curve, Fq);

    criterion_main!(compression_benches);
}

fn main() {
    // Run the benchmarks for all three sizes
    bench_compression_one::compression_benches();
    bench_compression_three::compression_benches();
    bench_compression_five::compression_benches();
}
