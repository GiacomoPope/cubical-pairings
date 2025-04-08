#![allow(non_snake_case)]
#![allow(non_upper_case_globals)]
#![allow(unused_macros)] // macro used for benchmarking upsets the linter
#![allow(unused_imports)]

macro_rules! define_tate_pairing_even_only_bench {
    ($Curve:ty, $Fq:ty) => {
        fn bench_tate_pairing_2exp(c: &mut Criterion) {
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

            // Make sure we're benchmarking something which works!!
            let (eQ2P2, _) = <$Fq>::decode(&hex::decode(tate_pairing_str).unwrap());
            let eQ2P2_test = E0.tate_pairing_2exp(&xP, &xQ, &xPQ, ea, &dA, dA_BITLEN);
            assert!(eQ2P2.equals(&eQ2P2_test) == u32::MAX);

            // Benchmark the Tate pairing
            let bench_id = format!(
                "Tate Pairing ({} bit field) order 2^{}",
                <$Fq>::BIT_LENGTH,
                ea
            );
            c.bench_function(&bench_id, |b| {
                b.iter(|| {
                    black_box(E0).tate_pairing_2exp(
                        &black_box(xP),
                        &black_box(xQ),
                        &black_box(xPQ),
                        black_box(ea),
                        &black_box(dA),
                        black_box(dA_BITLEN),
                    )
                })
            });
        }

        criterion_group! {
            name = tate_even_only_benches;
            config = Criterion::default().measurement_time(Duration::from_secs(10));
            targets = bench_tate_pairing_2exp
        }
    };
}

pub(crate) use define_tate_pairing_even_only_bench;
