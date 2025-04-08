#![allow(non_snake_case)]
#![allow(non_upper_case_globals)]
#![allow(unused_macros)] // macro used for benchmarking upsets the linter
#![allow(unused_imports)]

macro_rules! define_weil_pairing_even_only_bench {
    ($Curve:ty, $Fq:ty) => {
        fn bench_weil_pairing_2exp(c: &mut Criterion) {
            // Create E0
            // A supersingular elliptic curve of order (p + 1)
            let (A0, _) = <$Fq>::decode(&hex::decode(A0_str).unwrap());
            let E0 = <$Curve>::new(&A0);

            // Point of order 2^ea x(P)
            let (xP, _) = <$Fq>::decode(&hex::decode(PX_str).unwrap());

            // Point of order 2^ea x(Q)
            let (xQ, _) = <$Fq>::decode(&hex::decode(QX_str).unwrap());

            // Sum x(P + Q)
            let (xPQ, _) = <$Fq>::decode(&hex::decode(PmQX_str).unwrap());

            // Make sure we're benchmarking something which works!!
            let (eQP, _) = <$Fq>::decode(&hex::decode(weil_pairing_str).unwrap());
            let eQP_test = E0.weil_pairing_2exp(&xP, &xQ, &xPQ, ea);
            assert!(eQP.equals(&eQP_test) == u32::MAX);

            // Benchmark the Weil pairing
            let bench_id = format!(
                "Weil Pairing ({} bit field) order 2^{}",
                <$Fq>::BIT_LENGTH,
                ea
            );
            c.bench_function(&bench_id, |b| {
                b.iter(|| {
                    black_box(E0).weil_pairing_2exp(
                        &black_box(xP),
                        &black_box(xQ),
                        &black_box(xPQ),
                        black_box(ea),
                    )
                })
            });
        }

        criterion_group! {
            name = weil_even_only_benches;
            config = Criterion::default().measurement_time(Duration::from_secs(10));
            targets = bench_weil_pairing_2exp
        }
    };
}

pub(crate) use define_weil_pairing_even_only_bench;
