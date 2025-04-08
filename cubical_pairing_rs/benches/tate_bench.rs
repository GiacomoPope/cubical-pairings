#![allow(non_snake_case)]
#![allow(non_upper_case_globals)]
#![allow(unused_macros)] // macro used for benchmarking upsets the linter
#![allow(unused_imports)]

macro_rules! define_tate_pairing_bench {
    ($Curve:ty, $Fq:ty) => {
        fn bench_tate_pairing_2exp(c: &mut Criterion) {
            // Create E0
            // A supersingular elliptic curve of order 2^74 * 3^41
            let (A0, _) = <$Fq>::decode(&hex::decode(A0_str).unwrap());
            let E0 = <$Curve>::new(&A0);

            // Point of order 2^ea x(P)
            let (xP, _) = <$Fq>::decode(&hex::decode(P2X_str).unwrap());

            // Point of order 2^ea x(Q)
            let (xQ, _) = <$Fq>::decode(&hex::decode(Q2X_str).unwrap());

            // Difference point x(P - Q)
            let (xPQ, _) = <$Fq>::decode(&hex::decode(P2mQ2X_str).unwrap());

            // Make sure we're benchmarking something which works!!
            let (eQ2P2, _) = <$Fq>::decode(&hex::decode(tate_pairing_A_str).unwrap());
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

        fn bench_tate_pairing_odd(c: &mut Criterion) {
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

            // Make sure we're benchmarking something which works!!
            let (eQ3P3, _) = <$Fq>::decode(&hex::decode(tate_pairing_B_str).unwrap());
            let eQ3P3_test = E0.tate_pairing(&xP, &xQ, &xPQ, &B, B_BITLEN, &dB, dB_BITLEN);
            assert!(eQ3P3.equals(&eQ3P3_test) == u32::MAX);

            // Benchmark the Tate pairing
            let bench_id = format!(
                "Tate Pairing ({} bit field) order 3^{}",
                <$Fq>::BIT_LENGTH,
                eb
            );
            c.bench_function(&bench_id, |b| {
                b.iter(|| {
                    black_box(E0).tate_pairing(
                        &black_box(xP),
                        &black_box(xQ),
                        &black_box(xPQ),
                        &black_box(B),
                        black_box(B_BITLEN),
                        &black_box(dB),
                        black_box(dB_BITLEN),
                    )
                })
            });
        }

        fn bench_tate_pairing_even(c: &mut Criterion) {
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

            // Make sure we're benchmarking something which works!!
            let (eQP, _) = <$Fq>::decode(&hex::decode(tate_pairing_str).unwrap());
            let eQP_test = E0.tate_pairing(&xP, &xQ, &xPQ, &N, N_BITLEN, &dN, dN_BITLEN);
            assert!(eQP.equals(&eQP_test) == u32::MAX);

            // Benchmark the Tate pairing
            let bench_id = format!(
                "Tate Pairing ({} bit field) order 2^{} * 3^{}",
                <$Fq>::BIT_LENGTH,
                ea,
                eb
            );
            c.bench_function(&bench_id, |b| {
                b.iter(|| {
                    black_box(E0).tate_pairing(
                        &black_box(xP),
                        &black_box(xQ),
                        &black_box(xPQ),
                        &black_box(N),
                        black_box(N_BITLEN),
                        &black_box(dN),
                        black_box(dN_BITLEN),
                    )
                })
            });
        }

        criterion_group! {
            name = tate_benches;
            config = Criterion::default().measurement_time(Duration::from_secs(10));
            targets = bench_tate_pairing_2exp, bench_tate_pairing_odd, bench_tate_pairing_even
        }
    };
}

pub(crate) use define_tate_pairing_bench;
