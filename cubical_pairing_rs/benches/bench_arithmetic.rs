mod util;

macro_rules! define_fq_benchmarks {
    ($Fq:ty, $ext_degree: literal) => {
        use criterion::{black_box, criterion_group, criterion_main, Criterion};
        use std::time::Duration;

        fn benchmark_add(c: &mut Criterion) {
            let mut rng = crate::util::DRNG::new();

            let x = <$Fq>::rand(&mut rng);
            let y = <$Fq>::rand(&mut rng);

            let bench_id = format!(
                "Benchmarking x + y with Fp^k with p ~2^{}, k = {}",
                <$Fq>::BIT_LENGTH,
                $ext_degree
            );
            c.bench_function(&bench_id, |b| b.iter(|| black_box(x) + black_box(y)));
        }

        fn benchmark_sub(c: &mut Criterion) {
            let mut rng = crate::util::DRNG::new();

            let x = <$Fq>::rand(&mut rng);
            let y = <$Fq>::rand(&mut rng);

            let bench_id = format!(
                "Benchmarking x - y with Fp^k with p ~2^{}, k = {}",
                <$Fq>::BIT_LENGTH,
                $ext_degree
            );
            c.bench_function(&bench_id, |b| b.iter(|| black_box(x) - black_box(y)));
        }

        fn benchmark_mul(c: &mut Criterion) {
            let mut rng = crate::util::DRNG::new();

            let x = <$Fq>::rand(&mut rng);
            let y = <$Fq>::rand(&mut rng);

            let bench_id = format!(
                "Benchmarking x * y with Fp^k with p ~2^{}, k = {}",
                <$Fq>::BIT_LENGTH,
                $ext_degree
            );
            c.bench_function(&bench_id, |b| b.iter(|| black_box(x) * black_box(y)));
        }

        fn benchmark_sqr(c: &mut Criterion) {
            let mut rng = crate::util::DRNG::new();

            let x = <$Fq>::rand(&mut rng);

            let bench_id = format!(
                "Benchmarking x^2 with Fp^k with p ~2^{}, k = {}",
                <$Fq>::BIT_LENGTH,
                $ext_degree
            );
            c.bench_function(&bench_id, |b| b.iter(|| black_box(x).square()));
        }

        fn benchmark_invert(c: &mut Criterion) {
            let mut rng = crate::util::DRNG::new();

            let x = <$Fq>::rand(&mut rng);

            let bench_id = format!(
                "Benchmarking 1 / x with Fp^k with p ~2^{}, k = {}",
                <$Fq>::BIT_LENGTH,
                $ext_degree
            );
            c.bench_function(&bench_id, |b| b.iter(|| black_box(x).invert()));
        }

        fn benchmark_sqrt(c: &mut Criterion) {
            let mut rng = crate::util::DRNG::new();

            let x = <$Fq>::rand(&mut rng);
            let y = x * x;

            let bench_id = format!(
                "Benchmarking sqrt(x) with Fp^k with p ~2^{}, k = {}",
                <$Fq>::BIT_LENGTH,
                $ext_degree
            );
            c.bench_function(&bench_id, |b| b.iter(|| black_box(y).sqrt()));
        }

        criterion_group! {
            name = fq_benchmarks;
            config = Criterion::default().measurement_time(Duration::from_secs(3));
            targets = benchmark_add, benchmark_sub, benchmark_mul, benchmark_sqr, benchmark_invert, benchmark_sqrt
        }
    };
}

// mod bench_csidh_512 {
//     define_fq_benchmarks!(cubical_pairing::csidh_fields::Fp511::Fp, 1);
//     criterion_main!(fq_benchmarks);
// }

// mod bench_dctidh_2048 {
//     define_fq_benchmarks!(cubical_pairing::csidh_fields::Fp2047::Fp, 1);
//     criterion_main!(fq_benchmarks);
// }

mod bench_sqi_lvl1 {
    define_fq_benchmarks!(cubical_pairing::sqisign_fields::Fp248Ext::Fp2, 2);
    criterion_main!(fq_benchmarks);
}

fn main() {
    bench_sqi_lvl1::fq_benchmarks();
    // bench_csidh_512::fq_benchmarks();
    // bench_dctidh_2048::fq_benchmarks();
}
