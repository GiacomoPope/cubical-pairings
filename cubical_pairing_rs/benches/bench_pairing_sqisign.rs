#![allow(non_snake_case)]
#![allow(non_upper_case_globals)]

mod tate_bench_even_only;
mod weil_bench_even_only;

mod bench_sqisign_one {
    use cubical_pairing::kummer248::{Curve, Fq};
    use cubical_pairing::test_data::sqisign_one::*;

    use criterion::{black_box, criterion_group, criterion_main, Criterion};
    use std::time::Duration;

    crate::tate_bench_even_only::define_tate_pairing_even_only_bench!(Curve, Fq);
    crate::weil_bench_even_only::define_weil_pairing_even_only_bench!(Curve, Fq);

    criterion_main!(weil_even_only_benches, tate_even_only_benches);
}

mod bench_sqisign_three {
    use cubical_pairing::kummer383::{Curve, Fq};
    use cubical_pairing::test_data::sqisign_three::*;

    use criterion::{black_box, criterion_group, criterion_main, Criterion};
    use std::time::Duration;

    crate::tate_bench_even_only::define_tate_pairing_even_only_bench!(Curve, Fq);
    crate::weil_bench_even_only::define_weil_pairing_even_only_bench!(Curve, Fq);

    criterion_main!(weil_even_only_benches, tate_even_only_benches);
}

mod bench_sqisign_five {
    use cubical_pairing::kummer505::{Curve, Fq};
    use cubical_pairing::test_data::sqisign_five::*;

    use criterion::{black_box, criterion_group, criterion_main, Criterion};
    use std::time::Duration;

    crate::tate_bench_even_only::define_tate_pairing_even_only_bench!(Curve, Fq);
    crate::weil_bench_even_only::define_weil_pairing_even_only_bench!(Curve, Fq);

    criterion_main!(weil_even_only_benches, tate_even_only_benches);
}

fn main() {
    // Run the benchmarks for the first set of tests
    bench_sqisign_one::weil_even_only_benches();
    bench_sqisign_one::tate_even_only_benches();

    // Run the benchmarks for the first set of tests
    bench_sqisign_three::weil_even_only_benches();
    bench_sqisign_three::tate_even_only_benches();

    // Run the benchmarks for the first set of tests
    bench_sqisign_five::weil_even_only_benches();
    bench_sqisign_five::tate_even_only_benches();
}
