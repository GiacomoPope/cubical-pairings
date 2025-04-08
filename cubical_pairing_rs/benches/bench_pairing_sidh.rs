#![allow(non_snake_case)]
#![allow(non_upper_case_globals)]

mod tate_bench;
mod weil_bench;

mod bench_sidh_one {
    use cubical_pairing::kummer434::{Curve, Fq};
    use cubical_pairing::test_data::sidh_one::*;

    use criterion::{black_box, criterion_group, criterion_main, Criterion};
    use std::time::Duration;

    crate::tate_bench::define_tate_pairing_bench!(Curve, Fq);
    crate::weil_bench::define_weil_pairing_bench!(Curve, Fq);

    criterion_main!(weil_benches, tate_benches);
}

mod bench_sidh_three {
    use cubical_pairing::kummer610::{Curve, Fq};
    use cubical_pairing::test_data::sidh_three::*;

    use criterion::{black_box, criterion_group, criterion_main, Criterion};
    use std::time::Duration;

    crate::tate_bench::define_tate_pairing_bench!(Curve, Fq);
    crate::weil_bench::define_weil_pairing_bench!(Curve, Fq);

    criterion_main!(weil_benches, tate_benches);
}

mod bench_sidh_five {
    use cubical_pairing::kummer751::{Curve, Fq};
    use cubical_pairing::test_data::sidh_five::*;

    use criterion::{black_box, criterion_group, criterion_main, Criterion};
    use std::time::Duration;

    crate::tate_bench::define_tate_pairing_bench!(Curve, Fq);
    crate::weil_bench::define_weil_pairing_bench!(Curve, Fq);

    criterion_main!(weil_benches, tate_benches);
}

fn main() {
    // Run the benchmarks for all three sizes
    bench_sidh_one::weil_benches();
    bench_sidh_one::tate_benches();

    bench_sidh_three::weil_benches();
    bench_sidh_three::tate_benches();

    bench_sidh_five::weil_benches();
    bench_sidh_five::tate_benches();
}
