# Simpler and faster pairings from the Montgomery Ladder

Code accompanying the research paper "Simpler and faster pairings from the Montgomery Ladder" by Giacomo Pope, Krijn Reijnders, Damien Robert, Alessandro Sferlazza and Benjamin Smith.

This repository includes three distinct implementations:

- `cubical_pairing_rs`: a constant time and efficient implementation of cubical pairings for elliptic curves $E / \mathbb{F}_{p^2}$ for applications in isogeny-based cryptography
- `cubical_pairings_sage`: a flexible and highly commented SageMath implementation of cubical pairings for elliptic curves $E / \mathbb{F}_{q}$
- `operation_count_sage`: a detailed implementation with described optimisations in the paper which explicitly counts the number of $\mathbb{F}_{p}$ operations. Not intended for general use.

### Running the Rust code

The code was written and tested with `rustc 1.87.0-nightly`. To compute benchmarks of the Rust code, we recommend running

```
RUSTFLAGS="-C target-cpu=native" cargo bench
```

To run the tests, simply run

```
cargo test
```
