use criterion::{black_box, criterion_group, criterion_main, Criterion};
use sheriff_rs::KmerCounter;

fn benchmark_kmer_counter(c: &mut Criterion) {
    c.bench_function("kmer_counter_creation", |b| {
        b.iter(|| KmerCounter::new(black_box(21)))
    });

    c.bench_function("kmer_counter_len", |b| {
        let counter = KmerCounter::new(21);
        b.iter(|| counter.len())
    });
}

criterion_group!(benches, benchmark_kmer_counter);
criterion_main!(benches);
