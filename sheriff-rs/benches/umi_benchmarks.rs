use criterion::{black_box, criterion_group, criterion_main, Criterion};
use sheriff_rs::UmiDeduplicator;

fn benchmark_umi_deduplicator(c: &mut Criterion) {
    c.bench_function("umi_deduplicator_creation", |b| {
        b.iter(|| UmiDeduplicator::new(black_box(12)))
    });

    c.bench_function("umi_add_single", |b| {
        let mut dedup = UmiDeduplicator::new(12);
        b.iter(|| dedup.add(black_box("ACGTACGTACGT".to_string())))
    });

    c.bench_function("umi_contains", |b| {
        let mut dedup = UmiDeduplicator::new(12);
        dedup.add("ACGTACGTACGT".to_string());
        b.iter(|| dedup.contains(black_box("ACGTACGTACGT")))
    });
}

criterion_group!(benches, benchmark_umi_deduplicator);
criterion_main!(benches);
