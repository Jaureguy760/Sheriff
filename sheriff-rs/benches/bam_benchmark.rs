use criterion::{black_box, criterion_group, criterion_main, Criterion, BenchmarkId};
use sheriff_rs::kmer::{kmer_to_num, count_kmers};

fn bench_kmer_hashing(c: &mut Criterion) {
    let mut group = c.benchmark_group("kmer_hashing");

    // Benchmark different k-mer lengths
    for k in [3, 4, 5, 6].iter() {
        group.bench_with_input(BenchmarkId::from_parameter(k), k, |b, &k| {
            let sequence = b"AACGTACGTACGTACGT";
            b.iter(|| {
                for window in sequence.windows(k) {
                    black_box(kmer_to_num(window));
                }
            });
        });
    }

    group.finish();
}

fn bench_kmer_counting(c: &mut Criterion) {
    let mut group = c.benchmark_group("kmer_counting");

    // Benchmark counting with different sequence lengths
    let sequences = [
        ("short", b"AACGTACGTACGT".as_ref()),
        ("medium", b"AACGTACGTACGTAACGTACGTACGTAACGTACGTACGT".as_ref()),
        ("long", b"AACGTACGTACGTAACGTACGTACGTAACGTACGTACGTAACGTACGTACGTAACGTACGTACGT".as_ref()),
    ];

    for (name, seq) in sequences.iter() {
        group.bench_with_input(BenchmarkId::from_parameter(name), seq, |b, &seq| {
            b.iter(|| {
                black_box(count_kmers(seq, 4));
            });
        });
    }

    group.finish();
}

criterion_group!(benches, bench_kmer_hashing, bench_kmer_counting);
criterion_main!(benches);
