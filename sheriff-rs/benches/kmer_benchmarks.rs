use criterion::{black_box, criterion_group, criterion_main, Criterion, BenchmarkId};
use rustc_hash::FxHashSet;
use sheriff_rs::{kmer_to_num, match_kmer, match_kmer_rolling};

fn benchmark_kmer_to_num(c: &mut Criterion) {
    let mut group = c.benchmark_group("kmer_to_num");

    for k in [6, 8, 10, 12].iter() {
        let kmer = "ACGTACGTACGT"[..*k].as_bytes();
        group.bench_with_input(BenchmarkId::from_parameter(k), k, |b, _| {
            b.iter(|| kmer_to_num(black_box(kmer)))
        });
    }

    group.finish();
}

fn benchmark_match_kmer(c: &mut Criterion) {
    let mut group = c.benchmark_group("match_kmer");

    // Test different sequence lengths
    for seq_len in [100, 500, 1000, 5000].iter() {
        // Generate sequence
        let sequence: Vec<u8> = (0..*seq_len)
            .map(|i| match i % 4 {
                0 => b'A',
                1 => b'C',
                2 => b'G',
                _ => b'T',
            })
            .collect();

        // Create whitelist (5 k-mers from t7 barcode)
        let t7_barcode = b"GGGAGAGTAT";
        let k = 6;
        let mut whitelist = FxHashSet::default();
        for i in 0..=(t7_barcode.len() - k) {
            let hash = kmer_to_num(&t7_barcode[i..i + k]);
            whitelist.insert(hash);
        }

        group.bench_with_input(
            BenchmarkId::from_parameter(seq_len),
            seq_len,
            |b, _| {
                b.iter(|| {
                    match_kmer(
                        black_box(&sequence),
                        black_box(k),
                        black_box(&whitelist),
                    )
                })
            },
        );
    }

    group.finish();
}

fn benchmark_match_kmer_different_k(c: &mut Criterion) {
    let mut group = c.benchmark_group("match_kmer_k_size");

    let sequence: Vec<u8> = (0..1000)
        .map(|i| match i % 4 {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            _ => b'T',
        })
        .collect();

    for k in [4, 6, 8, 10, 12].iter() {
        // Create whitelist with k-mers
        let mut whitelist = FxHashSet::default();
        whitelist.insert(kmer_to_num(&b"AAAAAAAAAAAAA"[..*k]));
        whitelist.insert(kmer_to_num(&b"CCCCCCCCCCCCC"[..*k]));
        whitelist.insert(kmer_to_num(&b"GGGGGGGGGGGGG"[..*k]));
        whitelist.insert(kmer_to_num(&b"TTTTTTTTTTTTT"[..*k]));

        group.bench_with_input(BenchmarkId::from_parameter(k), k, |b, &k_size| {
            b.iter(|| {
                match_kmer(
                    black_box(&sequence),
                    black_box(k_size),
                    black_box(&whitelist),
                )
            })
        });
    }

    group.finish();
}

fn benchmark_match_kmer_real_data(c: &mut Criterion) {
    // Simulate real Sheriff read (198bp)
    let real_sequence = b"ATGCTGCTCAATGGAAAACAACCATAGGGAACTTCAGATCACATGTTAAGACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";

    let t7_barcode = b"GGGAGAGTAT";
    let k = 6;
    let mut whitelist = FxHashSet::default();
    for i in 0..=(t7_barcode.len() - k) {
        let hash = kmer_to_num(&t7_barcode[i..i + k]);
        whitelist.insert(hash);
    }

    c.bench_function("match_kmer_real_198bp", |b| {
        b.iter(|| {
            match_kmer(
                black_box(real_sequence),
                black_box(k),
                black_box(&whitelist),
            )
        })
    });
}

fn benchmark_match_kmer_rolling_vs_regular(c: &mut Criterion) {
    let mut group = c.benchmark_group("match_kmer_rolling_vs_regular");

    // Test on real Sheriff read length (198bp)
    let real_sequence = b"ATGCTGCTCAATGGAAAACAACCATAGGGAACTTCAGATCACATGTTAAGACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";

    // Build whitelist from t7 barcode
    let t7_barcode = b"GGGAGAGTAT";
    let k = 6;
    let mut whitelist = FxHashSet::default();
    for i in 0..=(t7_barcode.len() - k) {
        let hash = kmer_to_num(&t7_barcode[i..i + k]);
        whitelist.insert(hash);
    }

    group.bench_function("regular", |b| {
        b.iter(|| {
            match_kmer(
                black_box(real_sequence),
                black_box(k),
                black_box(&whitelist),
            )
        })
    });

    group.bench_function("rolling_hash", |b| {
        b.iter(|| {
            match_kmer_rolling(
                black_box(real_sequence),
                black_box(k),
                black_box(&whitelist),
            )
        })
    });

    group.finish();
}

fn benchmark_rolling_hash_different_k(c: &mut Criterion) {
    let mut group = c.benchmark_group("rolling_hash_k_size");

    let sequence: Vec<u8> = (0..1000)
        .map(|i| match i % 4 {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            _ => b'T',
        })
        .collect();

    for k in [4, 6, 8, 10, 12].iter() {
        // Create whitelist with k-mers
        let mut whitelist = FxHashSet::default();
        whitelist.insert(kmer_to_num(&b"AAAAAAAAAAAAA"[..*k]));
        whitelist.insert(kmer_to_num(&b"CCCCCCCCCCCCC"[..*k]));
        whitelist.insert(kmer_to_num(&b"GGGGGGGGGGGGG"[..*k]));
        whitelist.insert(kmer_to_num(&b"TTTTTTTTTTTTT"[..*k]));

        group.bench_with_input(BenchmarkId::from_parameter(k), k, |b, &k_size| {
            b.iter(|| {
                match_kmer_rolling(
                    black_box(&sequence),
                    black_box(k_size),
                    black_box(&whitelist),
                )
            })
        });
    }

    group.finish();
}

fn benchmark_rolling_hash_sequence_length(c: &mut Criterion) {
    let mut group = c.benchmark_group("rolling_hash_seq_length");

    let k = 6;
    let mut whitelist = FxHashSet::default();
    let t7_barcode = b"GGGAGAGTAT";
    for i in 0..=(t7_barcode.len() - k) {
        let hash = kmer_to_num(&t7_barcode[i..i + k]);
        whitelist.insert(hash);
    }

    for seq_len in [100, 200, 500, 1000, 5000].iter() {
        // Generate sequence
        let sequence: Vec<u8> = (0..*seq_len)
            .map(|i| match i % 4 {
                0 => b'A',
                1 => b'C',
                2 => b'G',
                _ => b'T',
            })
            .collect();

        group.bench_with_input(
            BenchmarkId::from_parameter(seq_len),
            seq_len,
            |b, _| {
                b.iter(|| {
                    match_kmer_rolling(
                        black_box(&sequence),
                        black_box(k),
                        black_box(&whitelist),
                    )
                })
            },
        );
    }

    group.finish();
}

criterion_group!(
    benches,
    benchmark_kmer_to_num,
    benchmark_match_kmer,
    benchmark_match_kmer_different_k,
    benchmark_match_kmer_real_data,
    benchmark_match_kmer_rolling_vs_regular,
    benchmark_rolling_hash_different_k,
    benchmark_rolling_hash_sequence_length
);
criterion_main!(benches);
