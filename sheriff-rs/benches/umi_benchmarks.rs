use criterion::{black_box, criterion_group, criterion_main, Criterion, BenchmarkId};
use sheriff_rs::{deduplicate_umis_unionfind, within_hamming_threshold};

fn benchmark_hamming_distance(c: &mut Criterion) {
    let mut group = c.benchmark_group("hamming_distance");

    for umi_len in [8, 10, 12, 16].iter() {
        let umi1 = "ACGTACGTACGTACGT"[..*umi_len].as_bytes();
        let umi2 = "ACGTACGTACGTACGG"[..*umi_len].as_bytes(); // 1 mismatch

        group.bench_with_input(
            BenchmarkId::from_parameter(umi_len),
            umi_len,
            |b, _| {
                b.iter(|| {
                    within_hamming_threshold(
                        black_box(umi1),
                        black_box(umi2),
                        black_box(1),
                    )
                })
            },
        );
    }

    group.finish();
}

fn benchmark_deduplicate_umis_small(c: &mut Criterion) {
    // Simulate real Sheriff cell (10-50 UMIs)
    let umis: Vec<&[u8]> = vec![
        b"ATCGATCG",
        b"ATCGATCC", // 1 mismatch from first
        b"GCGCGCGC",
        b"GCGCGCGT", // 1 mismatch from third
        b"TATATATA",
        b"TATATAAA", // 1 mismatch from fifth
        b"CGCGCGCG",
        b"AGAGAGAG",
        b"CTCTCTCT",
        b"GTGTGTGT",
    ];

    c.bench_function("deduplicate_umis_10", |b| {
        b.iter(|| deduplicate_umis_unionfind(black_box(&umis), black_box(1)))
    });
}

fn benchmark_deduplicate_umis_scaling(c: &mut Criterion) {
    let mut group = c.benchmark_group("deduplicate_umis_scaling");

    for n_umis in [10, 25, 50, 100, 200].iter() {
        // Generate UMIs with some duplicates
        let umis: Vec<Vec<u8>> = (0..*n_umis)
            .map(|i| {
                format!(
                    "{}{}{}{}{}{}{}{}",
                    ["A", "C", "G", "T"][i % 4],
                    ["A", "C", "G", "T"][(i / 4) % 4],
                    ["A", "C", "G", "T"][(i / 16) % 4],
                    ["A", "C", "G", "T"][(i / 64) % 4],
                    ["A", "C", "G", "T"][(i / 256) % 4],
                    ["A", "C", "G", "T"][(i / 1024) % 4],
                    ["A", "C", "G", "T"][(i / 4096) % 4],
                    ["A", "C", "G", "T"][(i / 16384) % 4],
                )
                .into_bytes()
            })
            .collect();

        let umi_refs: Vec<&[u8]> = umis.iter().map(|u| u.as_slice()).collect();

        group.bench_with_input(
            BenchmarkId::from_parameter(n_umis),
            n_umis,
            |b, _| {
                b.iter(|| {
                    deduplicate_umis_unionfind(black_box(&umi_refs), black_box(1))
                })
            },
        );
    }

    group.finish();
}

fn benchmark_deduplicate_umis_real_cell(c: &mut Criterion) {
    // Simulate real Sheriff cell with 49 UMIs (from benchmark_real_data.py)
    let umis: Vec<Vec<u8>> = (0..49)
        .map(|i| {
            format!(
                "{}{}{}{}{}{}{}{}{}{}",
                ["A", "C", "G", "T"][i % 4],
                ["A", "C", "G", "T"][(i / 4) % 4],
                ["A", "C", "G", "T"][(i / 16) % 4],
                ["A", "C", "G", "T"][(i / 64) % 4],
                ["A", "C", "G", "T"][(i / 256) % 4],
                ["A", "C", "G", "T"][(i / 1024) % 4],
                ["A", "C", "G", "T"][(i / 4096) % 4],
                ["A", "C", "G", "T"][(i / 16384) % 4],
                ["A", "C", "G", "T"][i % 4],
                ["A", "C", "G", "T"][(i / 4) % 4],
            )
            .into_bytes()
        })
        .collect();

    let umi_refs: Vec<&[u8]> = umis.iter().map(|u| u.as_slice()).collect();

    c.bench_function("deduplicate_umis_real_49", |b| {
        b.iter(|| deduplicate_umis_unionfind(black_box(&umi_refs), black_box(1)))
    });
}

fn benchmark_deduplicate_umis_threshold(c: &mut Criterion) {
    let mut group = c.benchmark_group("deduplicate_umis_threshold");

    let umis: Vec<Vec<u8>> = (0..50)
        .map(|i| {
            format!(
                "{}{}{}{}{}{}{}{}",
                ["A", "C", "G", "T"][i % 4],
                ["A", "C", "G", "T"][(i / 4) % 4],
                ["A", "C", "G", "T"][(i / 16) % 4],
                ["A", "C", "G", "T"][(i / 64) % 4],
                ["A", "C", "G", "T"][(i / 256) % 4],
                ["A", "C", "G", "T"][(i / 1024) % 4],
                ["A", "C", "G", "T"][(i / 4096) % 4],
                ["A", "C", "G", "T"][(i / 16384) % 4],
            )
            .into_bytes()
        })
        .collect();

    let umi_refs: Vec<&[u8]> = umis.iter().map(|u| u.as_slice()).collect();

    for threshold in [1, 2, 3].iter() {
        group.bench_with_input(
            BenchmarkId::from_parameter(threshold),
            threshold,
            |b, &t| {
                b.iter(|| {
                    deduplicate_umis_unionfind(black_box(&umi_refs), black_box(t))
                })
            },
        );
    }

    group.finish();
}

fn benchmark_deduplicate_umis_worst_case(c: &mut Criterion) {
    // Worst case: all UMIs are within threshold=1 (forms one giant cluster)
    let base_umi = b"AAAAAAAA";
    let umis: Vec<Vec<u8>> = (0..50)
        .map(|i| {
            let mut umi = base_umi.to_vec();
            // Change only 1 position to create threshold=1 neighbors
            umi[i % 8] = b"ACGT"[i % 4];
            umi
        })
        .collect();

    let umi_refs: Vec<&[u8]> = umis.iter().map(|u| u.as_slice()).collect();

    c.bench_function("deduplicate_umis_worst_case_50", |b| {
        b.iter(|| deduplicate_umis_unionfind(black_box(&umi_refs), black_box(1)))
    });
}

fn benchmark_deduplicate_umis_best_case(c: &mut Criterion) {
    // Best case: all UMIs are completely different (no clustering needed)
    let umis: Vec<&[u8]> = vec![
        b"AAAAAAAA",
        b"CCCCCCCC",
        b"GGGGGGGG",
        b"TTTTTTTT",
        b"ACACACAC",
        b"GTGTGTGT",
        b"AGAGAGAG",
        b"CTCTCTCT",
        b"ATATATAT",
        b"GCGCGCGC",
    ];

    c.bench_function("deduplicate_umis_best_case_10", |b| {
        b.iter(|| deduplicate_umis_unionfind(black_box(&umis), black_box(1)))
    });
}

criterion_group!(
    benches,
    benchmark_hamming_distance,
    benchmark_deduplicate_umis_small,
    benchmark_deduplicate_umis_scaling,
    benchmark_deduplicate_umis_real_cell,
    benchmark_deduplicate_umis_threshold,
    benchmark_deduplicate_umis_worst_case,
    benchmark_deduplicate_umis_best_case
);
criterion_main!(benches);
