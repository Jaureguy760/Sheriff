use criterion::{black_box, criterion_group, criterion_main, Criterion, BenchmarkId};
use sheriff_rs::{
    deduplicate_umis_unionfind,
    deduplicate_umis_unionfind_simd,
    deduplicate_umis_bktree,
    deduplicate_umis_adaptive,
    within_hamming_threshold,
};
use sheriff_rs::simd::{hamming_distance_simd, within_hamming_threshold_simd};

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

// SIMD Benchmarks

fn benchmark_hamming_distance_simd_vs_scalar(c: &mut Criterion) {
    let mut group = c.benchmark_group("hamming_simd_vs_scalar");

    for umi_len in [8, 12, 16, 32, 64].iter() {
        let umi1 = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"[..*umi_len].as_bytes();
        let umi2 = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCC"[..*umi_len].as_bytes();

        group.bench_with_input(
            BenchmarkId::new("scalar", umi_len),
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

        group.bench_with_input(
            BenchmarkId::new("simd", umi_len),
            umi_len,
            |b, _| {
                b.iter(|| {
                    within_hamming_threshold_simd(
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

fn benchmark_hamming_distance_full_simd_vs_scalar(c: &mut Criterion) {
    let mut group = c.benchmark_group("hamming_distance_full");

    for umi_len in [8, 12, 16, 32, 64].iter() {
        let umi1 = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"[..*umi_len].as_bytes();
        let umi2 = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCC"[..*umi_len].as_bytes();

        group.bench_with_input(
            BenchmarkId::new("scalar", umi_len),
            umi_len,
            |b, _| {
                b.iter(|| {
                    sheriff_rs::hamming_distance(
                        black_box(umi1),
                        black_box(umi2),
                    )
                })
            },
        );

        group.bench_with_input(
            BenchmarkId::new("simd", umi_len),
            umi_len,
            |b, _| {
                b.iter(|| {
                    hamming_distance_simd(
                        black_box(umi1),
                        black_box(umi2),
                    )
                })
            },
        );
    }

    group.finish();
}

fn benchmark_deduplicate_umis_simd_vs_scalar(c: &mut Criterion) {
    let mut group = c.benchmark_group("deduplicate_simd_vs_scalar");

    for n_umis in [10, 25, 50, 100].iter() {
        // Generate realistic UMIs (12bp)
        let umis: Vec<Vec<u8>> = (0..*n_umis)
            .map(|i| {
                format!(
                    "{}{}{}{}{}{}{}{}{}{}{}{}",
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
                    ["A", "C", "G", "T"][(i / 16) % 4],
                    ["A", "C", "G", "T"][(i / 64) % 4],
                )
                .into_bytes()
            })
            .collect();

        let umi_refs: Vec<&[u8]> = umis.iter().map(|u| u.as_slice()).collect();

        group.bench_with_input(
            BenchmarkId::new("scalar", n_umis),
            n_umis,
            |b, _| {
                b.iter(|| {
                    deduplicate_umis_unionfind(black_box(&umi_refs), black_box(1))
                })
            },
        );

        group.bench_with_input(
            BenchmarkId::new("simd", n_umis),
            n_umis,
            |b, _| {
                b.iter(|| {
                    deduplicate_umis_unionfind_simd(black_box(&umi_refs), black_box(1))
                })
            },
        );
    }

    group.finish();
}

fn benchmark_simd_early_exit(c: &mut Criterion) {
    let mut group = c.benchmark_group("simd_early_exit");

    // Test early exit with sequences that differ significantly
    let scenarios = vec![
        ("identical", b"ATCGATCGATCG".as_slice(), b"ATCGATCGATCG".as_slice()),
        ("1_mismatch", b"ATCGATCGATCG", b"TTCGATCGATCG"),
        ("2_mismatches", b"ATCGATCGATCG", b"TTCGATCGTTCG"),
        ("all_different", b"AAAAAAAAAAAA", b"TTTTTTTTTTTT"),
    ];

    for (name, umi1, umi2) in scenarios {
        group.bench_with_input(
            BenchmarkId::new("scalar", name),
            &(umi1, umi2),
            |b, (u1, u2)| {
                b.iter(|| {
                    within_hamming_threshold(
                        black_box(u1),
                        black_box(u2),
                        black_box(1),
                    )
                })
            },
        );

        group.bench_with_input(
            BenchmarkId::new("simd", name),
            &(umi1, umi2),
            |b, (u1, u2)| {
                b.iter(|| {
                    within_hamming_threshold_simd(
                        black_box(u1),
                        black_box(u2),
                        black_box(1),
                    )
                })
            },
        );
    }

    group.finish();
}

fn benchmark_simd_realistic_workload(c: &mut Criterion) {
    // Simulate a realistic Sheriff cell: 49 UMIs, 12bp each, threshold=1
    let umis: Vec<Vec<u8>> = (0..49)
        .map(|i| {
            format!(
                "{}{}{}{}{}{}{}{}{}{}{}{}",
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
                ["A", "C", "G", "T"][(i / 16) % 4],
                ["A", "C", "G", "T"][(i / 64) % 4],
            )
            .into_bytes()
        })
        .collect();

    let umi_refs: Vec<&[u8]> = umis.iter().map(|u| u.as_slice()).collect();

    let mut group = c.benchmark_group("realistic_workload_49_umis");

    group.bench_function("scalar", |b| {
        b.iter(|| deduplicate_umis_unionfind(black_box(&umi_refs), black_box(1)))
    });

    group.bench_function("simd", |b| {
        b.iter(|| deduplicate_umis_unionfind_simd(black_box(&umi_refs), black_box(1)))
    });

    group.finish();
}

// BK-Tree Benchmarks

fn benchmark_bktree_vs_bruteforce(c: &mut Criterion) {
    let mut group = c.benchmark_group("bktree_vs_bruteforce");

    // Test with varying UMI counts to find crossover point
    for n_umis in [10, 25, 50, 75, 100, 150, 200, 300, 500].iter() {
        // Generate diverse UMIs (realistic scenario)
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
            BenchmarkId::new("bruteforce", n_umis),
            n_umis,
            |b, _| {
                b.iter(|| {
                    deduplicate_umis_unionfind(black_box(&umi_refs), black_box(1))
                })
            },
        );

        group.bench_with_input(
            BenchmarkId::new("bktree", n_umis),
            n_umis,
            |b, _| {
                b.iter(|| {
                    deduplicate_umis_bktree(black_box(&umi_refs), black_box(1))
                })
            },
        );

        group.bench_with_input(
            BenchmarkId::new("adaptive", n_umis),
            n_umis,
            |b, _| {
                b.iter(|| {
                    deduplicate_umis_adaptive(black_box(&umi_refs), black_box(1))
                })
            },
        );
    }

    group.finish();
}

fn benchmark_bktree_worst_case(c: &mut Criterion) {
    // Worst case for BK-tree: all UMIs within threshold (dense clustering)
    // This tests when BK-tree degenerates to O(n²) like brute force
    let mut group = c.benchmark_group("bktree_worst_case");

    for n_umis in [50, 100, 200].iter() {
        let base_umi = b"AAAAAAAA";
        let umis: Vec<Vec<u8>> = (0..*n_umis)
            .map(|i| {
                let mut umi = base_umi.to_vec();
                // Change only 1 position to create threshold=1 neighbors
                umi[i % 8] = b"ACGT"[i % 4];
                umi
            })
            .collect();

        let umi_refs: Vec<&[u8]> = umis.iter().map(|u| u.as_slice()).collect();

        group.bench_with_input(
            BenchmarkId::new("bruteforce", n_umis),
            n_umis,
            |b, _| {
                b.iter(|| {
                    deduplicate_umis_unionfind(black_box(&umi_refs), black_box(1))
                })
            },
        );

        group.bench_with_input(
            BenchmarkId::new("bktree", n_umis),
            n_umis,
            |b, _| {
                b.iter(|| {
                    deduplicate_umis_bktree(black_box(&umi_refs), black_box(1))
                })
            },
        );
    }

    group.finish();
}

fn benchmark_bktree_best_case(c: &mut Criterion) {
    // Best case for BK-tree: all UMIs very different (sparse clustering)
    // This shows maximum benefit of BK-tree's O(n log n) average case
    let mut group = c.benchmark_group("bktree_best_case");

    for n_umis in [50, 100, 200, 500].iter() {
        // Generate UMIs that are all >1 mismatch apart
        let umis: Vec<Vec<u8>> = (0..*n_umis)
            .map(|i| {
                // Use modular arithmetic to ensure UMIs differ by >= 2 positions
                let offset = (i * 65537) % 256; // Prime number for good distribution
                format!(
                    "{}{}{}{}{}{}{}{}",
                    ["A", "C", "G", "T"][(offset / 1) % 4],
                    ["A", "C", "G", "T"][(offset / 4) % 4],
                    ["A", "C", "G", "T"][(offset / 16) % 4],
                    ["A", "C", "G", "T"][(offset / 64) % 4],
                    ["A", "C", "G", "T"][(offset / 256) % 4],
                    ["A", "C", "G", "T"][(offset / 1024) % 4],
                    ["A", "C", "G", "T"][(offset / 4096) % 4],
                    ["A", "C", "G", "T"][(offset / 16384) % 4],
                )
                .into_bytes()
            })
            .collect();

        let umi_refs: Vec<&[u8]> = umis.iter().map(|u| u.as_slice()).collect();

        group.bench_with_input(
            BenchmarkId::new("bruteforce", n_umis),
            n_umis,
            |b, _| {
                b.iter(|| {
                    deduplicate_umis_unionfind(black_box(&umi_refs), black_box(1))
                })
            },
        );

        group.bench_with_input(
            BenchmarkId::new("bktree", n_umis),
            n_umis,
            |b, _| {
                b.iter(|| {
                    deduplicate_umis_bktree(black_box(&umi_refs), black_box(1))
                })
            },
        );
    }

    group.finish();
}

fn benchmark_bktree_realistic_cells(c: &mut Criterion) {
    // Benchmark realistic Sheriff cell sizes
    let mut group = c.benchmark_group("bktree_realistic_cells");

    // Most Sheriff cells have 10-50 UMIs, some have 100-500
    for n_umis in [10, 20, 30, 40, 50, 100, 200].iter() {
        let umis: Vec<Vec<u8>> = (0..*n_umis)
            .map(|i| {
                format!(
                    "{}{}{}{}{}{}{}{}{}{}{}{}",
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
                    ["A", "C", "G", "T"][(i / 16) % 4],
                    ["A", "C", "G", "T"][(i / 64) % 4],
                )
                .into_bytes()
            })
            .collect();

        let umi_refs: Vec<&[u8]> = umis.iter().map(|u| u.as_slice()).collect();

        group.bench_with_input(
            BenchmarkId::new("bruteforce", n_umis),
            n_umis,
            |b, _| {
                b.iter(|| {
                    deduplicate_umis_unionfind(black_box(&umi_refs), black_box(1))
                })
            },
        );

        group.bench_with_input(
            BenchmarkId::new("bktree", n_umis),
            n_umis,
            |b, _| {
                b.iter(|| {
                    deduplicate_umis_bktree(black_box(&umi_refs), black_box(1))
                })
            },
        );

        group.bench_with_input(
            BenchmarkId::new("adaptive", n_umis),
            n_umis,
            |b, _| {
                b.iter(|| {
                    deduplicate_umis_adaptive(black_box(&umi_refs), black_box(1))
                })
            },
        );
    }

    group.finish();
}

criterion_group!(
    benches,
    benchmark_hamming_distance,
    benchmark_deduplicate_umis_small,
    benchmark_deduplicate_umis_scaling,
    benchmark_deduplicate_umis_real_cell,
    benchmark_deduplicate_umis_threshold,
    benchmark_deduplicate_umis_worst_case,
    benchmark_deduplicate_umis_best_case,
    // SIMD benchmarks
    benchmark_hamming_distance_simd_vs_scalar,
    benchmark_hamming_distance_full_simd_vs_scalar,
    benchmark_deduplicate_umis_simd_vs_scalar,
    benchmark_simd_early_exit,
    benchmark_simd_realistic_workload,
    // BK-tree benchmarks
    benchmark_bktree_vs_bruteforce,
    benchmark_bktree_worst_case,
    benchmark_bktree_best_case,
    benchmark_bktree_realistic_cells
);
criterion_main!(benches);
