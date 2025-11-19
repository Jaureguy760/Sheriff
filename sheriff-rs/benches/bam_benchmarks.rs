use criterion::{black_box, criterion_group, criterion_main, Criterion, BenchmarkId};
use sheriff_rs::bam::{BamProcessor, collect_stats, collect_stats_parallel, get_umi_and_barcode, process_records, process_records_parallel, group_by_cell_parallel};

const BAM_PATH: &str = "../example_data/barcode_headAligned_anno.sorted.edit_regions_200kb.bam";

fn benchmark_bam_open(c: &mut Criterion) {
    if !std::path::Path::new(BAM_PATH).exists() {
        eprintln!("Skipping benchmark: BAM file not found at {}", BAM_PATH);
        return;
    }

    c.bench_function("bam_processor_open", |b| {
        b.iter(|| {
            let processor = BamProcessor::new(black_box(BAM_PATH));
            assert!(processor.is_ok());
        })
    });
}

fn benchmark_bam_read_all_records(c: &mut Criterion) {
    if !std::path::Path::new(BAM_PATH).exists() {
        eprintln!("Skipping benchmark: BAM file not found");
        return;
    }

    c.bench_function("bam_read_all_records", |b| {
        b.iter(|| {
            let mut processor = BamProcessor::new(BAM_PATH).unwrap();
            let mut count = 0;

            process_records(&mut processor, |_record| {
                count += 1;
                Ok(())
            })
            .unwrap();

            black_box(count)
        })
    });
}

fn benchmark_bam_extract_tags(c: &mut Criterion) {
    if !std::path::Path::new(BAM_PATH).exists() {
        eprintln!("Skipping benchmark: BAM file not found");
        return;
    }

    c.bench_function("bam_extract_tags_all", |b| {
        b.iter(|| {
            let mut processor = BamProcessor::new(BAM_PATH).unwrap();
            let mut found_tags = 0;

            process_records(&mut processor, |record| {
                if let Some((_umi, _cb)) = get_umi_and_barcode(record) {
                    found_tags += 1;
                }
                Ok(())
            })
            .unwrap();

            black_box(found_tags)
        })
    });
}

fn benchmark_bam_collect_stats(c: &mut Criterion) {
    if !std::path::Path::new(BAM_PATH).exists() {
        eprintln!("Skipping benchmark: BAM file not found");
        return;
    }

    c.bench_function("bam_collect_stats_all", |b| {
        b.iter(|| {
            let stats = collect_stats(black_box(BAM_PATH)).unwrap();
            black_box(stats)
        })
    });
}

fn benchmark_bam_per_read_cost(c: &mut Criterion) {
    if !std::path::Path::new(BAM_PATH).exists() {
        eprintln!("Skipping benchmark: BAM file not found");
        return;
    }

    // First, count total reads
    let mut processor = BamProcessor::new(BAM_PATH).unwrap();
    let mut total_reads = 0;
    process_records(&mut processor, |_record| {
        total_reads += 1;
        Ok(())
    })
    .unwrap();

    c.bench_function("bam_per_read_amortized", |b| {
        b.iter(|| {
            let mut processor = BamProcessor::new(BAM_PATH).unwrap();
            let mut count = 0;

            process_records(&mut processor, |record| {
                if let Some((_umi, _cb)) = get_umi_and_barcode(record) {
                    count += 1;
                }
                Ok(())
            })
            .unwrap();

            black_box(count)
        });

        // Print per-read cost
        println!("\nTotal reads: {}", total_reads);
    });
}

fn benchmark_bam_scaling(c: &mut Criterion) {
    if !std::path::Path::new(BAM_PATH).exists() {
        eprintln!("Skipping benchmark: BAM file not found");
        return;
    }

    let mut group = c.benchmark_group("bam_scaling");

    // Benchmark processing different numbers of reads
    for &limit in [1000, 10000, 50000, 100000, 352535].iter() {
        group.bench_with_input(BenchmarkId::from_parameter(limit), &limit, |b, &limit| {
            b.iter(|| {
                let mut processor = BamProcessor::new(BAM_PATH).unwrap();
                let mut count = 0;

                process_records(&mut processor, |record| {
                    if count >= limit {
                        // Early exit (not ideal but works for benchmark)
                        return Ok(());
                    }

                    if let Some((_umi, _cb)) = get_umi_and_barcode(record) {
                        count += 1;
                    }
                    Ok(())
                })
                .ok(); // Ignore errors from early exit

                black_box(count)
            })
        });
    }

    group.finish();
}

fn benchmark_bam_tag_extraction_only(c: &mut Criterion) {
    if !std::path::Path::new(BAM_PATH).exists() {
        eprintln!("Skipping benchmark: BAM file not found");
        return;
    }

    // Pre-load a record to benchmark tag extraction in isolation
    c.bench_function("bam_tag_extraction_single", |b| {
        let mut processor = BamProcessor::new(BAM_PATH).unwrap();

        // Get first record
        let mut first_record = None;
        process_records(&mut processor, |record| {
            first_record = Some(record.clone());
            Err(sheriff_rs::bam::BamError::ParseError("stop".to_string()))
        })
        .ok();

        let record = first_record.unwrap();

        b.iter(|| {
            let tags = get_umi_and_barcode(black_box(&record));
            black_box(tags)
        })
    });
}

fn benchmark_bam_collect_stats_parallel(c: &mut Criterion) {
    if !std::path::Path::new(BAM_PATH).exists() {
        eprintln!("Skipping benchmark: BAM file not found");
        return;
    }

    let mut group = c.benchmark_group("bam_stats_sequential_vs_parallel");

    group.bench_function("sequential", |b| {
        b.iter(|| {
            let stats = collect_stats(black_box(BAM_PATH)).unwrap();
            black_box(stats)
        })
    });

    group.bench_function("parallel", |b| {
        b.iter(|| {
            let stats = collect_stats_parallel(black_box(BAM_PATH)).unwrap();
            black_box(stats)
        })
    });

    group.finish();
}

fn benchmark_bam_parallel_speedup(c: &mut Criterion) {
    if !std::path::Path::new(BAM_PATH).exists() {
        eprintln!("Skipping benchmark: BAM file not found");
        return;
    }

    c.bench_function("bam_parallel_stats", |b| {
        b.iter(|| {
            let stats = collect_stats_parallel(black_box(BAM_PATH)).unwrap();
            black_box(stats)
        })
    });
}

fn benchmark_bam_group_by_cell(c: &mut Criterion) {
    if !std::path::Path::new(BAM_PATH).exists() {
        eprintln!("Skipping benchmark: BAM file not found");
        return;
    }

    c.bench_function("bam_group_by_cell_parallel", |b| {
        b.iter(|| {
            let cell_map = group_by_cell_parallel(black_box(BAM_PATH)).unwrap();
            black_box(cell_map)
        })
    });
}

fn benchmark_bam_parallel_processing(c: &mut Criterion) {
    if !std::path::Path::new(BAM_PATH).exists() {
        eprintln!("Skipping benchmark: BAM file not found");
        return;
    }

    c.bench_function("bam_process_records_parallel", |b| {
        b.iter(|| {
            use std::sync::atomic::{AtomicUsize, Ordering};
            let count = AtomicUsize::new(0);

            process_records_parallel(BAM_PATH, |record| {
                if let Some((_umi, _cb)) = get_umi_and_barcode(record) {
                    count.fetch_add(1, Ordering::Relaxed);
                }
            })
            .unwrap();

            black_box(count.load(Ordering::Relaxed))
        })
    });
}

criterion_group!(
    benches,
    benchmark_bam_open,
    benchmark_bam_read_all_records,
    benchmark_bam_extract_tags,
    benchmark_bam_collect_stats,
    benchmark_bam_per_read_cost,
    benchmark_bam_scaling,
    benchmark_bam_tag_extraction_only,
    benchmark_bam_collect_stats_parallel,
    benchmark_bam_parallel_speedup,
    benchmark_bam_group_by_cell,
    benchmark_bam_parallel_processing
);
criterion_main!(benches);
