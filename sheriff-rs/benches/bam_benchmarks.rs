use criterion::{black_box, criterion_group, criterion_main, Criterion, BenchmarkId};
use sheriff_rs::{BamReader, BamRecord};

fn benchmark_bam_record_creation(c: &mut Criterion) {
    c.bench_function("bam_record_creation", |b| {
        b.iter(|| {
            BamRecord::new(
                black_box("read1".to_string()),
                black_box("ACGTACGTACGT".to_string()),
                black_box("IIIIIIIIIIII".to_string()),
            )
        })
    });
}

fn benchmark_bam_record_seq_len(c: &mut Criterion) {
    let record = BamRecord::new(
        "read1".to_string(),
        "ACGTACGTACGT".to_string(),
        "IIIIIIIIIIII".to_string(),
    );

    c.bench_function("bam_record_seq_len", |b| {
        b.iter(|| black_box(&record).seq_len())
    });
}

fn benchmark_bam_reader_creation(c: &mut Criterion) {
    c.bench_function("bam_reader_creation", |b| {
        b.iter(|| BamReader::new(black_box("/path/to/file.bam".to_string())))
    });
}

fn benchmark_bam_record_clone(c: &mut Criterion) {
    let record = BamRecord::new(
        "read1".to_string(),
        "ACGTACGTACGT".to_string(),
        "IIIIIIIIIIII".to_string(),
    );

    c.bench_function("bam_record_clone", |b| {
        b.iter(|| black_box(&record).clone())
    });
}

fn benchmark_bam_record_different_sizes(c: &mut Criterion) {
    let mut group = c.benchmark_group("bam_record_sizes");

    for seq_len in [50, 100, 200, 500].iter() {
        let seq: String = (0..*seq_len).map(|i| {
            match i % 4 {
                0 => 'A',
                1 => 'C',
                2 => 'G',
                _ => 'T',
            }
        }).collect();
        let qual: String = "I".repeat(*seq_len);

        group.bench_with_input(
            BenchmarkId::from_parameter(seq_len),
            seq_len,
            |b, _| {
                b.iter(|| {
                    BamRecord::new(
                        black_box("read1".to_string()),
                        black_box(seq.clone()),
                        black_box(qual.clone()),
                    )
                })
            },
        );
    }

    group.finish();
}

criterion_group!(
    benches,
    benchmark_bam_record_creation,
    benchmark_bam_record_seq_len,
    benchmark_bam_reader_creation,
    benchmark_bam_record_clone,
    benchmark_bam_record_different_sizes
);
criterion_main!(benches);
