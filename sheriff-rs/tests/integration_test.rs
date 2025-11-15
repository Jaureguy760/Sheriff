use sheriff_rs::bam_filter::{filter_bam_by_barcodes, load_whitelist};
use std::collections::HashSet;

#[test]
fn test_load_whitelist_integration() {
    use std::io::Write;

    // Create temporary whitelist file
    let whitelist_path = "/tmp/test_whitelist_integration.txt";
    let mut file = std::fs::File::create(whitelist_path).unwrap();
    writeln!(file, "AAACCTGAGAAACCAT-1").unwrap();
    writeln!(file, "AAACCTGAGAAACCGC-1").unwrap();
    writeln!(file, "AAACCTGAGAAACCTA-1").unwrap();

    // Load whitelist
    let whitelist = load_whitelist(whitelist_path).unwrap();

    assert_eq!(whitelist.len(), 3);
    assert!(whitelist.contains("AAACCTGAGAAACCAT-1"));
    assert!(whitelist.contains("AAACCTGAGAAACCGC-1"));
    assert!(whitelist.contains("AAACCTGAGAAACCTA-1"));

    // Cleanup
    std::fs::remove_file(whitelist_path).unwrap();
}

#[test]
#[ignore] // Ignore by default, requires test BAM files
fn test_filter_bam_integration() {
    // This test requires actual BAM files
    // Run with: cargo test --ignored -- test_filter_bam_integration

    let input_bam = "tests/data/test_input.bam";
    let output_bam = "/tmp/test_output_integration.bam";
    let whitelist_path = "tests/data/whitelist.txt";

    // Check if test files exist
    if !std::path::Path::new(input_bam).exists() {
        eprintln!("Skipping test: {} not found", input_bam);
        eprintln!("Run: samtools view -h -s 0.01 ../example_data/barcode_headAligned_anno.sorted.edit_regions_200kb.bam | samtools view -b > {}", input_bam);
        return;
    }

    // Load whitelist
    let whitelist = load_whitelist(whitelist_path).unwrap();

    // Filter BAM
    let result = filter_bam_by_barcodes(input_bam, output_bam, &whitelist).unwrap();

    println!("Integration Test Results:");
    println!("  Reads processed: {}", result.reads_processed);
    println!("  Reads kept: {}", result.reads_kept);
    println!("  Reads rejected: {}", result.reads_rejected);

    assert!(result.reads_processed > 0, "Should process some reads");
    assert!(result.reads_kept < result.reads_processed, "Should filter some reads");

    // Cleanup
    std::fs::remove_file(output_bam).unwrap();
}
