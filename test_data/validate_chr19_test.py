#!/usr/bin/env python3
"""
Validate that the single-chromosome test dataset works with Sheriff.

This validates:
1. BAM can be read with correct barcodes
2. Reference FASTA can be accessed (coordinate mapping)
3. GTF annotations can be loaded
4. Sheriff's chromosome name mapping works (hg38_19 -> 19)
"""

import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).parent.absolute()
SHERIFF_ROOT = SCRIPT_DIR.parent
sys.path.insert(0, str(SHERIFF_ROOT))

import pysam


def test_bam_and_barcodes():
    """Test BAM reading and barcode matching."""
    print("=== Testing BAM and Barcodes ===")

    bam_file = SCRIPT_DIR / "test_chr19.bam"
    barcodes_file = SCRIPT_DIR / "barcodes_chr19.txt"

    # Load barcodes
    with open(barcodes_file) as f:
        barcodes = set(line.strip() for line in f if line.strip())
    print(f"  Barcodes loaded: {len(barcodes)}")

    # Read BAM
    with pysam.AlignmentFile(str(bam_file), "rb") as bam:
        total_reads = 0
        matched_reads = 0
        chromosomes = set()

        for read in bam.fetch():
            total_reads += 1
            chromosomes.add(read.reference_name)
            if read.has_tag("CB"):
                cb = read.get_tag("CB")
                if cb in barcodes:
                    matched_reads += 1

            if total_reads >= 1000:  # Quick check
                break

    print(f"  Reads checked: {total_reads}")
    print(f"  Matched barcodes: {matched_reads}")
    print(f"  Chromosomes: {chromosomes}")

    # Should be hg38_19
    assert "hg38_19" in chromosomes, f"Expected hg38_19, got {chromosomes}"
    assert matched_reads > 0, "No reads matched barcodes"

    print("  [OK] BAM and barcodes work correctly")
    return True


def test_reference_access():
    """Test that reference FASTA can be accessed."""
    print("\n=== Testing Reference FASTA ===")

    # Check for both compressed and uncompressed
    fa_gz = SCRIPT_DIR / "chr19.fa.gz"
    fa = SCRIPT_DIR / "chr19.fa"

    if fa.exists():
        ref_path = fa
    elif fa_gz.exists():
        print("  Decompressing chr19.fa.gz...")
        import gzip
        with gzip.open(fa_gz, 'rb') as f_in:
            with open(fa, 'wb') as f_out:
                f_out.write(f_in.read())
        ref_path = fa
        # Index it
        pysam.faidx(str(ref_path))
    else:
        print("  [FAIL] No reference FASTA found")
        return False

    # Test access
    fasta = pysam.FastaFile(str(ref_path))
    print(f"  Reference file: {ref_path.name}")
    print(f"  Chromosomes: {fasta.references}")

    # Should be just "19"
    assert "19" in fasta.references, f"Expected '19', got {fasta.references}"

    # Get sequence at known position (from edit sites)
    # Edit site: hg38_19:10986350 -> reference uses "19"
    seq = fasta.fetch("19", 10986340, 10986360)  # 20bp around edit site
    print(f"  Sample sequence at 19:10986340-10986360: {seq}")
    assert len(seq) == 20, f"Expected 20bp, got {len(seq)}"

    fasta.close()
    print("  [OK] Reference FASTA accessible")
    return True


def test_chromosome_mapping():
    """Test Sheriff's chromosome name mapping."""
    print("\n=== Testing Chromosome Name Mapping ===")

    # Import Sheriff's mapping function
    from sheriff.count_t7 import reformat_chr_name

    # Create mock read
    class MockRead:
        reference_name = "hg38_19"

    mock_read = MockRead()
    mapped_name = reformat_chr_name(mock_read)

    print(f"  BAM chrom: {mock_read.reference_name}")
    print(f"  Mapped name: {mapped_name}")

    assert mapped_name == "19", f"Expected '19', got {mapped_name}"
    print("  [OK] Chromosome mapping works (hg38_19 -> 19)")
    return True


def test_gtf_loading():
    """Test GTF can be loaded."""
    print("\n=== Testing GTF Annotations ===")

    gtf_gz = SCRIPT_DIR / "chr19.gtf.gz"
    gtf = SCRIPT_DIR / "chr19.gtf"

    if gtf.exists():
        gtf_path = gtf
    elif gtf_gz.exists():
        print("  Decompressing chr19.gtf.gz...")
        import gzip
        with gzip.open(gtf_gz, 'rb') as f_in:
            with open(gtf, 'wb') as f_out:
                f_out.write(f_in.read())
        gtf_path = gtf
    else:
        print("  [FAIL] No GTF found")
        return False

    # Count features
    gene_count = 0
    with open(gtf_path) as f:
        for line in f:
            if not line.startswith("#"):
                if "\tgene\t" in line:
                    gene_count += 1
            if gene_count >= 10:  # Quick check
                break

    print(f"  GTF file: {gtf_path.name}")
    print(f"  Genes found (first 10): {gene_count}")

    assert gene_count > 0, "No genes found in GTF"
    print("  [OK] GTF annotations accessible")
    return True


def test_edit_sites():
    """Test edit sites file."""
    print("\n=== Testing Edit Sites ===")

    edit_sites_file = SCRIPT_DIR / "edit_sites_chr19.txt"

    if not edit_sites_file.exists():
        print("  [FAIL] edit_sites_chr19.txt not found")
        return False

    with open(edit_sites_file) as f:
        sites = [line.strip().split() for line in f if line.strip()]

    print(f"  Edit sites: {len(sites)}")
    for chrom, start, end in sites:
        print(f"    {chrom}:{start}-{end}")

    # All should be hg38_19
    for chrom, start, end in sites:
        assert chrom == "hg38_19", f"Expected hg38_19, got {chrom}"

    print("  [OK] Edit sites loaded")
    return True


def main():
    print("Sheriff Single-Chromosome Test Validation")
    print("=" * 50)
    print(f"Test data directory: {SCRIPT_DIR}")
    print()

    results = {}
    results["bam_barcodes"] = test_bam_and_barcodes()
    results["reference"] = test_reference_access()
    results["chrom_mapping"] = test_chromosome_mapping()
    results["gtf"] = test_gtf_loading()
    results["edit_sites"] = test_edit_sites()

    # Summary
    print("\n" + "=" * 50)
    print("Summary")
    print("=" * 50)

    all_pass = all(results.values())
    for test_name, passed in results.items():
        status = "PASS" if passed else "FAIL"
        print(f"  {test_name}: {status}")

    if all_pass:
        print("\n✅ All validations passed!")
        print("\nThis test dataset is ready for full pipeline testing.")
        print("Package size: ~23MB (compressed)")
        print("\nFiles:")
        print("  - test_chr19.bam (3.7MB) - 42k reads")
        print("  - chr19.fa.gz (16MB) - GRCh38 chr19")
        print("  - chr19.gtf.gz (2.9MB) - Gene annotations")
        print("  - barcodes_chr19.txt - 4,680 cell barcodes")
        print("  - edit_sites_chr19.txt - 4 edit sites")
        sys.exit(0)
    else:
        print("\n❌ Some validations failed!")
        sys.exit(1)


if __name__ == "__main__":
    main()
