#!/usr/bin/env python3
"""
Generate MD5 checksums for Sheriff test data validation.

Usage:
    python test_data/generate_checksums.py --generate  # Generate from current outputs
    python test_data/generate_checksums.py --verify    # Verify against expected
"""

import hashlib
import json
import argparse
from pathlib import Path

SCRIPT_DIR = Path(__file__).parent.absolute()
CHECKSUM_FILE = SCRIPT_DIR / "expected_checksums.json"


def md5_file(filepath):
    """Calculate MD5 checksum of a file."""
    hash_md5 = hashlib.md5()
    with open(filepath, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()


def md5_text_file(filepath):
    """Calculate MD5 of text file (normalized line endings)."""
    hash_md5 = hashlib.md5()
    with open(filepath, "r") as f:
        content = f.read()
        # Normalize line endings for cross-platform consistency
        content = content.replace('\r\n', '\n').rstrip('\n')
        hash_md5.update(content.encode('utf-8'))
    return hash_md5.hexdigest()


def generate_checksums(output_dir=None):
    """Generate checksums for test files and outputs."""
    checksums = {
        "_comment": "MD5 checksums for Sheriff test data validation",
        "_generated_by": "generate_checksums.py",
        "test_files": {},
        "expected_outputs": {}
    }

    # Checksum test input files
    test_files = [
        "test_200kb.bam",
        "test_200kb.bam.bai",
        "barcodes.txt",
        "edit_sites.txt",
        "blacklist.bed",
        "blacklist_seqs.txt"
    ]

    print("Generating checksums for test input files...")
    for fname in test_files:
        fpath = SCRIPT_DIR / fname
        if fpath.exists():
            if fname.endswith(('.bam', '.bai')):
                checksum = md5_file(fpath)
            else:
                checksum = md5_text_file(fpath)
            checksums["test_files"][fname] = checksum
            print(f"  {fname}: {checksum}")
        else:
            print(f"  {fname}: NOT FOUND")

    # Checksum output files (if output_dir provided)
    if output_dir:
        output_dir = Path(output_dir)
        print(f"\nGenerating checksums for outputs in {output_dir}...")

        for outfile in sorted(output_dir.glob("*")):
            if outfile.is_file():
                if outfile.suffix in ['.csv', '.txt', '.tsv']:
                    checksum = md5_text_file(outfile)
                else:
                    checksum = md5_file(outfile)
                checksums["expected_outputs"][outfile.name] = checksum
                print(f"  {outfile.name}: {checksum}")

    # Save checksums
    with open(CHECKSUM_FILE, "w") as f:
        json.dump(checksums, f, indent=2)

    print(f"\nChecksums saved to: {CHECKSUM_FILE}")
    return checksums


def verify_checksums(output_dir):
    """Verify outputs match expected checksums."""
    if not CHECKSUM_FILE.exists():
        print(f"ERROR: {CHECKSUM_FILE} not found")
        print("Run with --generate first to create expected checksums")
        return False

    with open(CHECKSUM_FILE) as f:
        expected = json.load(f)

    output_dir = Path(output_dir)
    all_pass = True

    print("Verifying test input files...")
    for fname, expected_md5 in expected.get("test_files", {}).items():
        fpath = SCRIPT_DIR / fname
        if not fpath.exists():
            print(f"  {fname}: MISSING")
            all_pass = False
            continue

        if fname.endswith(('.bam', '.bai')):
            actual_md5 = md5_file(fpath)
        else:
            actual_md5 = md5_text_file(fpath)

        if actual_md5 == expected_md5:
            print(f"  {fname}: OK")
        else:
            print(f"  {fname}: MISMATCH")
            print(f"    Expected: {expected_md5}")
            print(f"    Got:      {actual_md5}")
            all_pass = False

    print("\nVerifying output files...")
    for fname, expected_md5 in expected.get("expected_outputs", {}).items():
        if expected_md5 is None:
            print(f"  {fname}: SKIPPED (no expected checksum)")
            continue

        fpath = output_dir / fname
        if not fpath.exists():
            print(f"  {fname}: MISSING")
            all_pass = False
            continue

        if fname.endswith(('.csv', '.txt', '.tsv')):
            actual_md5 = md5_text_file(fpath)
        else:
            actual_md5 = md5_file(fpath)

        if actual_md5 == expected_md5:
            print(f"  {fname}: OK")
        else:
            print(f"  {fname}: MISMATCH")
            print(f"    Expected: {expected_md5}")
            print(f"    Got:      {actual_md5}")
            all_pass = False

    if all_pass:
        print("\n✅ All checksums verified!")
    else:
        print("\n❌ Checksum verification FAILED")

    return all_pass


def main():
    parser = argparse.ArgumentParser(description="Sheriff test data checksum utility")
    parser.add_argument("--generate", action="store_true", help="Generate checksums")
    parser.add_argument("--verify", action="store_true", help="Verify checksums")
    parser.add_argument("--output-dir", type=str, help="Directory with Sheriff outputs to verify")

    args = parser.parse_args()

    if args.generate:
        generate_checksums(args.output_dir)
    elif args.verify:
        if not args.output_dir:
            print("ERROR: --output-dir required for verification")
            return
        verify_checksums(args.output_dir)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
