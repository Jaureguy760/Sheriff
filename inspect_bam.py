#!/usr/bin/env python3
"""Inspect BAM file to see what tags are available"""

import pysam

bam_file = "/home/user/Sheriff/example_data/barcode_headAligned_anno.sorted.edit_regions_200kb.bam"

print("Inspecting first 10 reads from BAM...")
print("=" * 70)

with pysam.AlignmentFile(bam_file, "rb") as bam:
    for i, read in enumerate(bam):
        if i >= 10:
            break

        print(f"\nRead {i+1}:")
        print(f"  Query name: {read.query_name}")
        print(f"  Sequence length: {len(read.query_sequence) if read.query_sequence else 0}")
        print(f"  Tags: {read.get_tags()}")
