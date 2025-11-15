"""
High-performance BAM utilities with optional Rust acceleration

This module provides BAM filtering operations with automatic detection
of Rust acceleration. If the Rust module is available, it will be used
for 10-50x performance improvement. Otherwise, falls back to pure Python.
"""
import pysam
import tempfile
import os
from typing import Set, Dict
import warnings

# Try to import Rust acceleration
try:
    import sheriff_rs
    HAS_RUST = True
except ImportError:
    HAS_RUST = False
    warnings.warn(
        "Rust acceleration not available. Using pure Python implementation. "
        "For 10-50x speedup, install Rust acceleration:\n"
        "  cd sheriff-rs && maturin build --release --features python && "
        "pip install target/wheels/sheriff_rs-*.whl",
        UserWarning
    )


def filter_bam_by_barcodes(
    input_bam: str,
    output_bam: str,
    cell_barcodes: Set[str],
    use_rust: bool = True,
    use_parallel: bool = False,
    verbose: bool = True
) -> Dict[str, int]:
    """
    Filter BAM file by cell barcode whitelist.

    Uses Rust implementation if available (10-50x faster),
    falls back to Python pysam if not.

    Args:
        input_bam: Input BAM file path
        output_bam: Output BAM file path
        cell_barcodes: Set of allowed cell barcodes
        use_rust: Prefer Rust if available (default True)
        use_parallel: Use parallel processing (rayon par_bridge, default False)
                      Provides 2-5x speedup on medium/large files
        verbose: Print performance information (default True)

    Returns:
        Dict with keys: reads_processed, reads_kept, reads_rejected
    """
    if HAS_RUST and use_rust:
        if verbose:
            mode = "parallel" if use_parallel else "sequential"
            print(f"  Using Rust acceleration for BAM filtering ({mode})...")

        # Use direct barcode list (avoids temporary file overhead)
        barcodes_list = list(cell_barcodes)

        # Choose sequential or parallel version
        if use_parallel:
            result = sheriff_rs.filter_bam_by_barcodes_rust_parallel(
                input_bam, output_bam, barcodes_list
            )
        else:
            result = sheriff_rs.filter_bam_by_barcodes_rust(
                input_bam, output_bam, barcodes_list
            )

        if verbose:
            duration = result.get('duration_seconds', 0)
            if duration > 0:
                throughput = result['reads_processed'] / duration
            else:
                throughput = 0
            mode_str = "Rust-parallel" if use_parallel else "Rust"
            print(f"    Processed {result['reads_processed']:,} reads")
            print(f"    Kept {result['reads_kept']:,} reads ({100*result['reads_kept']/max(1, result['reads_processed']):.1f}%)")
            print(f"    Duration: {duration:.2f}s")
            print(f"    Throughput: {throughput:,.0f} reads/sec ({mode_str})")

        return result
    else:
        if verbose:
            if not HAS_RUST:
                print(f"  Using Python implementation (Rust not available)...")
            else:
                print(f"  Using Python implementation (use_rust=False)...")

        return _filter_bam_python(input_bam, output_bam, cell_barcodes, verbose=verbose)


def _filter_bam_python(
    input_bam: str,
    output_bam: str,
    cell_barcodes: Set[str],
    verbose: bool = True
) -> Dict[str, int]:
    """Python fallback implementation using pysam"""
    import time

    start_time = time.time()

    bam = pysam.AlignmentFile(input_bam, 'rb')
    out = pysam.AlignmentFile(output_bam, 'wb', template=bam)

    stats = {
        'reads_processed': 0,
        'reads_kept': 0,
        'reads_rejected': 0
    }

    for read in bam:
        stats['reads_processed'] += 1

        try:
            cb = read.get_tag('CB')
            if cb in cell_barcodes:
                out.write(read)
                stats['reads_kept'] += 1
            else:
                stats['reads_rejected'] += 1
        except KeyError:
            stats['reads_rejected'] += 1

    duration = time.time() - start_time

    bam.close()
    out.close()

    if verbose:
        throughput = stats['reads_processed'] / max(0.001, duration)
        print(f"    Processed {stats['reads_processed']:,} reads")
        print(f"    Kept {stats['reads_kept']:,} reads ({100*stats['reads_kept']/max(1, stats['reads_processed']):.1f}%)")
        print(f"    Throughput: {throughput:,.0f} reads/sec (Python)")

    stats['duration_seconds'] = duration

    return stats
