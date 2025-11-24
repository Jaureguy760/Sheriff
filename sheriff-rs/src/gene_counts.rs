use std::collections::HashMap;

use ahash::AHashSet;
use rust_htslib::bam::{self, record::Aux, Read};
use smallvec::SmallVec;

/// Compute unique gene UMI counts per cell from a barcoded BAM file (Rust).
///
/// Reads these tags per record:
/// - **CB**: cell barcode
/// - **GX/GN/gn**: gene identifier (multiple tag names are supported)
/// - **pN/UB**: UMI sequence (pipeline-dependent UMI tag)
///
/// For each read that has all required tags and whose barcode/gene appear in
/// the provided lists, the UMI is inserted into a per-(gene, cell) set so each
/// UMI is counted once. Returns a dense matrix `[genes][cells]` of `u32`
/// counts in the order of `gene_ids` (rows) and `barcodes` (cols).
///
/// Implementation notes:
/// - Uses byte-keyed hash maps for barcodes/genes to avoid per-record `String`
///   allocations on tag parsing.
/// - Uses `SmallVec<[u8; 24]>` inside an `AHashSet` for UMI storage to reduce
///   heap churn for short UMIs.
/// - Honors `max_reads` as an optional early-exit cap (0 = no limit).
///
/// Errors if the BAM cannot be opened/read; records missing expected tags are
/// skipped without failing the run.
pub fn gene_counts_per_cell(
    bam_path: &str,
    barcodes: &[String],
    gene_ids: &[String],
    max_reads: usize,
) -> Result<Vec<Vec<u32>>, Box<dyn std::error::Error>> {
    // Map barcodes/genes to indices using byte keys to avoid per-record String allocs.
    let mut barcode_to_idx: HashMap<Vec<u8>, usize> = HashMap::with_capacity(barcodes.len());
    for (i, bc) in barcodes.iter().enumerate() {
        barcode_to_idx.insert(bc.as_bytes().to_vec(), i);
    }

    let mut gene_to_idx: HashMap<Vec<u8>, usize> = HashMap::with_capacity(gene_ids.len());
    for (i, g) in gene_ids.iter().enumerate() {
        gene_to_idx.insert(g.as_bytes().to_vec(), i);
    }

    let n_genes = gene_ids.len();
    let n_cells = barcodes.len();

    // umi_sets[gene][cell] = unique UMI byte buffers
    type UmiBuf = SmallVec<[u8; 24]>;
    let mut umi_sets: Vec<Vec<AHashSet<UmiBuf>>> = vec![vec![AHashSet::new(); n_cells]; n_genes];

    let mut reader = bam::Reader::from_path(bam_path)?;

    for (read_idx, result) in reader.records().enumerate() {
        if max_reads > 0 && read_idx >= max_reads {
            break;
        }

        let record = result?;

        // Cell barcode as bytes
        let cb = match record.aux(b"CB") {
            Ok(Aux::String(s)) => s.as_bytes(),
            _ => continue,
        };

        let cell_idx = match barcode_to_idx.get(cb) {
            Some(i) => *i,
            None => continue,
        };

        // Gene id - try multiple tag names (GX, GN, gn)
        let gx_bytes = if let Ok(Aux::String(s)) = record.aux(b"GX") {
            s.as_bytes()
        } else if let Ok(Aux::String(s)) = record.aux(b"GN") {
            s.as_bytes()
        } else if let Ok(Aux::String(s)) = record.aux(b"gn") {
            s.as_bytes()
        } else {
            continue;
        };

        let gene_idx = match gene_to_idx.get(gx_bytes) {
            Some(i) => *i,
            None => continue,
        };

        // UMI tag - try multiple tag names (pN, UB)
        let umi_bytes = if let Ok(Aux::String(s)) = record.aux(b"pN") {
            s.as_bytes()
        } else if let Ok(Aux::String(s)) = record.aux(b"UB") {
            s.as_bytes()
        } else {
            continue;
        };

        // Copy UMI bytes into a smallvec to avoid heap churn for short UMIs.
        let mut buf: UmiBuf = SmallVec::new();
        buf.extend_from_slice(umi_bytes);
        umi_sets[gene_idx][cell_idx].insert(buf);
    }

    let counts: Vec<Vec<u32>> = umi_sets
        .into_iter()
        .map(|cells| cells.into_iter().map(|umis| umis.len() as u32).collect())
        .collect();

    Ok(counts)
}
