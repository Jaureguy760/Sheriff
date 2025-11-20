use std::collections::HashMap;
use ahash::AHashSet;
use rust_htslib::bam::{self, Read};

/// Compute gene UMI counts per cell using BAM tags (CB for barcode, GX for gene id, pN for UMI).
///
/// Returns a matrix[G][C] of u32 counts in gene_ids order and barcodes order.
pub fn gene_counts_per_cell(
    bam_path: &str,
    barcodes: &[String],
    gene_ids: &[String],
) -> Result<Vec<Vec<u32>>, Box<dyn std::error::Error>> {
    let mut barcode_to_idx = HashMap::with_capacity(barcodes.len());
    for (i, bc) in barcodes.iter().enumerate() {
        barcode_to_idx.insert(bc.as_str(), i);
    }

    let mut gene_to_idx = HashMap::with_capacity(gene_ids.len());
    for (i, g) in gene_ids.iter().enumerate() {
        gene_to_idx.insert(g.as_str(), i);
    }

    let n_genes = gene_ids.len();
    let n_cells = barcodes.len();

    // umi_sets[gene][cell] = unique UMIs
    let mut umi_sets: Vec<Vec<AHashSet<String>>> = vec![vec![AHashSet::new(); n_cells]; n_genes];

    let mut reader = bam::Reader::from_path(bam_path)?;

    for result in reader.records() {
        let record = result?;

        // Cell barcode
        let cb = match record.aux(b"CB") {
            Ok(bam::record::Aux::String(s)) => s.to_string(),
            _ => continue,
        };

        let cell_idx = match barcode_to_idx.get(cb.as_str()) {
            Some(i) => *i,
            None => continue,
        };

        // Gene id
        let gx = match record.aux(b"GX") {
            Ok(bam::record::Aux::String(s)) => s.to_string(),
            _ => continue,
        };

        let gene_idx = match gene_to_idx.get(gx.as_str()) {
            Some(i) => *i,
            None => continue,
        };

        // UMI tag
        let umi = match record.aux(b"pN") {
            Ok(bam::record::Aux::String(s)) => s.to_string(),
            _ => continue,
        };

        umi_sets[gene_idx][cell_idx].insert(umi);
    }

    // Convert to counts
    let counts: Vec<Vec<u32>> = umi_sets
        .into_iter()
        .map(|cells| cells.into_iter().map(|umis| umis.len() as u32).collect())
        .collect();

    Ok(counts)
}
