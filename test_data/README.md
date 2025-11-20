Sheriff chr19 Test Data (Fork Release)
======================================

This repo does **not** track large test files. Download the small, single–chromosome test bundle from the Jaureguy760 fork release and unpack it locally:

```bash
mkdir -p test_data && cd test_data
wget -O- https://github.com/Jaureguy760/Sheriff/releases/download/test-data-v1/test_data_chr19.tar.gz | tar xvzf -

# Optional: verify checksum
echo "eb2f9c4c0d52930ad318c1218ee552fae4ab4c2144b63bcf9f388b3b419f07f2  test_data_chr19.tar.gz" | sha256sum -c -
```

Contents (after extract):
- `test_chr19.bam` / `.bai` (~42k reads from hg38_19:58–59 Mb)
- `chr19.fa.gz`, `chr19.gtf.gz` (GRCh38 chr19)
- `barcodes_chr19.txt` (11,006 barcodes)
- `blacklist_seqs.txt`, `blacklist.bed`, `edit_sites_chr19.txt`, `expected_checksums.json`
- `README.md`, `validate_chr19_test.py`, `ci_validation.py`

Quick validation:
```bash
PYTHONPATH=. python validate_chr19_test.py
```

Example Rust-accelerated run:
```bash
PYTHONPATH=. sheriff run test_chr19.bam chr19.fa barcodes_chr19.txt /iblm/netapp/data4/bbalderson/edit_capture/data/t7_indel_calling/Homo_sapiens.GRCh38.110.gtf \
  --blacklist blacklist.bed --blacklist_seqs blacklist_seqs.txt \
  --edit_site_min_cells 1 --chunk 50 --cpu 4 -o out
```

Release asset checksum:
- `test_data_chr19.tar.gz`: `eb2f9c4c0d52930ad318c1218ee552fae4ab4c2144b63bcf9f388b3b419f07f2`
