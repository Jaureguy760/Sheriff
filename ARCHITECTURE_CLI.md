# Sheriff CLI Architecture Design

## Overview
Comprehensive redesign of Sheriff CLI to world-class standards with subcommands, config files, resume functionality, and plugin system.

## Design Principles

1. **Backward Compatibility**: Existing `sheriff get-t7-edits` continues to work
2. **Progressive Enhancement**: Features are optional, defaults match current behavior
3. **Professional UX**: Clear feedback, validation, error messages
4. **Extensibility**: Plugin system for custom filters/processors
5. **Production Ready**: Logging, profiling, checkpointing for large datasets

---

## Priority 2: YAML Config, Subcommands, Logging

### 1. Subcommand Architecture

**Current Structure:**
```
sheriff get-t7-edits [OPTIONS] bam ref barcodes gtf
```

**New Structure:**
```
sheriff [GLOBAL_OPTIONS] COMMAND [COMMAND_OPTIONS]

Global Options:
  --version, -V          Show version
  --config FILE          Load config from YAML
  --log-file FILE        Log to file
  --log-level LEVEL      Set log level (DEBUG, INFO, WARNING, ERROR)
  --profile              Enable profiling mode

Commands:
  run                    Run full Sheriff pipeline (default, replaces get-t7-edits)
  filter                 BAM filtering only (Rust accelerated)
  validate               Validate inputs without processing
  benchmark              Run performance benchmarks
  config                 Config file utilities
  interactive            Interactive mode with prompts
```

**Implementation Strategy:**
- Keep `get-t7-edits` as alias to `run` for backward compatibility
- Each subcommand is a separate function with its own options
- Global options handled via `app.callback()`
- Shared logic extracted to common functions

### 2. YAML Configuration

**File Format:** `sheriff-config.yaml`
```yaml
# Sheriff Configuration File
version: "1.2.0"

# Input files (required)
input:
  bam_file: "path/to/sample.bam"
  ref_file: "path/to/genome.fa"
  barcode_file: "path/to/barcodes.txt"
  gtf_file: "path/to/genes.gtf"

# Optional filter files
filters:
  blacklist_file: "path/to/blacklist.bed"
  whitelist_file: "path/to/whitelist.bed"
  blacklist_seqs: "path/to/seqs.txt"
  cnv_file: "path/to/cnv.tsv"

# Pipeline parameters
parameters:
  t7_barcode: "GGGAGAGTAT"
  k: 6
  edit_dist: 140
  stranded_edit_dist: 15
  edit_site_min_cells: 3
  nonbc_edit_dist: 1000
  ploidy: 2
  mrna_count_mode: "all"

# Processing options
processing:
  n_cpus: 8
  chunk_size_mb: 250
  verbosity: 1

# Output options
output:
  outdir: "results/"
  edit_site_rev_comp_filt: true
  uncorrected_gene_count: false

# Resume functionality
resume:
  enabled: false
  checkpoint_dir: ".sheriff_checkpoints/"
  checkpoint_interval: 5  # minutes

# Logging
logging:
  file: "sheriff.log"
  level: "INFO"
  format: "detailed"  # simple, detailed, json

# Profiling (optional)
profiling:
  enabled: false
  output: "sheriff_profile.json"
  memory_tracking: true
```

**Implementation:**
- Use PyYAML for parsing
- Pydantic for validation (add to dependencies)
- CLI args override config file values
- `sheriff config generate` creates template
- `sheriff config validate` checks config file

**Priority Order:**
1. CLI arguments (highest priority)
2. Config file values
3. Default values (lowest priority)

### 3. Structured Logging

**Architecture:**
```
Logger Hierarchy:
  sheriff (root)
    ├── sheriff.cli (CLI operations)
    ├── sheriff.pipeline (main pipeline)
    ├── sheriff.filtering (BAM filtering)
    ├── sheriff.kmer (k-mer matching)
    └── sheriff.io (file I/O)
```

**Features:**
- Dual output: console (INFO+) + file (DEBUG+)
- Structured format with timestamps, module, level
- Configurable via CLI or config file
- Log rotation for large files
- JSON format option for parsing

**Example Output:**
```
2025-11-15 10:23:45 [INFO] sheriff.cli: Starting Sheriff v1.2.0
2025-11-15 10:23:45 [INFO] sheriff.pipeline: Validating inputs...
2025-11-15 10:23:46 [INFO] sheriff.pipeline: ✓ All inputs valid
2025-11-15 10:23:46 [INFO] sheriff.filtering: Filtering BAM with Rust acceleration
2025-11-15 10:23:52 [INFO] sheriff.filtering: Processed 1.2M reads in 6.2s
```

---

## Priority 3: Resume, Interactive Mode, Rich UI

### 4. Checkpoint/Resume System

**Architecture:**
```
Checkpoint Format: JSON
Location: {outdir}/.sheriff_checkpoints/run_{timestamp}.json

{
  "version": "1.2.0",
  "timestamp": "2025-11-15T10:23:45",
  "status": "in_progress",
  "config": {...},  // Full configuration
  "progress": {
    "current_step": "kmer_matching",
    "steps_completed": ["validation", "bam_filtering"],
    "steps_remaining": ["kmer_matching", "edit_calling", "umi_counting"],
    "percent_complete": 40
  },
  "outputs": {
    "filtered_bam": "results/filtered.bam",
    "temp_files": [...]
  },
  "metrics": {
    "reads_processed": 1234567,
    "runtime_seconds": 125.3
  }
}
```

**Implementation:**
- Automatic checkpointing every 5 minutes (configurable)
- Checkpoint after each major step
- `--resume` flag auto-detects last checkpoint
- `--resume PATH` specifies checkpoint file
- Validates checkpoint compatibility (version, inputs)
- Cleans up old checkpoints on success

**Pipeline Steps** (for checkpointing):
1. Input validation
2. BAM filtering
3. K-mer matching
4. Edit site calling
5. UMI counting
6. Output generation

**Usage:**
```bash
# Normal run
sheriff run config.yaml

# If interrupted, resume:
sheriff run config.yaml --resume
# or
sheriff run --resume results/.sheriff_checkpoints/run_20251115_102345.json
```

### 5. Interactive Mode

**Flow:**
```
$ sheriff interactive

🔬 Sheriff Interactive Setup

Let's configure your CRISPR edit site analysis.

📁 Input Files
  BAM file: [browse] > /path/to/sample.bam
  ✓ File exists (2.3 GB, 12.5M reads)

  Reference genome: [browse] > /path/to/hg38.fa
  ✓ File exists (3.1 GB)

  Barcodes: [browse] > /path/to/barcodes.txt
  ✓ File exists (5,247 barcodes)

  GTF file: [browse] > /path/to/genes.gtf
  ✓ File exists (58,381 genes)

⚙️  Parameters
  T7 barcode [GGGAGAGTAT]: <enter for default>
  K-mer size [6]: <enter for default>
  Edit distance [140]: 200
  Min cells per edit [3]: 5

🖥️  Processing
  CPUs [1]: 8
  ✓ 8 CPUs available

📤 Output
  Output directory [current]: results/

💾 Save configuration?
  [y] Save to config file: sheriff-config.yaml

🚀 Ready to run!
  [r] Run now
  [q] Quit and save config

> r

Running Sheriff pipeline...
```

**Implementation:**
- Rich prompts with validation
- File browser integration (questionary library)
- Auto-detection of available resources
- Generates YAML config for reproducibility
- Can exit and save config without running

### 6. Rich Tables and Progress

**Features:**

**A. Progress Bars:**
```
Filtering BAM...
[████████████████░░░░] 75% | 9.4M/12.5M reads | 45s elapsed | ~15s remaining
```

**B. Summary Tables:**
```
╭─────────────────────── Sheriff Results ────────────────────────╮
│                                                                 │
│  📊 Pipeline Summary                                            │
│  ────────────────────────────────────────────────────────────  │
│  Input Reads              12,456,789                            │
│  Filtered Reads            8,234,123  (66%)                     │
│  T7 Barcoded Reads           45,678  (0.37%)                    │
│  Edit Sites Detected             23                             │
│  Cells with Edits             1,234                             │
│  Genes Quantified            18,452                             │
│                                                                 │
│  ⏱️  Performance                                                │
│  ────────────────────────────────────────────────────────────  │
│  Total Runtime              1h 23m                              │
│  BAM Filtering                12m (Rust)                        │
│  K-mer Matching               45m (Rust)                        │
│  UMI Counting                 26m                               │
│                                                                 │
│  📁 Output Files                                                │
│  ────────────────────────────────────────────────────────────  │
│  Edit Sites          results/edit_site_info.txt                 │
│  UMI Counts          results/umi_counts.mtx                     │
│  Gene Counts         results/gene_counts.mtx                    │
│  Log File            results/sheriff.log                        │
│                                                                 │
╰─────────────────────────────────────────────────────────────────╯
```

**C. Status Panels:**
```
╭─ Current Operation ──────────────────────────────╮
│ K-mer Matching                                   │
│ Chromosome: chr7                                 │
│ Matches found: 1,234                             │
│ Speed: 25,000 reads/s                            │
╰──────────────────────────────────────────────────╯
```

**Implementation:**
- Rich Progress for bars
- Rich Table for summary
- Rich Panel for status
- Rich Console for colored output
- Live updating display

---

## Priority 4: Profiling and Plugin System

### 7. Profiling Mode

**Enabled via:**
```bash
sheriff run config.yaml --profile
# or in config:
profiling:
  enabled: true
  output: "sheriff_profile.json"
```

**Metrics Collected:**
```json
{
  "version": "1.2.0",
  "timestamp": "2025-11-15T10:23:45",
  "system": {
    "os": "Linux 4.4.0",
    "cpu_count": 8,
    "memory_total_gb": 32.0,
    "python_version": "3.10.12"
  },
  "input": {
    "bam_size_gb": 2.3,
    "total_reads": 12456789
  },
  "performance": {
    "total_runtime_seconds": 4980.2,
    "steps": {
      "validation": {
        "runtime_seconds": 2.1,
        "memory_peak_mb": 125
      },
      "bam_filtering": {
        "runtime_seconds": 720.5,
        "reads_per_second": 17283,
        "memory_peak_mb": 856,
        "accelerator": "rust_parallel"
      },
      "kmer_matching": {
        "runtime_seconds": 2700.3,
        "matches_per_second": 461,
        "memory_peak_mb": 1240,
        "accelerator": "rust_iterative"
      },
      "umi_counting": {
        "runtime_seconds": 1560.0,
        "genes_per_second": 12,
        "memory_peak_mb": 2048
      }
    }
  },
  "speedup": {
    "vs_python_estimated": "4.5x",
    "rust_components": ["bam_filtering", "kmer_matching"]
  }
}
```

**Implementation:**
- Memory profiling with `memory_profiler`
- Time profiling with custom decorators
- Resource monitoring with `psutil`
- Export to JSON for analysis
- Optional: HTML report generation

### 8. Plugin System

**Architecture:**
```
Plugin Directory: ~/.sheriff/plugins/ or {project}/.sheriff/plugins/

Plugin Structure:
my_plugin/
  ├── __init__.py
  ├── plugin.yaml      # Plugin metadata
  └── filters.py       # Plugin implementation

plugin.yaml:
name: "my_custom_filter"
version: "1.0.0"
sheriff_version: ">=1.2.0"
type: "bam_filter"
entry_point: "filters:CustomFilter"
```

**Plugin Types:**
1. **BAM Filters**: Custom read filtering logic
2. **Edit Validators**: Custom edit site validation
3. **Output Formatters**: Custom output formats
4. **Reporters**: Custom analysis reports

**Plugin API:**
```python
# Example: Custom BAM filter plugin

from sheriff.plugin_api import BAMFilterPlugin

class CustomFilter(BAMFilterPlugin):
    """Filter reads by custom criteria."""

    name = "custom_quality_filter"
    version = "1.0.0"

    def filter_read(self, read) -> bool:
        """Return True to keep read, False to filter out."""
        # Custom logic
        return read.mapping_quality >= 30 and read.is_proper_pair

    def get_config(self):
        """Return plugin configuration schema."""
        return {
            "min_quality": {"type": int, "default": 30},
            "require_proper_pair": {"type": bool, "default": True}
        }
```

**Usage:**
```bash
# List plugins
sheriff plugins list

# Install plugin
sheriff plugins install my_plugin/

# Use in config
plugins:
  - name: "custom_quality_filter"
    enabled: true
    config:
      min_quality: 35
```

**Implementation:**
- Plugin discovery via directory scanning
- Dynamic loading with `importlib`
- Validation against Sheriff version
- Sandboxing for safety
- Plugin registry and management commands

---

## Implementation Plan

### Phase 1: Priority 2 (Core Infrastructure)
**Time: 2-3 hours**

1. **YAML Config Support** (45 min)
   - Add pydantic, PyYAML to dependencies
   - Create config schema with Pydantic
   - Implement config loading and merging
   - Add `sheriff config` subcommand

2. **Subcommand Restructure** (60 min)
   - Refactor to subcommand architecture
   - Implement: run, filter, validate, benchmark, config
   - Maintain backward compatibility with get-t7-edits

3. **Structured Logging** (45 min)
   - Set up logger hierarchy
   - Dual output (console + file)
   - Configurable levels and formats
   - Integrate into existing code

### Phase 2: Priority 3 (UX Enhancement)
**Time: 3-4 hours**

4. **Checkpoint/Resume** (90 min)
   - Design checkpoint format
   - Implement checkpoint saving
   - Implement resume detection and loading
   - Integrate into pipeline steps

5. **Interactive Mode** (60 min)
   - Add questionary for prompts
   - Implement interactive flow
   - File browser and validation
   - Config generation from interactive session

6. **Rich UI** (90 min)
   - Progress bars for long operations
   - Summary tables at completion
   - Status panels for live updates
   - Integrate into count_t7.py

### Phase 3: Priority 4 (Advanced Features)
**Time: 2-3 hours**

7. **Profiling Mode** (60 min)
   - Add psutil, memory_profiler
   - Implement metric collection
   - Decorators for timing
   - JSON export

8. **Plugin System** (90 min)
   - Define plugin API
   - Implement discovery and loading
   - Create example plugin
   - Plugin management commands

### Phase 4: Testing & Documentation
**Time: 1-2 hours**

9. **Testing** (60 min)
   - Test all new features
   - Backward compatibility tests
   - Example configs and plugins

10. **Documentation** (30 min)
    - Update README
    - CLI usage guide
    - Config file reference
    - Plugin development guide

---

## File Structure (New)

```
sheriff/
├── __init__.py
├── __main__.py              # CLI entry point (refactored)
├── cli/                     # NEW: CLI subcommands
│   ├── __init__.py
│   ├── run.py              # Main pipeline command
│   ├── filter.py           # BAM filtering command
│   ├── validate.py         # Validation command
│   ├── benchmark.py        # Benchmarking command
│   ├── config_cmd.py       # Config utilities
│   └── interactive.py      # Interactive mode
├── config/                  # NEW: Configuration
│   ├── __init__.py
│   ├── schema.py           # Pydantic models
│   └── loader.py           # Config loading/merging
├── checkpoint/              # NEW: Resume functionality
│   ├── __init__.py
│   ├── manager.py          # Checkpoint management
│   └── format.py           # Checkpoint format
├── plugin_api/              # NEW: Plugin system
│   ├── __init__.py
│   ├── base.py             # Base plugin classes
│   ├── loader.py           # Plugin loading
│   └── registry.py         # Plugin registry
├── profiling/               # NEW: Profiling
│   ├── __init__.py
│   ├── metrics.py          # Metric collection
│   └── reporter.py         # Report generation
├── logging_config.py        # NEW: Logging setup
├── count_t7.py             # Existing (with logging/progress added)
└── ... (existing files)

examples/
├── sheriff-config.yaml      # NEW: Example config
├── minimal-config.yaml      # NEW: Minimal config
└── plugins/                 # NEW: Example plugins
    └── quality_filter/

docs/
├── CLI_USAGE.md             # NEW: CLI guide
├── CONFIG_REFERENCE.md      # NEW: Config reference
└── PLUGIN_DEVELOPMENT.md    # NEW: Plugin guide
```

---

## Dependencies to Add

```toml
dependencies = [
    # Existing
    "pysam>=0.19.0",
    "typer>=0.9.0",
    "rich>=13.0.0",
    # ... existing deps

    # NEW
    "pyyaml>=6.0",           # Config files
    "pydantic>=2.0",         # Config validation
    "questionary>=2.0",      # Interactive prompts
    "psutil>=5.9.0",         # System monitoring
]

[project.optional-dependencies]
dev = [
    # Existing
    "pytest>=7.0",
    # ... existing dev deps

    # NEW
    "memory-profiler>=0.61",  # Memory profiling
]
```

---

## Success Metrics

**After Implementation:**

| Feature | Before | After |
|---------|--------|-------|
| CLI UX Score | 8/10 | 9.5/10 |
| Config Management | Manual | YAML + validation |
| Resumability | None | Full checkpoint system |
| Logging | Print statements | Structured logging |
| Profiling | Manual timing | Automated metrics |
| Extensibility | None | Plugin system |
| Interactive Setup | None | Full guided flow |
| Progress Feedback | Minimal | Rich progress bars |

**User Benefits:**
- ✅ Run Sheriff with single config file
- ✅ Resume interrupted runs without starting over
- ✅ Interactive mode for new users
- ✅ Professional logging for debugging
- ✅ Performance profiling for optimization
- ✅ Extensible via plugins
- ✅ Beautiful, informative output

---

## Backward Compatibility

All existing commands continue to work:
```bash
# Old way (still works)
sheriff get-t7-edits sample.bam ref.fa barcodes.txt genes.gtf

# New way (equivalent)
sheriff run sample.bam ref.fa barcodes.txt genes.gtf

# Or with config
sheriff run --config my-analysis.yaml
```

Zero breaking changes. Pure enhancement.
