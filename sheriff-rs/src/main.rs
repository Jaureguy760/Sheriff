use sheriff_rs::bam_filter::{filter_bam_by_barcodes, load_whitelist};
use clap::Parser;
use anyhow::Result;

#[derive(Parser)]
#[command(name = "sheriff-rs")]
#[command(about = "High-performance BAM filtering for Sheriff", long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Parser)]
enum Commands {
    /// Filter BAM file by cell barcode whitelist
    FilterBam {
        /// Input BAM file path
        #[arg(short, long)]
        input: String,

        /// Output BAM file path
        #[arg(short, long)]
        output: String,

        /// Cell barcode whitelist file (one barcode per line)
        #[arg(short, long)]
        whitelist: String,
    },
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    match cli.command {
        Commands::FilterBam { input, output, whitelist } => {
            eprintln!("Sheriff-rs BAM Filter");
            eprintln!("=====================");
            eprintln!();
            eprintln!("Loading whitelist from: {}", whitelist);

            let barcodes = load_whitelist(&whitelist)?;
            eprintln!("  ✓ Loaded {} barcodes", barcodes.len());

            eprintln!();
            eprintln!("Filtering BAM: {} -> {}", input, output);

            let start = std::time::Instant::now();
            let result = filter_bam_by_barcodes(&input, &output, &barcodes)?;
            let duration = start.elapsed();

            eprintln!("  ✓ Complete in {:.3}s", duration.as_secs_f64());
            eprintln!();

            // Output JSON for benchmarking
            let json = serde_json::json!({
                "reads_processed": result.reads_processed,
                "reads_kept": result.reads_kept,
                "reads_rejected": result.reads_rejected,
                "duration_seconds": duration.as_secs_f64(),
                "throughput_reads_per_sec": result.reads_processed as f64 / duration.as_secs_f64(),
            });
            println!("{}", serde_json::to_string_pretty(&json)?);

            eprintln!();
            eprintln!("Statistics:");
            eprintln!("  Reads processed: {:>12}", format_number(result.reads_processed));
            eprintln!("  Reads kept:      {:>12} ({:.1}%)",
                     format_number(result.reads_kept),
                     100.0 * result.reads_kept as f64 / result.reads_processed as f64);
            eprintln!("  Reads rejected:  {:>12} ({:.1}%)",
                     format_number(result.reads_rejected),
                     100.0 * result.reads_rejected as f64 / result.reads_processed as f64);
            eprintln!("  Throughput:      {:>12} reads/sec",
                     format_number((result.reads_processed as f64 / duration.as_secs_f64()) as usize));
        }
    }

    Ok(())
}

fn format_number(n: usize) -> String {
    n.to_string()
        .as_bytes()
        .rchunks(3)
        .rev()
        .map(std::str::from_utf8)
        .collect::<Result<Vec<&str>, _>>()
        .unwrap()
        .join(",")
}
