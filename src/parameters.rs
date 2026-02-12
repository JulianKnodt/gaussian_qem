pub use super::F;
use clap::Parser;

/// Mesh geometry decimation.
#[derive(Parser, Default, Debug)]
#[clap(group(
            clap::ArgGroup::new("target")
                .required(true)
                .args(&["ratio", "number"]),
        ))]
pub struct Args {
    /// Input Splat file.
    #[arg(short, long, required = true)]
    pub input: String,

    /// Output Splat file.
    #[arg(short, long, required = true)]
    pub output: String,

    /// Output stat file of run details.
    #[arg(long, default_value_t = String::new())]
    pub stats: String,

    /// Ratio of output splats to input splats (i.e. 0.5 = 50%)
    #[arg(short = 'r', long, group = "target")]
    pub ratio: Option<F>,

    /// Number of output splats
    #[arg(short = 'n', long, group = "target")]
    pub number: Option<usize>,

    /// Do not use delta cost (probably don't need to change).
    #[arg(long)]
    pub no_delta_cost: bool,

    /// Do not include spherical harmonics during decimation.
    #[arg(long)]
    pub omit_sph: bool,

    /// Number of neighbors to use for sampling
    #[arg(long, default_value_t = 20)]
    pub k_nearest: usize,

    /// [DEBUG] Test I/O without simplification
    #[arg(long, hide = true)]
    pub test_io: bool,

    /// Weight to use for preserving orientation of gaussians.
    /// (Probably doesn't need to be changed).
    #[arg(long, default_value_t = 1e-4)]
    pub rot_weight: F,

    /// How often to clean up the priority queue.
    /// Only affects run-time (probably doesn't need to be changed).
    #[arg(long, default_value_t = 200_000)]
    pub clean_pq_freq: usize,
}
