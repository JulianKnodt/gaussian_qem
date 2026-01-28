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
    /// Input mesh file.
    #[arg(short, long, required = true)]
    pub input: String,

    /// Output mesh file.
    #[arg(short, long, required = true)]
    pub output: String,

    /// Output stat file.
    #[arg(long)]
    pub stats: String,

    /// Approximate ratio of output/input tris.
    #[arg(short = 'r', long, group = "target")]
    pub ratio: Option<F>,

    /// Approximate number of output tris.
    #[arg(short = 'n', long, group = "target")]
    pub number: Option<usize>,

    /// Do not use delta cost.
    #[arg(long)]
    pub no_delta_cost: bool,

    /// Do not include spherical harmonics during decimation.
    #[arg(long)]
    pub omit_sph: bool,

    /// Number of neighbors to use for sampling
    #[arg(long, default_value_t = 20)]
    pub k_nearest: usize,
}
