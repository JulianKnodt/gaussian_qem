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

    /// Approximate ratio of output/input tris.
    #[arg(short = 'r', long, group = "target")]
    pub ratio: Option<F>,

    /// Approximate number of output tris.
    #[arg(short = 'n', long, group = "target")]
    pub number: Option<usize>,

    /// Do not use delta cost.
    #[arg(long)]
    pub no_delta_cost: bool,
}
