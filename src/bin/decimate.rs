use clap::Parser;

use gaussian_qem::{Args, simplify};

pub fn main() {
    let args = Args::parse();
    let mut scene = pars3d::load(&args.input).expect(&format!("Failed to open {}", args.input));
    let mut m = scene.into_flattened_mesh();
    // just the vertices and faces, fuse together identical positions
    let va = &mut m.vertex_attrs;
    simplify(
        &mut m.v,
        &mut m.vert_colors,
        &mut va.opacity,
        &mut va.scale,
        &mut va.rot,
        &mut va.sph_harmonic_coeff,
        &args,
    );

    m.repopulate_scene(&mut scene);
    pars3d::save(&args.output, &scene).expect(&format!("Failed to save to {}", args.output));
}
