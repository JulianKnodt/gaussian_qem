use std::collections::HashMap;
use std::collections::hash_map::Entry;

use std::io::Write;
use std::time::Instant;

use clap::Parser;

use gaussian_qem::{Args, simplify};
use pars3d::F;

pub fn main() -> std::io::Result<()> {
    let args = Args::parse();
    assert!(args.output.ends_with(".ply"), "Only PLY output supported");
    let mut scene = pars3d::load(&args.input).expect(&format!("Failed to open {}", args.input));
    assert_eq!(scene.meshes.len(), 1);
    let mut m = scene.meshes.pop().unwrap();
    m.n.clear();

    println!("[INFO]: Input ({}) has #{} vertices", args.input, m.v.len());
    // just the vertices and faces, fuse together identical positions
    let va = &mut m.vertex_attrs;

    let rmed_verts = dedup(&mut m.v, &mut va.opacity, |_dst, src, verts, opacity| {
        verts.swap_remove(src);
        opacity.swap_remove(src);
        m.vert_colors.swap_remove(src);
        va.scale.swap_remove(src);
        va.rot.swap_remove(src);
        va.sph_harmonic_coeff.swap_remove(src);
    });
    if rmed_verts != 0 {
        println!("[INFO]: Deduped {rmed_verts} splats.");
    }

    for s in va.scale.iter_mut() {
        // switch from log-scale to linear
        *s = (*s).map(F::exp);
    }

    let start = Instant::now();
    // sort all the scale and rotations of the input
    let num_verts = if args.test_io {
        m.v.len()
    } else {
        simplify(
            &mut m.v,
            &mut m.vert_colors,
            &mut va.opacity,
            &mut va.scale,
            &mut va.rot,
            &mut va.sph_harmonic_coeff,
            &args,
        )
    };

    let elapsed = start.elapsed();

    m.v.truncate(num_verts);
    m.vert_colors.truncate(num_verts);
    va.truncate(num_verts);

    for s in va.scale.iter_mut() {
        *s = (*s).map(F::ln);
    }

    println!(
        "[INFO]: Output ({}) has #{} vertices",
        args.output,
        m.v.len()
    );

    let p: pars3d::ply::Ply = m.into();
    let out = std::fs::File::create(&args.output)?;
    let out = std::io::BufWriter::new(out);
    //p.write(out, pars3d::ply::FormatKind::Ascii)?;
    p.write(out, pars3d::ply::FormatKind::BinLil)?;
    //.expect(&format!("Failed to save to {}", args.output));

    if !args.stats.is_empty() {
        let mut f = std::fs::File::create(args.stats)?;
        write!(
            f,
            r#"{{
  "elapsed_secs": {}
}}"#,
            elapsed.as_secs_f32()
        )?;
    }

    Ok(())
}

pub fn dedup(
    vs: &mut Vec<[F; 3]>,
    opacity: &mut Vec<F>,
    mut merge: impl FnMut(usize, usize, &mut Vec<[F; 3]>, &mut Vec<F>),
) -> usize {
    let mut all_verts = HashMap::new();
    let mut rmed = 0;
    let mut vi = 0;
    while vi < vs.len() {
        let v = vs[vi];
        let key = v.map(F::to_bits);
        match all_verts.entry(key) {
            Entry::Vacant(v) => {
                v.insert(vi);
                vi += 1;
            }
            Entry::Occupied(o) => {
                merge(*o.get(), vi, vs, opacity);
                rmed += 1;
            }
        }
    }

    rmed
}
