use clap::Parser;

use gaussian_qem::{Args, simplify};
use pars3d::F;

use pars3d::quat::{quat_from_standard, quat_rot};

pub fn main() -> std::io::Result<()> {
    let args = Args::parse();
    assert!(args.output.ends_with(".ply"), "Only PLY output supported");
    let mut scene = pars3d::load(&args.input).expect(&format!("Failed to open {}", args.input));
    assert_eq!(scene.meshes.len(), 1);
    let mut m = scene.meshes.pop().unwrap();
    println!("[INFO]: Input has #{} vertices", m.v.len());
    // just the vertices and faces, fuse together identical positions
    let va = &mut m.vertex_attrs;

    /*
    for o in va.opacity.iter_mut() {
        *o = sigmoid(*o);
    }
    */

    let mut shuffle_order = vec![[0, 1, 2]; m.v.len()];

    for (vi, s) in va.scale.iter_mut().enumerate() {
        // switch from log-scale to linear
        *s = (*s).map(F::exp);

        // sort scales
        let rot = va.rot[vi];
        let mut bases = [
            quat_rot([1., 0., 0.], rot),
            quat_rot([0., 1., 0.], rot),
            quat_rot([0., 0., 1.], rot),
        ];

        for [i, j] in [[0, 1], [1, 2], [0, 1]] {
            if s[i] > s[j] {
                s.swap(i, j);
                bases.swap(i, j);
                shuffle_order[vi].swap(i, j);
            }
        }
        assert!(0. < s[0]);
        assert!(s[0] < s[1]);
        assert!(s[1] < s[2]);

        va.rot[vi] = quat_from_standard(bases[0], bases[1]);
    }

    // sort all the scale and rotations of the input

    let num_verts = simplify(
        &mut m.v,
        &mut m.vert_colors,
        &mut va.opacity,
        &mut va.scale,
        &mut va.rot,
        &mut va.sph_harmonic_coeff,
        &args,
    );

    m.v.truncate(num_verts);
    m.vert_colors.truncate(num_verts);

    va.opacity.truncate(num_verts);
    va.scale.truncate(num_verts);
    va.rot.truncate(num_verts);
    va.sph_harmonic_coeff.truncate(num_verts);

    for (vi, s) in va.scale.iter_mut().enumerate() {
        *s = (*s).map(F::ln);

        let mut out = [0.; 3];
        for (src, &dst) in shuffle_order[vi].iter().enumerate() {
            out[dst] = s[src];
        }
        *s = out;

        let rot = va.rot[vi];
        let bases = [
            quat_rot([1., 0., 0.], rot),
            quat_rot([0., 1., 0.], rot),
            quat_rot([0., 0., 1.], rot),
        ];
        let mut new_bases = [[0.; 3]; 3];
        for (src, &dst) in shuffle_order[vi].iter().enumerate() {
            new_bases[dst] = bases[src];
        }

        va.rot[vi] = quat_from_standard(new_bases[0], new_bases[1]);
    }

    /*
    for o in va.opacity.iter_mut() {
        *o = inv_sigmoid(*o);
    }
    */

    let p: pars3d::ply::Ply = m.into();
    let out = std::fs::File::create(&args.output)?;
    let out = std::io::BufWriter::new(out);
    //p.write(out, pars3d::ply::FormatKind::Ascii)
    p.write(out, pars3d::ply::FormatKind::BinLil)
        .expect(&format!("Failed to save to {}", args.output));
    Ok(())
}

pub fn sigmoid(x: F) -> F {
    1. / (1. + (-x).exp())
}

pub fn inv_sigmoid(y: F) -> F {
    (y / (1. - y)).ln()
}
/*
*/
