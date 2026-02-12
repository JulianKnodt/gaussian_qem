use ordered_float::NotNan;
use pars3d::mesh::SphHarmonicCoeff;
use pars3d::quat;
use priority_queue::PriorityQueue;

use super::{
    F, add, kmul,
    manifold::CollapsibleManifold,
    quadric::{AttrWeights, GN, Quadric, QuadricAccumulator},
    sub,
};

use union_find::UnionFindOp;

use super::parameters::Args;

/// In-place simplification of planar faces of a mesh.
/// Returns how many faces were removed
pub fn simplify(
    v: &mut [[F; 3]],
    vc: &mut [[F; 3]],
    op: &mut [F],
    scale: &mut [[F; 3]],
    rot: &mut [[F; 4]],
    sphs: &mut [[SphHarmonicCoeff; 3]],
    args: &Args,
) -> usize {
    let target = match (args.number, args.ratio) {
        (None, None) => unreachable!(),
        (Some(t), None) => t,
        (None, Some(r)) => (v.len() as F * r).floor() as usize,
        (Some(_), Some(_)) => todo!(),
    };

    let get_attrs = |i: usize| {
        let mut data = [0.; GN];
        data[0] = vc[i][0];
        data[1] = vc[i][1];
        data[2] = vc[i][2];
        data[3] = op[i];

        let mut idx = 4;
        for sh in sphs[i] {
            for s in sh.flat_iter() {
                data[idx] = s;
                idx += 1;
            }
        }

        data
    };

    let mut attr_ws = AttrWeights::<GN>::default();
    if args.omit_sph {
        attr_ws.ws[4..].fill(0.);
    }
    let attr_ws = attr_ws;

    // normalize all vertices to [-1., 1]
    use std::array::from_fn;
    // Normalize the geometry of this mesh to lay in the unit box.
    let [l, h] = v
        .iter()
        .copied()
        .fold([[F::INFINITY; 3], [F::NEG_INFINITY; 3]], |[l, h], n| {
            [from_fn(|i| l[i].min(n[i])), from_fn(|i| h[i].max(n[i]))]
        });
    let center = kmul(0.5, add(l, h));
    for v in v.iter_mut() {
        *v = sub(*v, center);
    }
    let largest_val = v
        .iter()
        .copied()
        .fold(0. as F, |m, [v0, v1, v2]| m.max(v0).max(v1).max(v2));
    let pos_scale = if largest_val == 0. {
        1.
    } else {
        largest_val.recip()
    };
    for v in v.iter_mut() {
        *v = kmul(pos_scale, *v);
    }
    for s in scale.iter_mut() {
        *s = kmul(pos_scale, *s);
    }

    let mut m = CollapsibleManifold::new_with(v.len(), |vi| {
        let q = Quadric::new_gaussian(v[vi], scale[vi], rot[vi]);
        (q, v[vi], scale[vi], rot[vi])
    });

    let mut curr_costs = vec![0.; v.len()];

    macro_rules! update_cost_of_edge {
        ($e0:expr, $e1: expr) => {{
            let e0 = $e0;
            let e1 = $e1;

            let mut q_acc = QuadricAccumulator::default();
            let &(q0, _, a_s, a_r) = m.get(e0);
            let &(q1, _, b_s, b_r) = m.get(e1);
            let q01 = q0 + q1;
            q_acc += q0;
            q_acc += q1;

            let p = q_acc.point_with_volume();
            debug_assert!(p.iter().copied().all(F::is_finite));

            let attrs = q01.attributes(p, attr_ws);

            // TODO rethink this cost here?
            // maybe make it be the size of the smallest eigenvalue
            let total_cost = q01.cost_attrib(p, attrs, attr_ws).max(0.)
                - if args.no_delta_cost {
                    0.
                } else {
                    curr_costs[e0] + curr_costs[e1]
                };

            //assert!(total_cost >= 0., "{total_cost:?}");

            let (rot_cost, _, _) = best_rot(&q01.a, a_s, a_r, b_s, b_r);
            assert!(rot_cost > 0.);

            let total_cost = total_cost + args.rot_weight * rot_cost;

            unsafe { NotNan::new(-total_cost).unwrap_unchecked() }
        }};
    }

    let bt = ball_tree::BallTree::new(v.to_vec(), (0..v.len()).collect::<Vec<_>>());
    let mut q = bt.query();

    use indicatif::ProgressIterator;
    for (vi, &v) in v.iter().enumerate().progress() {
        let nbrs = q
            .nn(&v)
            .filter(|&(_, _, &adj)| adj != vi)
            .take(args.k_nearest);
        for (num, (_, dist, &adj)) in nbrs.enumerate() {
            if dist > 0.04 && num >= 15 {
                break;
            }
            m.add_edge(vi, adj);
        }
    }

    // for each vertex, compute approximate surrounding delaunay tets
    let mut nbrs = vec![];
    let mut nbr_idxs = vec![];
    let mut tets = vec![];
    let mut buf = vec![];
    use pars3d::geom_processing::delaunay::{
        BowyerWatsonSettings, SuperSimplexStrategy, bowyer_watson_3d,
    };
    let settings = BowyerWatsonSettings {
        super_simplex_strat: SuperSimplexStrategy::AABB,
        do_not_check_flips: true,
    };

    for vi in (0..v.len()).progress() {
        nbrs.clear();
        nbr_idxs.clear();
        tets.clear();

        nbrs.push(v[vi]);
        nbr_idxs.push(vi);
        for adj in m.vertex_adj(vi) {
            nbrs.push(v[adj]);
            nbr_idxs.push(adj);
        }
        let s = scale[vi].iter().copied().map(F::abs).sum::<F>();

        bowyer_watson_3d(&nbrs, &mut tets, &mut buf, settings);

        tets.retain(|t: &[usize; 4]| t.contains(&0));
        if tets.is_empty() {
            continue;
        }

        let total_vol = tets
            .iter()
            .map(|t| pars3d::signed_tet_vol(t.map(|vi| unsafe { *nbrs.get_unchecked(vi) })).abs())
            .inspect(|&v| debug_assert!(v > 0.))
            .sum::<F>();

        for &vis in &tets {
            let ps = vis.map(|vi| unsafe { *nbrs.get_unchecked(vi) });
            let tet_vol = pars3d::signed_tet_vol(ps).abs();
            assert!(tet_vol > 0., "{tet_vol}");
            let attrs = vis.map(|vi| get_attrs(unsafe { *nbr_idxs.get_unchecked(vi) }));
            let qa = Quadric::tet_attribs(ps, attrs, attr_ws);

            unsafe { m.data.get_unchecked_mut(vi) }.0 += qa * s * (tet_vol / total_vol);
        }
    }

    let mut pq = PriorityQueue::new();
    let it = m
        .edges_ord()
        .map(|[e0, e1]| ([e0, e1], update_cost_of_edge!(e0, e1)));
    pq.extend(it);

    //let mut did_update = vec![];
    let pbar = indicatif::ProgressBar::new(m.vertices.len() as u64);
    let mut num_hit = 0;
    while let Some(([e0, e1], _qem)) = pq.pop() {
        assert!(e0 < e1);
        if m.is_deleted(e0) || m.is_deleted(e1) {
            continue;
        }
        let l = m.vertices.len();
        pbar.set_position(l as u64);
        if l <= target {
            break;
        }

        let mut q_acc = QuadricAccumulator::default();
        let q0 = m.get(e0).0;
        let q1 = m.get(e1).0;
        q_acc += q0;
        q_acc += q1;
        let pos = q_acc.point_with_volume();
        let q01 = q0 + q1;

        // -- Commit

        m.merge(e0, e1, |(_, _, a_s, a_r), (_, _, b_s, b_r)| {
            curr_costs[e1] = q01
                .cost_attrib(pos, q01.attributes(pos, attr_ws), attr_ws)
                .max(0.);

            let (_, best_s, best_r) = best_rot(&q01.a, *a_s, *a_r, *b_s, *b_r);
            (q01, pos, best_s, best_r)
        });
        debug_assert!(m.is_deleted(e0));
        debug_assert!(!m.is_deleted(e1));

        debug_assert_eq!(m.get_new_vertex(e1), e1);
        for adj in m.vertex_adj(e1) {
            let adj_e @ [l, h] = std::cmp::minmax(e1, adj);
            pq.push(adj_e, update_cost_of_edge!(l, h));
        }

        num_hit += 1;
        if num_hit == 400_000 {
            num_hit = 0;
            pq.retain(|&[e0, e1], _| !m.is_deleted(e0) && !m.is_deleted(e1));
        }
    }

    //const BASES: [[F; 3]; 3] = [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]];

    for (vi, (_, &(q, p, s, r))) in m.vertices().enumerate() {
        assert!(p.into_iter().all(F::is_finite));
        v[vi] = p;
        let attrs = q.attributes(p, attr_ws);

        let mut set_attrs = |i: usize, attrs: [F; GN]| {
            assert!(attrs.into_iter().all(F::is_finite));
            vc[i][0] = attrs[0];
            vc[i][1] = attrs[1];
            vc[i][2] = attrs[2];
            op[i] = attrs[3];

            let mut idx = 4;
            for sh in &mut sphs[i] {
                for s in sh.flat_iter_mut() {
                    *s = attrs[idx];
                    idx += 1;
                }
            }
        };

        set_attrs(vi, attrs);

        // eigenvalues, eigenvectors (scale, basis)
        rot[vi] = r;
        scale[vi] = s;
    }

    // denormalize all output vertices
    let inv_pos_scale = pos_scale.recip();
    for v in v.iter_mut() {
        *v = add(kmul(inv_pos_scale, *v), center);
    }
    for s in scale.iter_mut() {
        *s = kmul(inv_pos_scale, *s);
    }

    m.vertices.len()
}

fn covariance([sx, sy, sz]: [F; 3], q: [F; 4]) -> [[F; 3]; 3] {
    let l = [[sx, 0., 0.], [0., sy, 0.], [0., 0., sz]];
    let r = quat::quat_to_mat(q);
    let rl = pars3d::matmul(r, l);
    let cov = pars3d::matmul(rl, pars3d::transpose(rl));
    assert_eq!(cov[0][1], cov[1][0]);
    assert_eq!(cov[0][2], cov[2][0]);
    assert_eq!(cov[1][2], cov[2][1]);
    cov
}

pub fn frobenius(a: [[F; 3]; 3], b: [[F; 3]; 3]) -> F {
    let mut sum = 0.;
    for i in 0..3 {
        for j in 0..3 {
            //sum += a[i][j] * b[i][j];
            sum += (a[i][j] - b[i][j]).abs()
        }
    }
    return sum;
}

pub fn best_rot(
    q: &super::sym::SymMatrix3,
    a_s: [F; 3],
    a_r: [F; 4],
    b_s: [F; 3],
    b_r: [F; 4],
) -> (F, [F; 3], [F; 4]) {
    let a_rm = covariance(a_s, a_r);
    let b_rm = covariance(b_s, b_r);

    let (es, vs) = q.eigen();

    let mut best_v = F::INFINITY;
    let mut best_r = [0.; 4];
    let mut best_s = [0.; 3];

    const IJKS: [[usize; 3]; 6] = [
        [0, 1, 2],
        [0, 2, 1],
        [1, 0, 2],
        [1, 2, 0],
        [2, 0, 1],
        [2, 1, 0],
    ];
    for axes in IJKS {
        assert_ne!(axes[0], axes[1]);
        assert_ne!(axes[1], axes[2]);
        assert_ne!(axes[0], axes[2]);

        for f in 0..8 {
            let mat = std::array::from_fn(|i| {
                let is_neg = ((f >> i) & 1) == 1;
                vs[axes[i]].map(|v| if is_neg { -v } else { v })
            });
            let s = axes.map(|a| es[a]);
            let rot = quat::quat_from_mat(mat);
            let cov = covariance(s, rot);
            let v = frobenius(cov, a_rm) + frobenius(cov, b_rm);

            if v < best_v {
                best_v = v;
                best_s = s;
                best_r = rot;
            }
        }
    }

    (best_v, best_s, best_r)
}
