use ordered_float::NotNan;
use pars3d::mesh::SphHarmonicCoeff;
use priority_queue::PriorityQueue;

use super::{
    F, add, kmul,
    manifold::CollapsibleManifold,
    quadric::{AttrWeights, GN, Quadric, QuadricAccumulator},
    sub,
};

use std::ops::Neg;

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

    let attr_ws = AttrWeights::<GN>::default();
    //attr_ws.ws[0..3].fill(1e-2);
    //attr_ws.ws[3] = 1e-2;

    /*{
        // TODO temp testing block
        let is = [0,1,2,3];
        let qs = is.map(|i| Quadric::<GN>::new_gaussian(v[i], scale[i], rot[i]));
        let attrs = is.map(get_attrs);
        let qat = Quadric::tet_attribs(is.map(|i| v[i]), attrs, attr_ws);
        let mut qa = QuadricAccumulator::default();
        let q0 = qs[0] + qat * scale[0].iter().copied().map(F::abs).sum::<F>();
        qa += q0;
        /*
        println!("{:?} {:?}", qa.point_with_volume(), v[0]);
        let out_attribs = q0.attributes(v[0], attr_ws);
        println!("{:?}", &out_attribs[0..4]);
        println!("{:?}", &attrs[0][0..4]);
        println!();
        */
        use pars3d::length;

        let (evalues, [v0, v1, v2]) = q0.a.eigen_sorted();
        println!("new {:?}\nold {:?}", evalues, scale[0]);

        let quat_rot = pars3d::quat::quat_from_standard(v0.map(Neg::neg), v1.map(Neg::neg));
        println!("rot\nnew: {quat_rot:?}\nog : {:?}", rot[0]);
        println!();
        todo!();
    }*/

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

    let mut m = CollapsibleManifold::new_with(v.len(), |vi| {
        let q = Quadric::new_gaussian(v[vi], scale[vi], rot[vi]);
        (q, v[vi])
    });

    let mut curr_costs = vec![0.; v.len()];

    macro_rules! update_cost_of_edge {
        ($e0:expr, $e1: expr) => {{
            let e0 = $e0;
            let e1 = $e1;

            let mut q_acc = QuadricAccumulator::default();
            let &(q0, _) = m.get(e0);
            let &(q1, _) = m.get(e1);
            let q01 = q0 + q1;
            q_acc += q0;
            q_acc += q1;

            let p = q_acc.point_with_volume();
            debug_assert!(p.iter().copied().all(F::is_finite));

            let attrs = q01.attributes(p, attr_ws);

            let total_cost = q01.cost_attrib(p, attrs, attr_ws).max(0.)
                - if args.no_delta_cost {
                    0.
                } else {
                    curr_costs[e0] + curr_costs[e1]
                };

            unsafe { NotNan::new(-total_cost).unwrap_unchecked() }
        }};
    }

    /*
    {   // TODO temp testing block
        let mut qa = QuadricAccumulator::default();
        let q0 = m.data[0].0;
        qa += q0;
        println!("{:?} {:?}", qa.point_with_volume(), v[0]);
        todo!();
    }
    */

    let bt = ball_tree::BallTree::new(v.to_vec(), (0..v.len()).collect::<Vec<_>>());
    let mut q = bt.query();

    use indicatif::ProgressIterator;
    for (vi, &v) in v.iter().enumerate().progress() {
        let nbrs = q.nn(&v).filter(|&(_, _, &adj)| adj != vi).take(16);
        for (num, (_, dist, &adj)) in nbrs.enumerate() {
            if dist > 0.02 && num >= 12 {
                break;
            }
            let [e0, e1] = std::cmp::minmax(vi, adj);
            m.add_edge(e0, e1);
        }
        //break; // tmp
    }

    // for each vertex, compute delaunay tets
    timed_regions::TimerStruct! {
      struct DelaunayTimer {
        delaunay,
        compute1,
        compute2,
      }
    };

    let mut timer = DelaunayTimer::default();

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

        timer.delaunay.start();
        bowyer_watson_3d(&nbrs, &mut tets, &mut buf, settings);
        timer.delaunay.stop();

        tets.retain(|t: &[usize; 4]| t.contains(&0));
        assert!(!tets.is_empty());

        let total_vol = tets
            .iter()
            .map(|t| pars3d::signed_tet_vol(t.map(|vi| unsafe { *nbrs.get_unchecked(vi) })).abs())
            .inspect(|&v| debug_assert!(v > 0.))
            .sum::<F>();

        timer.compute1.start();
        for &vis in &tets {
            let ps = vis.map(|vi| unsafe { *nbrs.get_unchecked(vi) });
            let tet_vol = pars3d::signed_tet_vol(ps).abs();
            timer.compute2.start();
            let attrs = vis.map(|vi| get_attrs(unsafe { *nbr_idxs.get_unchecked(vi) }));
            timer.compute2.stop();
            let qa = Quadric::tet_attribs(ps, attrs, attr_ws);
            m.data[vi].0 += qa * s * (tet_vol / total_vol);
        }
        timer.compute1.stop();

        /*{
            // TODO temp testing block
            let attr0 = get_attrs(0);
            let mut qa = QuadricAccumulator::default();
            let q0 = m.data[0].0;
            qa += q0;
            println!("{:?} {:?}", qa.point_with_volume(), v[0]);
            let out_attribs = q0.attributes(v[0], attr_ws);
            println!("new attr {:?}", &out_attribs[0..4]);
            println!("original attr {:?}", &attr0[0..4]);

            let (evalues, [v0, v1, _v2]) = q0.a.eigen();
            let quat_rot = pars3d::quat::quat_from_standard(v0, v1);
            println!("scale {evalues:?} {:?}", scale[0]);
            println!("rot {quat_rot:?} {:?}", rot[0]);
            todo!();
        }*/
    }

    let mut pq = PriorityQueue::new();
    let it = m
        .edges_ord()
        .map(|[e0, e1]| ([e0, e1], update_cost_of_edge!(e0, e1)));
    pq.extend(it);

    timer.print();

    timed_regions::TimerStruct! {
      struct Timer {
        quadric_acc,
        merge,
        update0,
        update1,
        pq,
      }
    };
    let mut timer = Timer::default();

    //let mut did_update = vec![];
    let pbar = indicatif::ProgressBar::new(m.vertices.len() as u64);
    timer.pq.start();
    while let Some(([e0, e1], _qem)) = pq.pop() {
        timer.pq.stop();
        assert!(e0 < e1);
        if m.is_deleted(e0) || m.is_deleted(e1) {
            timer.pq.start();
            continue;
        }
        let l = m.vertices.len();
        pbar.set_position(l as u64);
        if l <= target {
            break;
        }

        timer.quadric_acc.start();
        let mut q_acc = QuadricAccumulator::default();
        let q0 = m.get(e0).0;
        let q1 = m.get(e1).0;
        q_acc += q0;
        q_acc += q1;
        let pos = q_acc.point_with_volume();
        let q01 = q0 + q1;
        timer.quadric_acc.stop();

        // -- Commit

        timer.merge.start();
        m.merge(e0, e1, |_, _| {
            curr_costs[e1] = q01
                .cost_attrib(pos, q01.attributes(pos, attr_ws), attr_ws)
                .max(0.);
            (q01, pos)
        });
        debug_assert!(m.is_deleted(e0));
        debug_assert!(!m.is_deleted(e1));
        timer.merge.stop();

        timer.update0.start();
        debug_assert_eq!(m.get_new_vertex(e1), e1);
        for adj in m.vertex_adj(e1) {
            let adj_e @ [l, h] = std::cmp::minmax(e1, adj);
            let prio = update_cost_of_edge!(l, h);
            pq.push(adj_e, prio);
        }
        timer.update0.stop();

        timer.pq.start();
    }

    timer.print();

    for (vi, &(q, p)) in m.vertices() {
        op[vi] = m.merged_vertices(vi).map(|og| op[og]).sum::<F>();
        /*
           .max_by(F::total_cmp)
           .unwrap();
        */

        assert!(p.into_iter().all(F::is_finite));
        v[vi] = p;
        let attrs = q.attributes(p, attr_ws);

        let mut set_attrs = |i: usize, attrs: [F; GN]| {
            vc[i][0] = attrs[0];
            vc[i][1] = attrs[1];
            vc[i][2] = attrs[2];
            //op[i] = attrs[3];

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
        let (es, [v0, v1, _v2]) = q.a.eigen_sorted();
        let quat_rot = pars3d::quat::quat_from_standard(v0.map(Neg::neg), v1.map(Neg::neg));
        rot[vi] = quat_rot;
        scale[vi] = es;
    }

    // denormalize all output vertices
    let inv_pos_scale = pos_scale.recip();
    for v in v.iter_mut() {
        *v = add(kmul(inv_pos_scale, *v), center);
    }

    m.vertices.len()
}
