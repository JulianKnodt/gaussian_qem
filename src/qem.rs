use ordered_float::NotNan;
use pars3d::mesh::SphHarmonicCoeff;
use priority_queue::PriorityQueue;

use super::{
    F, add, kmul,
    manifold::CollapsibleManifold,
    quadric::{AttrWeights, Quadric, QuadricAccumulator},
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

    let attr_ws = AttrWeights::<49>::default();
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
            let [e0, e1] = std::cmp::minmax($e1, $e0);
            let mut q_acc = QuadricAccumulator::default();
            let &(q0, _) = m.get(e0);
            q_acc += q0;
            let &(q1, _) = m.get(e1);
            q_acc += q1;
            let p = q_acc.point_with_volume();
            debug_assert!(p.iter().copied().all(F::is_finite));

            let q01 = q0 + q1;
            let attrs = q01.attributes(p, attr_ws);

            let total_cost =
                (q0 + q1).cost_attrib(p, attrs, attr_ws).max(0.) - curr_costs[e0] - curr_costs[e1];

            NotNan::new(-total_cost).unwrap()
        }};
    }

    let bt = ball_tree::BallTree::new(v.to_vec(), (0..v.len()).collect::<Vec<_>>());
    let mut q = bt.query();

    let mut pq = PriorityQueue::new();
    for (vi, &v) in v.iter().enumerate() {
        for (_, dist, &adj) in q.nn(&v).take(50) {
            if dist > 0.01 {
                break;
            }
            let [e0, e1] = std::cmp::minmax(vi, adj);
            pq.push([e0, e1], update_cost_of_edge!(e0, e1));
            m.add_edge(e0, e1);
        }
    }
    drop(bt);

    let mut did_update = vec![];
    while let Some(([e0, e1], _qem)) = pq.pop() {
        debug_assert!(e0 < e1);
        if m.is_deleted(e0) || m.is_deleted(e1) {
            continue;
        }
        if m.vertices.len() <= target {
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

        m.merge(e0, e1, |_, _| {
            curr_costs[e1] = q01
                .cost_attrib(pos, q01.attributes(pos, attr_ws), attr_ws)
                .max(0.);
            (q01, pos)
        });
        assert!(m.is_deleted(e0));
        assert!(!m.is_deleted(e1));

        did_update.clear();
        let e_dst = m.get_new_vertex(e1);
        for adj in m.vertex_adj(e_dst) {
            let prio = update_cost_of_edge!(e_dst, adj);
            let adj_e = std::cmp::minmax(e_dst, adj);
            pq.push(adj_e, prio);
            did_update.push(adj_e);
        }

        for adj in m.vertex_adj(e_dst) {
            for adj2 in m.vertex_adj(adj) {
                let adj_e = std::cmp::minmax(adj, adj2);
                if adj2 == e_dst || did_update.contains(&adj_e) {
                    continue;
                }

                did_update.push(adj_e);
                let prio = update_cost_of_edge!(adj, adj2);
                pq.push(adj_e, prio);
                continue;
            }
        }
    }

    for (vi, &(_, p)) in m.vertices() {
        v[vi] = p;
    }

    // denormalize all output vertices
    let inv_pos_scale = pos_scale.recip();
    for v in v.iter_mut() {
        *v = add(kmul(inv_pos_scale, *v), center);
    }

    m.vertices.len()
}
