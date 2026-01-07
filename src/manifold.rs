use union_find::{UnionFind, UnionFindOp};

/// A mesh representation which is suitable for collapsing vertices.
/// It can associate data with each vertex, and each edge.
/// Associated edge data is oriented.
#[derive(Debug, Clone)]
pub struct CollapsibleManifold<T> {
    pub(crate) vertices: UnionFind<u32>,

    pub edges: Vec<Vec<u32>>,

    pub data: Vec<T>,
}

impl<T> CollapsibleManifold<T> {
    pub fn new_with(size: usize, f: impl Fn(usize) -> T) -> Self {
        let mut data = Vec::with_capacity(size);
        for i in 0..size {
            data.push(f(i));
        }
        Self {
            vertices: UnionFind::new_u32(size),
            edges: vec![vec![]; size],

            data,
        }
    }
}

impl<T> CollapsibleManifold<T> {
    pub fn get_new_vertex(&self, old: usize) -> usize {
        self.vertices.find(old)
    }
    pub fn vertices(&self) -> impl Iterator<Item = (usize, &T)> + '_ {
        (0..self.vertices.capacity())
            .filter(|&vi| !self.is_deleted(vi))
            .map(|vi| (vi, &self.data[vi]))
    }
    pub fn is_deleted(&self, vi: usize) -> bool {
        !self.vertices.is_root(vi)
    }

    /// Adds an edge. For faces, should call `add_face`.
    pub fn add_edge(&mut self, v0: usize, v1: usize) {
        if v0 == v1 {
            return;
        }
        if !self.edges[v0].contains(&(v1 as u32)) {
            self.edges[v0].push(v1 as u32);
        }
        if !self.edges[v1].contains(&(v0 as u32)) {
            self.edges[v1].push(v0 as u32);
        }

        // note that this is not using the mapping since edges should only be added ahead of
        // time.
        self.edges[v0].sort_unstable_by_key(|&dst| dst);
        self.edges[v1].sort_unstable_by_key(|&dst| dst);

        self.edges[v0].dedup_by_key(|&mut dst| dst);
        self.edges[v1].dedup_by_key(|&mut dst| dst);
    }

    /// Returns adjacent vertices (should always be in sorted order)
    pub fn vertex_adj(&self, v: usize) -> impl Iterator<Item = usize> + '_ {
        self.edges[v].iter().map(|&dst| {
            assert!(!self.is_deleted(dst as usize));
            self.vertices.find(dst as usize)
        })
    }

    /// Returns whether two vertices v0 and v1 are adjacent.
    /// v0 and v1 can be merged into other vertices.
    #[inline]
    pub fn is_adj(&self, v0: usize, v1: usize) -> bool {
        let v0 = self.vertices.find(v0);
        let v1 = self.vertices.find(v1);
        self.edges[v0]
            .iter()
            .any(|&dst| self.vertices.find(dst as usize) == v1)
    }

    /// Merges v0 into v1.
    pub fn merge(&mut self, src: usize, dst: usize, mut merge: impl FnMut(&T, &T) -> T)
    where
        T: Clone,
    {
        assert!(src < dst);
        debug_assert!(!self.is_deleted(src));
        debug_assert!(!self.is_deleted(dst));
        debug_assert!(self.is_adj(src, dst));

        self.vertices.union(src, dst);

        let [data_src, data_dst] = unsafe { self.data.get_disjoint_unchecked_mut([src, dst]) };
        let new_data = merge(data_src, data_dst);
        // data_src should no longer be accessed
        //*data_src = new_data.clone();
        *data_dst = new_data;

        let [src_e, dst_e] = unsafe { self.edges.get_disjoint_unchecked_mut([src, dst]) };
        let pos = dst_e.iter().position(|&v| v == src as u32).unwrap();
        assert_eq!(src as u32, dst_e.swap_remove(pos));

        for v in dst_e.iter_mut() {
            *v = self.vertices.find(*v as usize) as u32;
        }
        let mut src_e = std::mem::take(src_e);
        let pos = src_e.iter().position(|&v| v == dst as u32).unwrap();
        assert_eq!(dst as u32, src_e.swap_remove(pos));

        let curr_dst_len = dst_e.len();
        for v in src_e {
            let v = self.vertices.find(v as usize) as u32;
            if !dst_e[0..curr_dst_len].contains(&v) {
                dst_e.push(v);
            }
        }

        let tmp = std::mem::take(&mut self.edges[dst]);
        for &adj in &tmp {
            let adj = adj as usize;
            assert_eq!(self.vertices.find(adj), adj);

            let adj_e = unsafe { self.edges.get_unchecked_mut(adj) };
            let mut i = 0;
            while i < adj_e.len() {
                if self.vertices.find(adj_e[i] as usize) == dst {
                    adj_e.swap_remove(i);
                } else {
                    i += 1;
                }
            }
            adj_e.push(dst as u32);
        }
        self.edges[dst] = tmp;
    }

    pub fn get(&self, v: usize) -> &T {
        unsafe { self.data.get_unchecked(self.vertices.find(v)) }
    }

    /// All edges in this manifold mesh with v0-v1 in sorted order.
    pub fn edges_ord(&self) -> impl Iterator<Item = [usize; 2]> + '_ {
        self.edges.iter().enumerate().flat_map(|(src, dsts)| {
            dsts.iter()
                .filter(move |&&dst| src < dst as usize)
                .map(move |&dst| [src, dst as usize])
        })
    }
}
