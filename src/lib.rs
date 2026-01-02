#![allow(incomplete_features)]
#![feature(generic_const_exprs)]
#![feature(cmp_minmax)]

pub type F = pars3d::F;

pub mod quadric;
pub mod svd;
pub mod sym;

pub mod vec;
use vec::*;

mod manifold;

mod parameters;
pub use parameters::Args;

mod qem;
pub use qem::simplify;
