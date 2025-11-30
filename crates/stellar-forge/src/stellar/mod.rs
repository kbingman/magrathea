pub(crate) mod sampling;
pub(crate) mod spectral_class;
pub(crate) mod stellar_color;
pub(crate) mod stellar_object;

#[cfg(test)]
mod stellar_color_test;
#[cfg(test)]
mod stellar_object_test;

pub use stellar_object::{
    BlackHole, GiantStar, MainSequenceStar, NeutronStar, StellarObject, WhiteDwarf,
};
