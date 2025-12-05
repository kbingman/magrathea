pub mod composition;
pub mod generation;
pub mod metadata;
pub mod planet;
pub mod planet_class;
pub mod planet_type;
pub mod sampling;
pub mod system;

// Re-export key generation functions at crate root
pub use generation::{
    from_star, from_star_with_id, generate_planetary_system, generate_planetary_system_named,
    generate_planetary_system_random,
};

// Re-export metadata types
pub use metadata::{GenerationMethod, SystemMetadata};

// Re-export planet types
pub use planet::HostStar;

#[cfg(test)]
pub mod composition_test;
#[cfg(test)]
pub mod generation_test;
#[cfg(test)]
pub mod planet_class_test;
#[cfg(test)]
pub mod planet_test;
#[cfg(test)]
pub mod planet_type_test;
#[cfg(test)]
pub mod sampling_test;
#[cfg(test)]
pub mod system_test;
