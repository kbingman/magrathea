//! Stochastic variation for planetary properties
//!
//! Provides deterministic pseudo-random variation to create natural diversity
//! in planetary properties while maintaining reproducibility.
//!
//! The key insight is that the same planet (defined by mass and orbital distance)
//! should always produce the same "random" variations, making results reproducible
//! across runs.

use std::hash::{Hash, Hasher};

use rand::Rng;
use rand_chacha::ChaCha8Rng;

/// Provides stochastic variation to planetary properties to create natural diversity
///
/// Uses deterministic seeding based on planet properties so the same planet
/// always gets the same variation (reproducible results).
///
/// # Example
/// ```
/// use planetary::variation::PlanetaryVariation;
///
/// // Create variation generator for Earth-like planet
/// let mut var = PlanetaryVariation::new_from_planet(1.0, 1.0);
///
/// // Vary a composition fraction by ±10%
/// let iron_fraction = var.vary_composition(0.32, 0.10);
/// assert!(iron_fraction >= 0.0 && iron_fraction <= 1.0);
///
/// // Same planet always produces same sequence
/// let mut var2 = PlanetaryVariation::new_from_planet(1.0, 1.0);
/// let iron_fraction2 = var2.vary_composition(0.32, 0.10);
/// assert_eq!(iron_fraction, iron_fraction2);
/// ```
pub struct PlanetaryVariation {
    rng: ChaCha8Rng,
}

impl PlanetaryVariation {
    /// Creates a new variation generator seeded from planet properties
    ///
    /// The seed is deterministic based on mass and orbital distance,
    /// ensuring reproducible results for the same planet.
    ///
    /// # Arguments
    /// * `mass_earth` - Planet mass in Earth masses
    /// * `semi_major_axis_au` - Orbital semi-major axis in AU
    pub fn new_from_planet(mass_earth: f64, semi_major_axis_au: f64) -> Self {
        // Hash the bit representations of both floats for proper mixing
        let mut hasher = std::collections::hash_map::DefaultHasher::new();
        mass_earth.to_bits().hash(&mut hasher);
        semi_major_axis_au.to_bits().hash(&mut hasher);
        let seed = hasher.finish();

        Self::from_seed(seed)
    }

    /// Creates a variation generator from an explicit seed
    ///
    /// Useful when you want full control over the random sequence.
    pub fn from_seed(seed: u64) -> Self {
        use rand::SeedableRng;
        Self {
            rng: ChaCha8Rng::seed_from_u64(seed),
        }
    }

    /// Varies a composition fraction within realistic bounds
    ///
    /// # Arguments
    /// * `base_value` - The baseline composition fraction (0.0 to 1.0)
    /// * `variation_factor` - How much to vary (e.g., 0.1 = ±10%)
    ///
    /// # Returns
    /// A value within realistic bounds, clamped to [0.0, 1.0]
    pub fn vary_composition(&mut self, base_value: f64, variation_factor: f64) -> f64 {
        let variation = self.rng.random_range(-variation_factor..=variation_factor);
        (base_value + variation).clamp(0.0, 1.0)
    }

    /// Varies a physical property (mass, radius, etc.) within a percentage range
    ///
    /// # Arguments
    /// * `base_value` - The baseline value
    /// * `percent_variation` - Percentage variation (e.g., 5.0 = ±5%)
    ///
    /// # Returns
    /// A varied value within the specified percentage range
    pub fn vary_property(&mut self, base_value: f64, percent_variation: f64) -> f64 {
        let factor = percent_variation / 100.0;
        let variation = self.rng.random_range(-factor..=factor);
        base_value * (1.0 + variation)
    }

    /// Varies albedo within realistic bounds for planet type
    ///
    /// # Arguments
    /// * `base_albedo` - The baseline albedo value
    /// * `min_albedo` - Minimum realistic albedo
    /// * `max_albedo` - Maximum realistic albedo
    ///
    /// # Returns
    /// A varied albedo clamped to realistic bounds
    pub fn vary_albedo(&mut self, base_albedo: f64, min_albedo: f64, max_albedo: f64) -> f64 {
        let variation = self.rng.random_range(-0.05..=0.05);
        (base_albedo + variation).clamp(min_albedo, max_albedo)
    }

    /// Adds subtle color variation to planet colors
    ///
    /// # Arguments
    /// * `base_color` - RGB tuple (0.0 to 1.0)
    /// * `variation` - Amount of variation (typically 0.02-0.05)
    ///
    /// # Returns
    /// Varied RGB tuple with values clamped to [0.0, 1.0]
    pub fn vary_color(&mut self, base_color: (f64, f64, f64), variation: f64) -> (f64, f64, f64) {
        let (r, g, b) = base_color;
        (
            (r + self.rng.random_range(-variation..=variation)).clamp(0.0, 1.0),
            (g + self.rng.random_range(-variation..=variation)).clamp(0.0, 1.0),
            (b + self.rng.random_range(-variation..=variation)).clamp(0.0, 1.0),
        )
    }

    /// Generates a random offset for orbital eccentricity effects
    ///
    /// Useful for varying temperature or other orbital-dependent properties.
    pub fn orbital_variation(&mut self) -> f64 {
        self.rng.random_range(-0.05..=0.05)
    }

    /// Generates a small random integer in a range (inclusive)
    pub fn random_int(&mut self, min: i32, max: i32) -> i32 {
        self.rng.random_range(min..=max)
    }

    /// Generates a random boolean with given probability
    ///
    /// # Arguments
    /// * `probability` - Probability of true (0.0 to 1.0)
    pub fn random_bool(&mut self, probability: f64) -> bool {
        self.rng.random_range(0.0..1.0) < probability
    }

    /// Generates a random value from 0.0 to 1.0
    pub fn random_unit(&mut self) -> f64 {
        self.rng.random_range(0.0..1.0)
    }

    /// Generates a random value in a range
    pub fn random_range(&mut self, min: f64, max: f64) -> f64 {
        self.rng.random_range(min..=max)
    }
}
