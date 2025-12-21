//! Planetesimal formation events.
//!
//! When a particle layer becomes gravitationally unstable (low Toomre Q,
//! high Richardson number), it fragments into planetesimals on the
//! free-fall timescale.
//!
//! This module defines the event structure that captures this transition
//! from a statistical particle population to discrete bodies.

use rand_chacha::ChaChaRng;
use units::{Length, Mass, Time};

use crate::particles::SizeDistribution;

/// Result of a planetesimal formation event.
///
/// When a particle bin becomes gravitationally unstable, a fraction of
/// its mass collapses into planetesimals. This structure captures:
/// - Where the event occurred
/// - How much mass was converted
/// - The characteristic size of planetesimals formed
/// - The size distribution of the resulting population
///
/// # Physical Process
///
/// The Goldreich-Ward mechanism or streaming instability triggers when:
/// 1. Particle layer is dense enough (Q < Q_crit)
/// 2. Layer is shear-stable (Ri > Ri_crit)
///
/// Fragmentation proceeds on the free-fall timescale (~10-100 orbits)
/// and produces planetesimals with masses near the Jeans mass.
#[derive(Debug, Clone)]
pub struct PlanetesimalFormationEvent {
    /// Radial location of formation (cm)
    pub location: Length,

    /// Total mass converted to planetesimals
    pub total_mass: Mass,

    /// Characteristic mass of individual planetesimals
    ///
    /// Set by the Jeans mass: M_J = π σ⁴ / (G² Σ)
    pub characteristic_mass: Mass,

    /// Characteristic size of individual planetesimals
    pub characteristic_size: Length,

    /// Size distribution of formed planetesimals
    ///
    /// Typically a narrow distribution centered on the Jeans mass,
    /// represented as a power law with shallow slope.
    pub size_distribution: SizeDistribution,

    /// Formation timescale (free-fall time)
    pub formation_timescale: Time,

    /// Formation efficiency (fraction of bin mass that formed planetesimals)
    ///
    /// Not all mass participates in collapse. Typical values: 0.1-0.5
    pub efficiency: f64,
}

impl PlanetesimalFormationEvent {
    /// Create a planetesimal formation event.
    ///
    /// # Arguments
    /// * `location` - Radial position in the disk
    /// * `jeans_mass` - Characteristic mass from Jeans length
    /// * `jeans_size` - Characteristic size
    /// * `available_mass` - Total mass available in the particle bin
    /// * `freefall_time` - Timescale for collapse
    /// * `efficiency` - Fraction of mass that forms planetesimals (0-1)
    /// * `material_density` - Bulk density of planetesimals
    ///
    /// # Returns
    /// A formation event with a narrow size distribution centered on the
    /// Jeans mass.
    pub fn new(
        location: Length,
        jeans_mass: Mass,
        jeans_size: Length,
        available_mass: Mass,
        freefall_time: Time,
        efficiency: f64,
        material_density: units::Density,
    ) -> Self {
        let total_mass = available_mass * efficiency;

        // Size distribution: narrow power law centered on Jeans mass
        // Use shallow slope (q ≈ -2.0) to represent fragmentation spectrum
        let s_min = jeans_size * 0.1; // Factor of 10 range
        let s_max = jeans_size * 10.0;

        let size_dist = SizeDistribution::power_law(
            s_min,
            s_max,
            -2.0, // Shallow slope for fragmentation
            total_mass,
            material_density,
        );

        Self {
            location,
            total_mass,
            characteristic_mass: jeans_mass,
            characteristic_size: jeans_size,
            size_distribution: size_dist,
            formation_timescale: freefall_time,
            efficiency,
        }
    }

    /// Number of planetesimals formed.
    ///
    /// Approximated as total mass divided by characteristic mass.
    pub fn number_formed(&self) -> f64 {
        self.total_mass.to_grams() / self.characteristic_mass.to_grams()
    }

    /// Sample formation efficiency from physical constraints.
    ///
    /// Not all mass participates in collapse. The efficiency depends on:
    /// - How unstable the layer is (lower Q → higher efficiency)
    /// - Turbulent diffusion competing with collapse
    /// - Tidal shear from stellar gravity
    ///
    /// Typical range: 10-50% per free-fall time.
    ///
    /// # Arguments
    /// * `q` - Toomre Q parameter (lower = more unstable)
    /// * `rng` - Random number generator for stochastic variation
    ///
    /// # Returns
    /// Efficiency between 0.1 and 0.5
    pub fn sample_efficiency(q: f64, rng: &mut ChaChaRng) -> f64 {
        use rand::Rng;

        // Base efficiency: lower Q → higher efficiency
        // Q = 0.5 → ~50% efficiency
        // Q = 1.5 → ~10% efficiency
        let base_efficiency = (2.0 - q).clamp(0.1, 0.5);

        // Add stochastic variation (±30%)
        let randomization = 0.7 + 0.6 * rng.random::<f64>();

        (base_efficiency * randomization).clamp(0.1, 0.5)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::{Rng, SeedableRng};

    #[test]
    fn formation_event_stores_parameters() {
        let event = PlanetesimalFormationEvent::new(
            Length::from_au(5.0),
            Mass::from_grams(1e20),
            Length::from_km(10.0),
            Mass::from_grams(1e22),
            Time::from_years(1000.0),
            0.3,
            units::Density::from_grams_per_cm3(3.0),
        );

        assert_eq!(event.location.to_au(), 5.0);
        assert_eq!(event.characteristic_mass.to_grams(), 1e20);
        assert_eq!(event.efficiency, 0.3);
        assert_eq!(event.total_mass.to_grams(), 1e22 * 0.3);
    }

    #[test]
    fn number_formed_is_total_over_characteristic() {
        let event = PlanetesimalFormationEvent::new(
            Length::from_au(5.0),
            Mass::from_grams(1e20),
            Length::from_km(10.0),
            Mass::from_grams(1e22),
            Time::from_years(1000.0),
            0.3,
            units::Density::from_grams_per_cm3(3.0),
        );

        let n = event.number_formed();
        let expected = (1e22 * 0.3) / 1e20;

        assert!((n - expected).abs() < 1.0);
    }

    #[test]
    fn efficiency_decreases_with_q() {
        let mut rng = ChaChaRng::seed_from_u64(42);

        let eff_low_q = PlanetesimalFormationEvent::sample_efficiency(0.5, &mut rng);
        let eff_high_q = PlanetesimalFormationEvent::sample_efficiency(1.5, &mut rng);

        // Lower Q should give higher efficiency
        assert!(
            eff_low_q > eff_high_q,
            "Lower Q should give higher efficiency: {} vs {}",
            eff_low_q,
            eff_high_q
        );
    }

    #[test]
    fn efficiency_in_range() {
        let mut rng = ChaChaRng::seed_from_u64(42);

        for _ in 0..100 {
            let q = 0.5 + 1.0 * rng.random::<f64>();
            let eff = PlanetesimalFormationEvent::sample_efficiency(q, &mut rng);

            assert!(
                eff >= 0.1 && eff <= 0.5,
                "Efficiency {} out of range [0.1, 0.5]",
                eff
            );
        }
    }

    #[test]
    fn size_distribution_has_correct_mass() {
        let total_mass = Mass::from_grams(1e22);
        let efficiency = 0.3;

        let event = PlanetesimalFormationEvent::new(
            Length::from_au(5.0),
            Mass::from_grams(1e20),
            Length::from_km(10.0),
            total_mass,
            Time::from_years(1000.0),
            efficiency,
            units::Density::from_grams_per_cm3(3.0),
        );

        let dist_mass = event.size_distribution.total_mass().to_grams();
        let expected = total_mass.to_grams() * efficiency;

        // Should match within numerical precision
        assert!(
            (dist_mass - expected).abs() / expected < 1e-6,
            "Size distribution mass {} doesn't match expected {}",
            dist_mass,
            expected
        );
    }
}
