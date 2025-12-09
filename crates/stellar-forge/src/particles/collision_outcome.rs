//! Collision outcome physics for dust particles.
//!
//! Determines what happens when two particles collide based on their
//! relative velocity, sizes, and material properties.
//!
//! # Physics
//!
//! Collision outcomes depend critically on impact velocity:
//!
//! - **v < v_bounce (~0.01-0.1 m/s)**: Perfect sticking
//! - **v_bounce < v < v_frag (~1-10 m/s)**: Probabilistic sticking/bouncing
//! - **v > v_frag**: Fragmentation and erosion
//!
//! These thresholds depend on particle composition, size, and porosity.
//!
//! # References
//!
//! - Güttler et al. (2010) - Laboratory collision experiments
//! - Windmark et al. (2012) - Bouncing and fragmentation barriers
//! - Birnstiel et al. (2012) - Fragmentation-limited dust evolution

use rand::Rng;
use rand_chacha::ChaChaRng;
use units::Velocity;

/// Outcome of a collision between two particles.
///
/// The outcome depends on the impact velocity and material properties.
#[derive(Debug, Clone, PartialEq)]
pub enum CollisionOutcome {
    /// Particles stick together perfectly.
    ///
    /// All mass is combined into a single larger particle.
    PerfectSticking,

    /// Particles stick together but with reduced efficiency.
    ///
    /// Only a fraction of the mass actually sticks. The rest bounces or
    /// is redistributed.
    PartialSticking {
        /// Fraction of mass that sticks (0 to 1)
        efficiency: f64,
    },

    /// Particles bounce without mass transfer.
    ///
    /// No growth occurs. This creates a "bouncing barrier" where particles
    /// cannot grow beyond a certain size.
    Bouncing,

    /// Mass is transferred from one particle to another.
    ///
    /// Occurs at intermediate velocities where smaller particles can
    /// erode the surface of larger ones.
    MassTransfer {
        /// True if mass goes from larger to smaller particle
        from_larger: bool,
        /// Fraction of particle mass transferred
        fraction: f64,
    },

    /// Particles fragment into smaller pieces.
    ///
    /// Occurs at high velocities. The resulting fragments follow a
    /// power-law size distribution.
    Fragmentation {
        /// Power-law exponent for fragment distribution
        ///
        /// Typical values: -1.8 to -2.5 (steeper than collisional cascade)
        power_law_exponent: f64,
    },

    /// Surface erosion of target particle.
    ///
    /// High-velocity small projectile chips off mass from larger target.
    Erosion {
        /// Fraction of target mass lost
        mass_loss_fraction: f64,
    },
}

/// Material properties that determine collision physics.
#[derive(Debug, Clone)]
pub struct MaterialProperties {
    /// Fragmentation velocity threshold (m/s)
    ///
    /// Above this velocity, collisions result in fragmentation.
    ///
    /// Typical values:
    /// - Silicates (basalt): ~1-3 m/s
    /// - Ices (water): ~10-20 m/s
    /// - Fluffy aggregates: ~0.1-1 m/s
    pub fragmentation_velocity: f64,

    /// Bouncing velocity threshold (m/s)
    ///
    /// Below this, particles stick. Above, they may bounce.
    ///
    /// Typical values:
    /// - Compact silicates: ~0.01-0.1 m/s
    /// - Icy grains: ~0.1-1 m/s
    /// - Fractal aggregates: ~0.001-0.01 m/s
    pub bouncing_velocity: f64,

    /// Erosion threshold (m/s)
    ///
    /// High-velocity small particles can erode larger targets.
    pub erosion_velocity: f64,
}

impl MaterialProperties {
    /// Compact silicate particles (basalt, olivine).
    ///
    /// Based on Güttler et al. (2010) experiments.
    pub fn compact_silicates() -> Self {
        Self {
            fragmentation_velocity: 2.0, // ~2 m/s
            bouncing_velocity: 0.05,     // ~5 cm/s
            erosion_velocity: 5.0,       // ~5 m/s
        }
    }

    /// Water ice particles.
    ///
    /// Ice is stronger than silicates at disk temperatures.
    pub fn water_ice() -> Self {
        Self {
            fragmentation_velocity: 15.0, // ~15 m/s
            bouncing_velocity: 0.5,       // ~50 cm/s
            erosion_velocity: 20.0,       // ~20 m/s
        }
    }

    /// Fluffy, porous aggregates.
    ///
    /// Low-density fractal structures stick easily but fragment easily.
    pub fn fluffy_aggregates() -> Self {
        Self {
            fragmentation_velocity: 0.5, // ~0.5 m/s
            bouncing_velocity: 0.01,     // ~1 cm/s
            erosion_velocity: 1.0,       // ~1 m/s
        }
    }
}

impl CollisionOutcome {
    /// Determine the collision outcome based on impact velocity.
    ///
    /// Uses probabilistic transitions between regimes to capture the
    /// stochastic nature of particle collisions.
    ///
    /// # Arguments
    /// * `impact_velocity` - Relative velocity of collision
    /// * `size_ratio` - Ratio of larger to smaller particle size
    /// * `material` - Material properties
    /// * `rng` - Random number generator for probabilistic outcomes
    ///
    /// # Example
    /// ```
    /// use stellar_forge::particles::{CollisionOutcome, MaterialProperties};
    /// use units::Velocity;
    /// use rand::SeedableRng;
    /// use rand_chacha::ChaChaRng;
    ///
    /// let mut rng = ChaChaRng::seed_from_u64(42);
    /// let material = MaterialProperties::compact_silicates();
    ///
    /// // Low velocity: sticking
    /// let outcome = CollisionOutcome::sample(
    ///     Velocity::from_meters_per_sec(0.01),
    ///     1.0,
    ///     &material,
    ///     &mut rng,
    /// );
    /// assert_eq!(outcome, CollisionOutcome::PerfectSticking);
    ///
    /// // High velocity: fragmentation
    /// let outcome = CollisionOutcome::sample(
    ///     Velocity::from_meters_per_sec(5.0),
    ///     1.0,
    ///     &material,
    ///     &mut rng,
    /// );
    /// assert!(matches!(outcome, CollisionOutcome::Fragmentation { .. }));
    /// ```
    pub fn sample(
        impact_velocity: Velocity,
        size_ratio: f64,
        material: &MaterialProperties,
        rng: &mut ChaChaRng,
    ) -> Self {
        let v = impact_velocity.to_meters_per_sec();
        let v_bounce = material.bouncing_velocity;
        let v_frag = material.fragmentation_velocity;
        let v_erosion = material.erosion_velocity;

        // Erosion regime: small fast projectile hitting large target
        if size_ratio > 10.0 && v > v_erosion {
            let mass_loss = Self::erosion_mass_loss(v, v_erosion, rng);
            return Self::Erosion {
                mass_loss_fraction: mass_loss,
            };
        }

        // Fragmentation regime: high velocity
        if v > v_frag {
            // Power-law exponent from impact experiments
            // Typical range: -1.8 to -2.5
            let exponent = -2.0 - rng.random::<f64>() * 0.5;
            return Self::Fragmentation {
                power_law_exponent: exponent,
            };
        }

        // Intermediate regime: probabilistic sticking/bouncing
        if v > v_bounce {
            // Sticking probability decreases linearly from v_bounce to v_frag
            let stick_prob = (v_frag - v) / (v_frag - v_bounce);

            if rng.random::<f64>() < stick_prob {
                // Partial sticking with decreasing efficiency
                let efficiency = 0.5 + 0.5 * rng.random::<f64>();
                Self::PartialSticking { efficiency }
            } else {
                // Bouncing: no mass transfer
                Self::Bouncing
            }
        } else {
            // Low velocity regime: perfect sticking
            Self::PerfectSticking
        }
    }

    /// Calculate erosion mass loss fraction.
    ///
    /// Based on impact crater scaling laws.
    fn erosion_mass_loss(v: f64, v_threshold: f64, rng: &mut ChaChaRng) -> f64 {
        // Erosion efficiency increases with velocity
        // Typical values: 0.01 - 0.3 of target mass
        let base_efficiency = 0.05;
        let velocity_factor = (v / v_threshold).powi(2);
        let randomization = 0.5 + rng.random::<f64>();

        (base_efficiency * velocity_factor * randomization).min(0.3)
    }

    /// Returns true if this outcome results in net mass growth.
    pub fn is_growth(&self) -> bool {
        matches!(self, Self::PerfectSticking | Self::PartialSticking { .. })
    }

    /// Returns true if this outcome results in fragmentation.
    pub fn is_fragmentation(&self) -> bool {
        matches!(self, Self::Fragmentation { .. } | Self::Erosion { .. })
    }

    /// Returns true if this outcome is a bouncing collision.
    pub fn is_bouncing(&self) -> bool {
        matches!(self, Self::Bouncing)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;

    #[test]
    fn low_velocity_sticks() {
        let mut rng = ChaChaRng::seed_from_u64(42);
        let material = MaterialProperties::compact_silicates();

        // Well below bouncing threshold
        let outcome = CollisionOutcome::sample(
            Velocity::from_meters_per_sec(0.001),
            1.0,
            &material,
            &mut rng,
        );

        assert_eq!(outcome, CollisionOutcome::PerfectSticking);
    }

    #[test]
    fn high_velocity_fragments() {
        let mut rng = ChaChaRng::seed_from_u64(42);
        let material = MaterialProperties::compact_silicates();

        // Well above fragmentation threshold
        let outcome = CollisionOutcome::sample(
            Velocity::from_meters_per_sec(10.0),
            1.0,
            &material,
            &mut rng,
        );

        assert!(outcome.is_fragmentation());
    }

    #[test]
    fn intermediate_velocity_probabilistic() {
        let mut rng = ChaChaRng::seed_from_u64(42);
        let material = MaterialProperties::compact_silicates();

        // Between bouncing and fragmentation thresholds
        let v_test = (material.bouncing_velocity + material.fragmentation_velocity) / 2.0;

        // Sample multiple times to check we get varied outcomes
        let mut outcomes = Vec::new();
        for _ in 0..100 {
            let outcome = CollisionOutcome::sample(
                Velocity::from_meters_per_sec(v_test),
                1.0,
                &material,
                &mut rng,
            );
            outcomes.push(outcome);
        }

        // Should see both sticking and bouncing
        let has_sticking = outcomes.iter().any(|o| o.is_growth());
        let has_bouncing = outcomes.iter().any(|o| o.is_bouncing());

        assert!(has_sticking, "Should see some sticking events");
        assert!(has_bouncing, "Should see some bouncing events");
    }

    #[test]
    fn erosion_requires_large_size_ratio() {
        let mut rng = ChaChaRng::seed_from_u64(42);
        let material = MaterialProperties::compact_silicates();

        // High velocity but equal sizes: fragmentation
        let outcome = CollisionOutcome::sample(
            Velocity::from_meters_per_sec(10.0),
            1.0,
            &material,
            &mut rng,
        );
        assert!(matches!(outcome, CollisionOutcome::Fragmentation { .. }));

        // High velocity and large size ratio: erosion
        let outcome = CollisionOutcome::sample(
            Velocity::from_meters_per_sec(10.0),
            20.0,
            &material,
            &mut rng,
        );
        assert!(matches!(outcome, CollisionOutcome::Erosion { .. }));
    }

    #[test]
    fn ice_has_higher_thresholds() {
        let silicates = MaterialProperties::compact_silicates();
        let ice = MaterialProperties::water_ice();

        assert!(ice.fragmentation_velocity > silicates.fragmentation_velocity);
        assert!(ice.bouncing_velocity > silicates.bouncing_velocity);
    }

    #[test]
    fn fluffy_aggregates_fragile() {
        let fluffy = MaterialProperties::fluffy_aggregates();
        let compact = MaterialProperties::compact_silicates();

        assert!(fluffy.fragmentation_velocity < compact.fragmentation_velocity);
        assert!(fluffy.bouncing_velocity < compact.bouncing_velocity);
    }
}
