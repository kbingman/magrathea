//! Discrete body representation for planetesimals, embryos, and protoplanets.
//!
//! Bodies are tracked individually once they become massive enough or when
//! gravitational interactions become important. They use Cartesian coordinates
//! for integration (compatible with nbody crate) with orbital elements stored
//! for physics calculations.

use nalgebra::{Point2, Vector2};
use planetary::composition::Composition;
use serde::{Deserialize, Serialize};
use units::{Length, Mass, MassRate, Time, Velocity};

use super::envelope;
use super::orbital_elements::cartesian_to_orbital_elements;
use crate::disk::DiskModel;

/// Gravitational constant in AU³ M☉⁻¹ year⁻²
const G: f64 = 39.478417;

// Re-export BodyId from nbody for convenience
pub use nbody::body::BodyId;

/// Default BodyId for deserialization
fn default_body_id() -> BodyId {
    BodyId(0)
}

/// Gas envelope state for protoplanets.
///
/// Tracks the evolution of a protoplanet's gas envelope through different
/// phases of gas accretion.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum EnvelopeState {
    /// No envelope captured yet
    None,

    /// Hydrostatic envelope in equilibrium
    ///
    /// Envelope mass grows slowly via Kelvin-Helmholtz contraction.
    /// Core mass < critical mass.
    Hydrostatic { envelope_mass: Mass },

    /// Runaway gas accretion phase
    ///
    /// Core mass exceeded critical mass. Envelope mass grows rapidly
    /// until gas supply is exhausted or gap is opened.
    Runaway {
        envelope_mass: Mass,
        accretion_rate: MassRate,
    },

    /// Final state (gas disk dispersed)
    Final { envelope_mass: Mass },
}

/// A discrete body in a protoplanetary disk simulation.
///
/// Represents planetesimals, embryos, and protoplanets as individual objects.
/// Uses Cartesian coordinates for integration, with orbital elements stored
/// and updated periodically for physics calculations and output.
///
/// # Integration with nbody
///
/// Implements `Massive` trait to work with Barnes-Hut tree gravity and
/// symplectic integrators from the nbody crate.
///
/// # Examples
///
/// ```
/// use protodisk::bodies::DiscreteBody;
/// use units::Mass;
/// use nalgebra::{Point2, Vector2};
/// use planetary::composition::Composition;
///
/// let body = DiscreteBody::new(
///     Mass::from_earth_masses(0.1),
///     Mass::zero(),
///     Point2::new(5.0, 0.0),
///     Vector2::new(0.0, 2.8),
///     Mass::from_solar_masses(1.0),
///     Composition::earth_like(),
/// );
///
/// let r_hill = body.hill_radius(Mass::from_solar_masses(1.0));
/// ```
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DiscreteBody {
    /// Unique identifier (skipped in serialization, assigned by system)
    #[serde(skip, default = "default_body_id")]
    pub id: BodyId,

    // ===== Cartesian State (source of truth for integration) =====
    /// Position in AU (heliocentric)
    pub position: Point2<f64>,

    /// Velocity in AU/year
    pub velocity: Vector2<f64>,

    // ===== Stored Orbital Elements (updated periodically) =====
    /// Semi-major axis (AU)
    pub semi_major_axis: Length,

    /// Eccentricity (dimensionless)
    pub eccentricity: f64,

    /// Inclination (radians, 0 for 2D simulations)
    pub inclination: f64,

    // ===== Physical Properties =====
    /// Core mass (rock + metal)
    pub core_mass: Mass,

    /// Envelope mass (H/He gas)
    pub envelope_mass: Mass,

    /// Physical radius for collisions
    pub physical_radius: Length,

    /// Bulk composition (mass fractions)
    pub composition: Composition,

    // ===== Planet Formation State =====
    /// Gas envelope state
    pub envelope_state: EnvelopeState,
}

impl DiscreteBody {
    /// Create a new discrete body.
    ///
    /// Computes initial orbital elements from Cartesian state.
    ///
    /// # Arguments
    /// * `core_mass` - Initial core mass
    /// * `envelope_mass` - Initial envelope mass (typically zero)
    /// * `position` - Initial position (AU)
    /// * `velocity` - Initial velocity (AU/year)
    /// * `stellar_mass` - Mass of central star (for orbital elements)
    /// * `composition` - Bulk composition
    ///
    /// # Returns
    /// New DiscreteBody with computed orbital elements
    pub fn new(
        core_mass: Mass,
        envelope_mass: Mass,
        position: Point2<f64>,
        velocity: Vector2<f64>,
        stellar_mass: Mass,
        composition: Composition,
    ) -> Self {
        let elements = cartesian_to_orbital_elements(position, velocity, stellar_mass);

        // Estimate physical radius from mass-radius relation
        let total_mass = core_mass + envelope_mass;
        let physical_radius = Self::estimate_radius(total_mass, &composition);

        Self {
            id: BodyId(0), // Will be assigned by system
            position,
            velocity,
            semi_major_axis: elements.semi_major_axis,
            eccentricity: elements.eccentricity,
            inclination: elements.inclination,
            core_mass,
            envelope_mass,
            physical_radius,
            composition,
            envelope_state: if envelope_mass > Mass::zero() {
                EnvelopeState::Hydrostatic { envelope_mass }
            } else {
                EnvelopeState::None
            },
        }
    }

    /// Total mass (core + envelope).
    pub fn total_mass(&self) -> Mass {
        self.core_mass + self.envelope_mass
    }

    /// Update stored orbital elements from current Cartesian state.
    ///
    /// Should be called after nbody integration or whenever elements are needed.
    ///
    /// # Arguments
    /// * `stellar_mass` - Mass of central star
    pub fn update_orbital_elements(&mut self, stellar_mass: Mass) {
        let elements = cartesian_to_orbital_elements(self.position, self.velocity, stellar_mass);
        self.semi_major_axis = elements.semi_major_axis;
        self.eccentricity = elements.eccentricity;
        self.inclination = elements.inclination;
    }

    /// Update from an nbody::Body after integration.
    ///
    /// Syncs position and velocity from the nbody integrator and updates
    /// stored orbital elements.
    ///
    /// # Arguments
    /// * `body` - Body from nbody integration
    /// * `stellar_mass` - Mass of central star
    pub fn update_from_nbody(&mut self, body: &nbody::body::Body, stellar_mass: Mass) {
        self.position = body.position;
        self.velocity = body.velocity;
        self.update_orbital_elements(stellar_mass);
    }

    /// Convert to nbody::Body for integration.
    ///
    /// # Returns
    /// Body compatible with nbody crate integrators
    pub fn to_nbody_body(&self) -> nbody::body::Body {
        nbody::body::Body {
            id: self.id,
            mass: self.total_mass().to_solar_masses(),
            radius: self.physical_radius.to_au(),
            position: self.position,
            velocity: self.velocity,
        }
    }

    // ===== Physics Methods =====

    /// Hill radius (gravitational sphere of influence).
    ///
    /// The Hill radius is the approximate radius within which the body's
    /// gravity dominates over the star's tidal forces.
    ///
    /// R_H = a × (M / 3M_star)^(1/3)
    ///
    /// # Arguments
    /// * `stellar_mass` - Mass of central star
    ///
    /// # Returns
    /// Hill radius in AU
    pub fn hill_radius(&self, stellar_mass: Mass) -> Length {
        let mass_ratio = self.total_mass().to_solar_masses() / stellar_mass.to_solar_masses();
        let factor = (mass_ratio / 3.0).powf(1.0 / 3.0);
        Length::from_au(self.semi_major_axis.to_au() * factor)
    }

    /// Escape velocity from the body's surface.
    ///
    /// v_esc = sqrt(2GM/R)
    ///
    /// # Returns
    /// Escape velocity in AU/year
    pub fn escape_velocity(&self) -> Velocity {
        let mu = G * self.total_mass().to_solar_masses();
        let r = self.physical_radius.to_au();
        Velocity::from_au_per_year((2.0 * mu / r).sqrt())
    }

    /// Gravitational focusing factor for accretion.
    ///
    /// Accounts for gravitational attraction enhancing collision cross-section.
    /// F_g = 1 + (v_esc / σ)²
    ///
    /// Where σ is the velocity dispersion of the accreted material.
    ///
    /// # Arguments
    /// * `velocity_dispersion` - Velocity dispersion of particle population
    ///
    /// # Returns
    /// Focusing factor (≥ 1)
    pub fn gravitational_focusing_factor(&self, velocity_dispersion: Velocity) -> f64 {
        let v_esc = self.escape_velocity().to_au_per_year();
        let sigma = velocity_dispersion.to_au_per_year();
        1.0 + (v_esc / sigma).powi(2)
    }

    /// Orbital frequency (mean motion).
    ///
    /// n = sqrt(GM/a³)
    ///
    /// # Arguments
    /// * `stellar_mass` - Mass of central star
    ///
    /// # Returns
    /// Angular frequency in radians/year
    pub fn orbital_frequency(&self, stellar_mass: Mass) -> f64 {
        let mu = G * stellar_mass.to_solar_masses();
        let a = self.semi_major_axis.to_au();
        (mu / a.powi(3)).sqrt()
    }

    /// Orbital period.
    ///
    /// T = 2π / n
    ///
    /// # Arguments
    /// * `stellar_mass` - Mass of central star
    ///
    /// # Returns
    /// Orbital period in years
    pub fn orbital_period(&self, stellar_mass: Mass) -> f64 {
        std::f64::consts::TAU / self.orbital_frequency(stellar_mass)
    }

    // ===== Growth Physics Methods =====

    /// Accretion rate from a particle population.
    ///
    /// Calculates the rate at which this body accretes solids from a
    /// particle bin, using the particle-in-a-box approximation with
    /// gravitational focusing.
    ///
    /// dM/dt = √3/2 × Σ_p × Ω × R² × F_g
    ///
    /// where:
    /// - Σ_p is the particle surface density
    /// - Ω is the orbital frequency
    /// - R is the physical radius
    /// - F_g is the gravitational focusing factor
    ///
    /// # Arguments
    /// * `particle_surface_density` - Surface density of solids
    /// * `particle_velocity_dispersion` - Random velocity of particles
    /// * `stellar_mass` - Mass of central star
    ///
    /// # Returns
    /// Accretion rate in solar masses per year
    ///
    /// # References
    /// - Armitage (2010) - Astrophysics of Planet Formation, §6.2
    /// - Lissauer (1993) - Particle accretion in the protoplanetary nebula
    pub fn accretion_rate_from_particles(
        &self,
        particle_surface_density: units::SurfaceDensity,
        particle_velocity_dispersion: units::Velocity,
        stellar_mass: Mass,
    ) -> units::MassRate {
        let sigma_p = particle_surface_density.to_grams_per_cm2();
        let omega = self.orbital_frequency(stellar_mass); // radians/year
        let r_cm = self.physical_radius.to_cm();
        let f_g = self.gravitational_focusing_factor(particle_velocity_dispersion);

        // dM/dt = √3/2 × Σ_p × Ω × R² × F_g
        // In CGS: g/year = (g/cm²) × (1/year) × cm² × (dimensionless)
        let dm_dt_cgs = (3.0_f64.sqrt() / 2.0) * sigma_p * omega * r_cm * r_cm * f_g;

        units::MassRate::from_grams_per_year(dm_dt_cgs)
    }

    /// Width of the feeding zone.
    ///
    /// The radial region from which a body can efficiently accrete material.
    /// Typically approximated as ~2.5 Hill radii on either side of the orbit.
    ///
    /// Δa ≈ 2.5 × R_H
    ///
    /// # Arguments
    /// * `stellar_mass` - Mass of central star
    ///
    /// # Returns
    /// Feeding zone width in AU
    ///
    /// # References
    /// - Lissauer (1987) - Timescales for planetary accretion
    pub fn feeding_zone_width(&self, stellar_mass: Mass) -> Length {
        let r_hill = self.hill_radius(stellar_mass);
        r_hill * 2.5
    }

    /// Isolation mass for this orbit.
    ///
    /// The maximum mass a body can reach by clearing its feeding zone
    /// of all available solids. Represents the transition to oligarchic growth.
    ///
    /// M_iso = 2π × a × Δa × Σ
    ///
    /// where:
    /// - a is the semi-major axis
    /// - Δa is the feeding zone width
    /// - Σ is the surface density of solids
    ///
    /// # Arguments
    /// * `surface_density` - Surface density of available solids
    /// * `stellar_mass` - Mass of central star
    ///
    /// # Returns
    /// Isolation mass in solar masses
    ///
    /// # References
    /// - Lissauer (1993) - Planet formation
    /// - Kokubo & Ida (1998) - Oligarchic growth of protoplanets
    pub fn isolation_mass(
        &self,
        surface_density: units::SurfaceDensity,
        stellar_mass: Mass,
    ) -> Mass {
        let a = self.semi_major_axis.to_au();
        let delta_a = self.feeding_zone_width(stellar_mass).to_au();
        let sigma = surface_density.to_solar_masses_per_au2();

        // M_iso = 2π × a × Δa × Σ
        let m_iso = 2.0 * std::f64::consts::PI * a * delta_a * sigma;

        Mass::from_solar_masses(m_iso)
    }

    // ===== Gas Envelope Evolution =====

    /// Evolve gas envelope for one timestep.
    ///
    /// Handles envelope state transitions:
    /// - None → Hydrostatic (when core mass sufficient and gas available)
    /// - Hydrostatic → Runaway (when core exceeds critical mass)
    /// - Runaway → Final (when gas disk disperses or gap opens)
    ///
    /// # Arguments
    /// * `disk` - Gas disk model
    /// * `core_accretion_rate` - Rate of planetesimal accretion (for luminosity)
    /// * `dt` - Timestep in years
    /// * `opacity` - Grain opacity in cm²/g (typically 0.01-10)
    ///
    /// # Physics
    ///
    /// Before critical mass: envelope grows via Kelvin-Helmholtz contraction
    /// After critical mass: runaway accretion limited by disk supply or gap
    pub fn evolve_envelope<D: DiskModel>(
        &mut self,
        disk: &D,
        core_accretion_rate: MassRate,
        dt: Time,
        opacity: f64,
    ) {
        let location = self.semi_major_axis;
        let gas_sigma = disk.surface_density(location);

        // For Runaway state, compute gap opening and local gas mass before mutable borrow
        let (opens_gap, local_gas_mass) =
            if matches!(self.envelope_state, EnvelopeState::Runaway { .. }) {
                let stellar_mass = disk.stellar_mass();
                (
                    self.opens_gap(disk, stellar_mass),
                    self.local_gas_mass(disk),
                )
            } else {
                (false, Mass::zero())
            };

        match &mut self.envelope_state {
            EnvelopeState::None => {
                // Check if conditions allow envelope capture
                if envelope::can_capture_envelope(self.core_mass, gas_sigma) {
                    // Start with small envelope mass
                    let initial_envelope = Mass::from_earth_masses(0.01);
                    self.envelope_mass = initial_envelope;
                    self.envelope_state = EnvelopeState::Hydrostatic {
                        envelope_mass: initial_envelope,
                    };
                }
            }

            EnvelopeState::Hydrostatic { envelope_mass } => {
                // Check if core exceeded critical mass
                let m_crit = envelope::critical_core_mass(core_accretion_rate, opacity);

                if self.core_mass > m_crit {
                    // Transition to runaway accretion!
                    let initial_rate = envelope::supply_limited_accretion_rate(disk, location);
                    self.envelope_state = EnvelopeState::Runaway {
                        envelope_mass: *envelope_mass,
                        accretion_rate: initial_rate,
                    };
                } else {
                    // Continue slow Kelvin-Helmholtz growth
                    let growth_rate = envelope::kelvin_helmholtz_growth_rate(
                        self.core_mass,
                        *envelope_mass,
                        self.physical_radius,
                        core_accretion_rate,
                    );

                    let dm = growth_rate.to_solar_masses_per_year() * dt.to_years();
                    *envelope_mass = *envelope_mass + Mass::from_solar_masses(dm);
                    self.envelope_mass = *envelope_mass;
                }
            }

            EnvelopeState::Runaway {
                envelope_mass,
                accretion_rate,
            } => {
                // Use precomputed values (computed before the match to avoid borrow conflicts)
                if opens_gap {
                    // Transition to final state (gap stops accretion)
                    self.envelope_state = EnvelopeState::Final {
                        envelope_mass: *envelope_mass,
                    };
                } else {
                    // Continue runaway accretion
                    // Rate limited by disk supply
                    *accretion_rate = envelope::supply_limited_accretion_rate(disk, location);

                    let dm = accretion_rate.to_solar_masses_per_year() * dt.to_years();
                    *envelope_mass = *envelope_mass + Mass::from_solar_masses(dm);
                    self.envelope_mass = *envelope_mass;

                    // Also check if we've consumed most of the local gas
                    if self.envelope_mass > local_gas_mass * 0.5 {
                        // Consumed most local gas, stop
                        self.envelope_state = EnvelopeState::Final {
                            envelope_mass: *envelope_mass,
                        };
                    }
                }
            }

            EnvelopeState::Final { .. } => {
                // No further evolution (gas disk dispersed or gap opened)
            }
        }

        // Update physical radius after envelope growth
        self.physical_radius = Self::estimate_radius(self.total_mass(), &self.composition);
    }

    /// Check if this body opens a gap in the gas disk.
    ///
    /// Gap opening requires the tidal torques to exceed viscous refilling.
    /// Two criteria must be met:
    ///
    /// **Thermal criterion**: q > (h/r)³ where q = M_p / M_*
    /// **Viscous criterion**: q > 40α × (h/r)⁵
    ///
    /// # Arguments
    /// * `disk` - Gas disk model
    /// * `stellar_mass` - Mass of central star
    ///
    /// # Returns
    /// true if both gap-opening criteria are satisfied
    ///
    /// # References
    /// - Crida et al. (2006) - "On the width and shape of gaps"
    /// - Kanagawa et al. (2015) - "Mass constraint for gap opening"
    pub fn opens_gap<D: DiskModel>(&self, disk: &D, stellar_mass: Mass) -> bool {
        let location = self.semi_major_axis;
        let h_over_r = disk.aspect_ratio(location);
        let alpha = disk.alpha();
        let q = self.total_mass().to_solar_masses() / stellar_mass.to_solar_masses();

        // Thermal criterion
        let thermal = q > h_over_r.powi(3);

        // Viscous criterion
        let viscous = q > 40.0 * alpha * h_over_r.powi(5);

        thermal && viscous
    }

    /// Estimate local gas mass available for accretion.
    ///
    /// Approximates the gas mass within the Hill sphere.
    ///
    /// # Arguments
    /// * `disk` - Gas disk model
    ///
    /// # Returns
    /// Local gas mass in solar masses
    pub(crate) fn local_gas_mass<D: DiskModel>(&self, disk: &D) -> Mass {
        let location = self.semi_major_axis;
        let r_hill = self.hill_radius(disk.stellar_mass());
        let sigma_gas = disk.surface_density(location);

        // Mass = π × R_H² × Σ
        let area_cm2 = std::f64::consts::PI * r_hill.to_cm().powi(2);
        let mass_g = area_cm2 * sigma_gas.to_grams_per_cm2();

        Mass::from_grams(mass_g)
    }

    /// Estimate physical radius from mass and composition.
    ///
    /// Uses simple mass-radius relations for different compositions.
    ///
    /// # Arguments
    /// * `mass` - Total mass
    /// * `composition` - Bulk composition
    ///
    /// # Returns
    /// Estimated physical radius
    pub(crate) fn estimate_radius(mass: Mass, composition: &Composition) -> Length {
        let m_earth = mass.to_earth_masses();

        // Terrestrial mass-radius relation: R ∝ M^0.27
        let r_rock = m_earth.powf(0.27);

        // Gas envelope inflation (simplified)
        let gas_fraction = composition.h_he_gas;
        let inflation = 1.0 + 2.0 * gas_fraction;

        Length::from_earth_radii(r_rock * inflation)
    }
}

impl DiscreteBody {
    /// Create a lightweight representation for Barnes-Hut tree calculations.
    ///
    /// The `Massive` trait requires `Copy`, which DiscreteBody can't implement.
    /// This creates a simple struct that can be used with BHTree.
    pub fn as_massive(&self) -> MassiveBody {
        MassiveBody {
            position: self.position,
            mass: self.total_mass(),
        }
    }
}

/// Lightweight representation of a body for Barnes-Hut tree calculations.
///
/// The `Massive` trait requires `Copy`, so we use this simple struct
/// for tree operations and keep the full `DiscreteBody` for storage.
#[derive(Debug, Clone, Copy)]
pub struct MassiveBody {
    pub position: Point2<f64>,
    pub mass: Mass,
}

impl nbody::arena_bhtree::Massive for MassiveBody {
    fn position(&self) -> Point2<f64> {
        self.position
    }

    fn mass(&self) -> Mass {
        self.mass
    }

    fn mass_solar(&self) -> f64 {
        self.mass.to_solar_masses()
    }
}
