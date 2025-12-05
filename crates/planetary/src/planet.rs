//! Planet representation with two-tier classification
//!
//! Combines the physical mass regime (PlanetClass) with observable
//! expression (PlanetType) and includes envelope structure modeling
//! from Lopez & Fortney (2014).

use rand::Rng;
use serde::{Deserialize, Serialize};
use units::{Length, Mass};

#[cfg(feature = "tsify")]
use tsify_next::Tsify;

use crate::composition::Composition;
use crate::planet_class::PlanetClass;
use crate::planet_type::PlanetType;

/// Host star properties for planet characterization
///
/// Groups stellar parameters needed for calculating planet properties
/// like equilibrium temperature and type classification.
#[derive(Debug, Clone, Copy)]
pub struct HostStar {
    /// Stellar luminosity in solar luminosities (L☉)
    pub luminosity: f64,
    /// Stellar mass in solar masses (M☉)
    pub mass: f64,
}

impl HostStar {
    /// Create a new host star context
    pub fn new(luminosity: f64, mass: f64) -> Self {
        Self { luminosity, mass }
    }

    /// Solar values (L = 1 L☉, M = 1 M☉)
    pub fn solar() -> Self {
        Self {
            luminosity: 1.0,
            mass: 1.0,
        }
    }
}

/// A fully characterized planet
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
#[cfg_attr(feature = "tsify", derive(Tsify))]
#[cfg_attr(feature = "tsify", tsify(into_wasm_abi, from_wasm_abi))]
pub struct Planet {
    /// Planet mass
    pub mass: Mass,
    /// Planet radius
    pub radius: Length,
    /// Orbital semi-major axis
    pub semi_major_axis: Length,
    /// Orbital eccentricity
    pub eccentricity: f64,
    /// Orbital inclination (radians)
    pub inclination: f64,
    /// Bulk composition
    pub composition: Composition,
    /// Mass-based physical regime
    pub class: PlanetClass,
    /// Observable expression
    pub planet_type: PlanetType,
    /// Equilibrium temperature (K)
    pub equilibrium_temp: f64,
}

impl Planet {
    /// Create a planet with full characterization
    ///
    /// Automatically determines class from mass and type from environment.
    pub fn new(
        mass: Mass,
        radius: Length,
        semi_major_axis: Length,
        eccentricity: f64,
        inclination: f64,
        composition: Composition,
        host: HostStar,
    ) -> Self {
        let class = PlanetClass::from_earth_masses(mass.to_earth_masses());
        let sma_au = semi_major_axis.to_au();
        let equilibrium_temp = calculate_equilibrium_temp(sma_au, host.luminosity);
        let incident_flux = host.luminosity / sma_au.powi(2);

        let planet_type = PlanetType::from_environment(
            class,
            &composition,
            equilibrium_temp,
            incident_flux,
            mass.to_earth_masses(),
            host.mass,
            sma_au,
        );

        Self {
            mass,
            radius,
            semi_major_axis,
            eccentricity,
            inclination,
            composition,
            class,
            planet_type,
            equilibrium_temp,
        }
    }

    /// Create a planet from mass and orbital parameters, calculating radius
    ///
    /// Uses the appropriate M-R relation based on composition.
    pub fn from_mass(
        mass: Mass,
        semi_major_axis: Length,
        eccentricity: f64,
        inclination: f64,
        composition: Composition,
        host: HostStar,
        rng: &mut impl Rng,
    ) -> Self {
        let mass_earth = mass.to_earth_masses();
        let class = PlanetClass::from_earth_masses(mass_earth);

        // Calculate radius - use envelope model if H/He present
        let radius_earth = if composition.h_he_gas > 0.01 {
            radius_with_envelope(mass_earth, composition.h_he_gas, rng)
        } else {
            class.radius_from_mass(mass_earth, true, rng)
        };

        let radius = Length::from_earth_radii(radius_earth);

        Self::new(
            mass,
            radius,
            semi_major_axis,
            eccentricity,
            inclination,
            composition,
            host,
        )
    }

    /// Returns bulk density in g/cm³
    pub fn density(&self) -> f64 {
        let mass_kg = self.mass.to_kg();
        let radius_m = self.radius.to_m();
        let volume_m3 = (4.0 / 3.0) * std::f64::consts::PI * radius_m.powi(3);
        (mass_kg / volume_m3) / 1000.0
    }

    /// Returns surface gravity in m/s²
    pub fn surface_gravity(&self) -> f64 {
        const G: f64 = 6.674e-11;
        G * self.mass.to_kg() / self.radius.to_m().powi(2)
    }

    /// Returns escape velocity in m/s
    pub fn escape_velocity(&self) -> f64 {
        const G: f64 = 6.674e-11;
        (2.0 * G * self.mass.to_kg() / self.radius.to_m()).sqrt()
    }

    /// Returns orbital period in years
    pub fn orbital_period(&self, stellar_mass_solar: f64) -> f64 {
        let a_au = self.semi_major_axis.to_au();
        (a_au.powi(3) / stellar_mass_solar).sqrt()
    }

    /// Returns incident flux in Earth units (F⊕)
    pub fn incident_flux(&self, stellar_luminosity: f64) -> f64 {
        stellar_luminosity / self.semi_major_axis.to_au().powi(2)
    }

    /// Check if planet is in the habitable zone
    pub fn in_habitable_zone(&self, stellar_luminosity: f64) -> bool {
        let flux = self.incident_flux(stellar_luminosity);

        (0.36..=1.1).contains(&flux) // Conservative HZ bounds
    }
}

/// Calculate equilibrium temperature in Kelvin
///
/// T_eq = 278 × (L/a²)^0.25 K
pub fn calculate_equilibrium_temp(semi_major_axis_au: f64, stellar_luminosity: f64) -> f64 {
    278.0 * (stellar_luminosity / semi_major_axis_au.powi(2)).powf(0.25)
}

/// Calculate radius for planet with H/He envelope (sub-Neptune structure)
///
/// Uses simplified layered model from Lopez & Fortney (2014):
/// 1. Calculate core radius from core mass using rocky planet scaling
/// 2. Add envelope layer using isothermal H/He structure
///
/// # Arguments
/// * `total_mass_earth` - Total planet mass in Earth masses
/// * `envelope_fraction` - Fraction of total mass in H/He envelope
/// * `rng` - Random number generator for scatter
///
/// # Returns
/// Planet radius in Earth radii
pub fn radius_with_envelope(
    total_mass_earth: f64,
    envelope_fraction: f64,
    rng: &mut impl Rng,
) -> f64 {
    // Split into core and envelope masses
    let core_mass_earth = total_mass_earth * (1.0 - envelope_fraction);
    let envelope_mass_earth = total_mass_earth * envelope_fraction;

    // Calculate core radius using rocky planet scaling (R ∝ M^0.27)
    // Chen et al. (2017) Earth-like composition fit
    let r_core = 1.07 * core_mass_earth.powf(0.27);

    // Calculate envelope contribution using scaling from Lopez & Fortney (2014)
    // ΔR/R_core scales with (M_env/M_core)^α where α ~ 0.57
    let envelope_scaling = 0.57;
    let envelope_contribution = if core_mass_earth > 0.0 {
        r_core * (envelope_mass_earth / core_mass_earth).powf(envelope_scaling)
    } else {
        0.0
    };

    // Total radius = core + envelope
    let base_radius = r_core + envelope_contribution;

    // Add observational scatter (~10% for sub-Neptunes)
    let scatter = rng.random_range(0.9..1.1);
    base_radius * scatter
}

// =============================================================================
// Factory functions for Solar System analogs
// =============================================================================

/// Create an Earth analog
pub fn earth_analog() -> Planet {
    Planet::new(
        Mass::from_earth_masses(1.0),
        Length::from_earth_radii(1.0),
        Length::from_au(1.0),
        0.017,
        0.0,
        Composition::earth_like(),
        HostStar::solar(),
    )
}

/// Create a Jupiter analog
pub fn jupiter_analog() -> Planet {
    Planet::new(
        Mass::from_earth_masses(317.8),
        Length::from_earth_radii(11.2),
        Length::from_au(5.2),
        0.049,
        0.022,
        Composition::gas_giant(),
        HostStar::solar(),
    )
}

/// Create a Neptune analog
pub fn neptune_analog() -> Planet {
    Planet::new(
        Mass::from_earth_masses(17.1),
        Length::from_earth_radii(3.88),
        Length::from_au(30.1),
        0.009,
        0.031,
        Composition::ice_giant(),
        HostStar::solar(),
    )
}
