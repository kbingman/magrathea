use crate::length::AU_TO_CM;
use crate::time::SECONDS_PER_YEAR;
use std::{
    f32::consts::TAU,
    ops::{Add, Div, Mul, Sub},
};

// Convert between AU/year and cm/sec
pub const AU_YEAR_TO_CM_SEC: f64 = AU_TO_CM / SECONDS_PER_YEAR;

// Orbital velocity scaling constants
// Earth's orbital velocity: ~29.78 km/s = 6.28 AU/year in Magrathea units
pub const EARTH_ORBITAL_VELOCITY_AU_PER_YEAR: f32 = TAU;

// Unit conversion from normalized orbital velocity to AU/year
// This converts from sqrt(GM/r) in natural units to AU/year matching real orbital velocities
pub const ORBITAL_VELOCITY_SCALE_FACTOR: f32 = EARTH_ORBITAL_VELOCITY_AU_PER_YEAR;

/// Calculate circular orbital velocity in AU/year units for the Magrathea system
///
/// This function converts from the normalized Kepler velocity sqrt(GM/r) to the proper
/// AU/year units used throughout the Magrathea codebase. The result matches real-world
/// orbital velocities when converted to km/s.
///
/// # Arguments
/// * `stellar_mass` - Mass of the central star in solar masses
/// * `radius` - Orbital radius in AU
///
/// # Returns
/// Circular orbital velocity in AU/year
///
/// # Examples
/// ```
/// use units::velocity::circular_orbital_velocity;
///
/// // Earth's orbital velocity around a 1 solar mass star at 1 AU
/// let earth_v = circular_orbital_velocity(1.0, 1.0);
/// assert_eq!(earth_v, 6.2831855); // AU/year, equivalent to ~30 km/s
///
/// // Jupiter's orbital velocity
/// let jupiter_v = circular_orbital_velocity(1.0, 5.2);
/// // Should be ~13.07 AU/year, equivalent to ~13 km/s
/// ```
pub fn circular_orbital_velocity(stellar_mass: f32, radius: f32) -> f32 {
    (stellar_mass / radius).sqrt() * ORBITAL_VELOCITY_SCALE_FACTOR
}

#[derive(Debug, Clone, Copy, PartialEq, PartialOrd, serde::Serialize, serde::Deserialize)]
pub struct Velocity(f64); // Base unit: AU/year

impl Velocity {
    pub fn from_au_per_year(value: f64) -> Self {
        Self(value)
    }

    pub fn from_cm_per_sec(value: f64) -> Self {
        Self(value / AU_YEAR_TO_CM_SEC)
    }

    pub fn from_meters_per_sec(value: f64) -> Self {
        Self::from_cm_per_sec(value * 100.0)
    }

    pub fn to_au_per_year(&self) -> f64 {
        self.0
    }

    pub fn to_cm_per_sec(&self) -> f64 {
        self.0 * AU_YEAR_TO_CM_SEC
    }

    pub fn to_meters_per_sec(&self) -> f64 {
        self.to_cm_per_sec() / 100.0
    }
}

impl Add for Velocity {
    type Output = Velocity;

    fn add(self, rhs: Velocity) -> Velocity {
        Velocity(self.0 + rhs.0)
    }
}

impl Sub for Velocity {
    type Output = Velocity;

    fn sub(self, rhs: Velocity) -> Velocity {
        Velocity(self.0 - rhs.0)
    }
}

impl Mul<f64> for Velocity {
    type Output = Velocity;

    fn mul(self, rhs: f64) -> Velocity {
        Velocity(self.0 * rhs)
    }
}

impl Div<f64> for Velocity {
    type Output = Velocity;

    fn div(self, rhs: f64) -> Velocity {
        Velocity(self.0 / rhs)
    }
}
