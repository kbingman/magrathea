use serde::{Deserialize, Serialize};
use std::ops::{Add, Div, Mul, Sub};

pub const AU_TO_CM: f64 = 1.496e13;
pub const AU_TO_M: f64 = 1.496e11;
pub const AU_TO_KM: f64 = 1.496e8;
pub const AU_TO_EARTH_RADIUS: f64 = 23481.4;
pub const MICRON_TO_CM: f64 = 1e-4;

/// Jupiter radius in Earth radii: 1 R_J = 11.209 R⊕
pub const JUPITER_TO_EARTH_RADII: f64 = 11.209;
/// AU to Jupiter radii
pub const AU_TO_JUPITER_RADIUS: f64 = AU_TO_EARTH_RADIUS / JUPITER_TO_EARTH_RADII;

/// Solar radius in AU: 1 R☉ = 0.00465047 AU
pub const SOLAR_RADIUS_AU: f64 = 1.0 / 215.032;
/// AU to solar radii
pub const AU_TO_SOLAR_RADIUS: f64 = 1.0 / SOLAR_RADIUS_AU;

/// A physical length quantity using f64 precision.
///
/// The `Length` struct represents length values with astronomical units (AU) as the base unit.
/// This is the natural choice for stellar and planetary system calculations.
///
/// # Examples
///
/// ```rust
/// use units::Length;
///
/// // Create lengths using different units
/// let earth_orbit = Length::from_au(1.0);
/// let earth_radius = Length::from_earth_radii(1.0);
/// let distance = Length::from_km(1000.0);
///
/// // Convert between units
/// let orbit_in_km = earth_orbit.to_km();
/// ```
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd, Deserialize, Serialize)]
#[serde(transparent)]
pub struct Length(f64); // Base unit: AU

impl Length {
    /// Creates a zero length value
    pub fn zero() -> Self {
        Self(0.0)
    }

    /// Creates a new `Length` from a value in astronomical units.
    pub fn from_au(value: f64) -> Self {
        Self(value)
    }

    /// Creates a new `Length` from a value in Earth radii.
    pub fn from_earth_radii(value: f64) -> Self {
        Self(value / AU_TO_EARTH_RADIUS)
    }

    /// Creates a new `Length` from a value in Jupiter radii.
    pub fn from_jupiter_radii(value: f64) -> Self {
        Self(value / AU_TO_JUPITER_RADIUS)
    }

    /// Creates a new `Length` from a value in solar radii.
    pub fn from_solar_radii(value: f64) -> Self {
        Self(value * SOLAR_RADIUS_AU)
    }

    /// Creates a new `Length` from a value in kilometers.
    pub fn from_km(value: f64) -> Self {
        Self(value / AU_TO_KM)
    }

    /// Creates a new `Length` from a value in meters.
    pub fn from_meters(value: f64) -> Self {
        Self(value / AU_TO_M)
    }

    /// Creates a new `Length` from a value in centimeters.
    pub fn from_cm(value: f64) -> Self {
        Self(value / AU_TO_CM)
    }

    /// Creates a new `Length` from a value in microns.
    pub fn from_microns(value: f64) -> Self {
        Self((value * MICRON_TO_CM) / AU_TO_CM)
    }

    /// Returns the length in astronomical units.
    pub fn to_au(&self) -> f64 {
        self.0
    }

    /// Converts the length to Earth radii.
    pub fn to_earth_radii(&self) -> f64 {
        self.0 * AU_TO_EARTH_RADIUS
    }

    /// Converts the length to Jupiter radii.
    pub fn to_jupiter_radii(&self) -> f64 {
        self.0 * AU_TO_JUPITER_RADIUS
    }

    /// Converts the length to solar radii.
    pub fn to_solar_radii(&self) -> f64 {
        self.0 * AU_TO_SOLAR_RADIUS
    }

    /// Converts the length to kilometers.
    pub fn to_km(&self) -> f64 {
        self.0 * AU_TO_KM
    }

    /// Converts the length to meters.
    pub fn to_m(&self) -> f64 {
        self.0 * AU_TO_M
    }

    /// Converts the length to centimeters.
    pub fn to_cm(&self) -> f64 {
        self.0 * AU_TO_CM
    }

    /// Converts the length to microns.
    pub fn to_microns(&self) -> f64 {
        self.0 * AU_TO_CM / MICRON_TO_CM
    }

    /// Returns the minimum of two lengths.
    pub fn min(self, other: Self) -> Self {
        if self.0 < other.0 {
            self
        } else {
            other
        }
    }

    /// Returns the maximum of two lengths.
    pub fn max(self, other: Self) -> Self {
        if self.0 > other.0 {
            self
        } else {
            other
        }
    }

    /// Raise to integer power (returns dimensionless f64 for dimensional consistency)
    pub fn powi(&self, n: i32) -> f64 {
        self.0.powi(n)
    }

    /// Power function
    pub fn powf(&self, n: f64) -> f64 {
        self.0.powf(n)
    }

    /// Natural logarithm
    pub fn ln(&self) -> f64 {
        self.0.ln()
    }

    /// Square root
    pub fn sqrt(&self) -> f64 {
        self.0.sqrt()
    }
}

impl Add for Length {
    type Output = Length;

    fn add(self, rhs: Length) -> Length {
        Length(self.0 + rhs.0)
    }
}

impl Sub for Length {
    type Output = Length;

    fn sub(self, rhs: Length) -> Length {
        Length(self.0 - rhs.0)
    }
}

impl Mul<f64> for Length {
    type Output = Length;

    fn mul(self, rhs: f64) -> Length {
        Length(self.0 * rhs)
    }
}

impl Mul for Length {
    type Output = Length;

    fn mul(self, rhs: Self) -> Length {
        Length(self.0 * rhs.0)
    }
}

impl Div<f64> for Length {
    type Output = Length;

    fn div(self, rhs: f64) -> Length {
        Length(self.0 / rhs)
    }
}

/// Division of Length by Length returns a dimensionless ratio
impl Div for Length {
    type Output = f64;

    fn div(self, rhs: Self) -> f64 {
        self.0 / rhs.0
    }
}

impl Div<&Length> for Length {
    type Output = f64;

    fn div(self, rhs: &Length) -> f64 {
        self.0 / rhs.0
    }
}

/// Allow f64 * Length (commutative multiplication)
impl Mul<Length> for f64 {
    type Output = Length;

    fn mul(self, rhs: Length) -> Length {
        rhs * self
    }
}
