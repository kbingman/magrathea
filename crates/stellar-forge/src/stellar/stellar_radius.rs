//! Stellar radius type with solar radii serialization.

use serde::{Deserialize, Serialize};
use units::Length;

/// A stellar radius that serializes to solar radii.
///
/// This is a thin wrapper around `Length` that ensures consistent
/// serialization format for stellar object radii across the API.
///
/// # Example
/// ```
/// use stellar_forge::stellar::stellar_radius::StellarRadius;
/// use units::Length;
///
/// let radius = StellarRadius::new(Length::from_solar_radii(1.0));
/// assert_eq!(radius.to_solar_radii(), 1.0);
/// ```
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct StellarRadius(Length);

impl StellarRadius {
    /// Create a new stellar radius from a Length.
    pub fn new(length: Length) -> Self {
        Self(length)
    }

    /// Create a stellar radius from solar radii.
    pub fn from_solar_radii(value: f64) -> Self {
        Self(Length::from_solar_radii(value))
    }

    /// Get the radius in solar radii.
    pub fn to_solar_radii(&self) -> f64 {
        self.0.to_solar_radii()
    }

    /// Create a stellar radius from kilometers.
    pub fn from_km(value: f64) -> Self {
        Self(Length::from_km(value))
    }

    /// Get the underlying Length.
    pub fn as_length(&self) -> Length {
        self.0
    }

    /// Get the radius in kilometers.
    pub fn to_km(&self) -> f64 {
        self.0.to_km()
    }
}

impl From<Length> for StellarRadius {
    fn from(length: Length) -> Self {
        Self::new(length)
    }
}

impl From<StellarRadius> for Length {
    fn from(radius: StellarRadius) -> Self {
        radius.0
    }
}

impl Serialize for StellarRadius {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        self.to_solar_radii().serialize(serializer)
    }
}

impl<'de> Deserialize<'de> for StellarRadius {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        let solar_radii = f64::deserialize(deserializer)?;
        Ok(Self::from_solar_radii(solar_radii))
    }
}
