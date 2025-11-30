use serde::{Deserialize, Serialize};
use std::ops::{Add, Div, Mul, Sub};

/// A physical volume density (mass per volume) quantity using f64 precision.
///
/// The `VolumeDensity` struct represents bulk density with grams per cubic centimeter
/// as the base unit, following the CGS convention used in astrophysics.
///
/// Typical planetary bulk densities:
/// - Iron: ~7.9 g/cm³
/// - Rock (silicates): ~3.3 g/cm³
/// - Water ice: ~0.9 g/cm³
/// - Hydrogen/Helium gas: ~0.001-0.1 g/cm³ (depends on pressure)
/// - Earth: ~5.5 g/cm³ (differentiated iron core + rocky mantle)
/// - Jupiter: ~1.3 g/cm³ (compressed gas)
///
/// # Examples
///
/// ```rust
/// use units::volume_density::VolumeDensity;
///
/// // Earth's bulk density
/// let earth_density = VolumeDensity::from_grams_per_cm3(5.5);
///
/// // Pure iron core
/// let iron = VolumeDensity::from_grams_per_cm3(7.9);
///
/// // Gas giant
/// let jupiter = VolumeDensity::from_grams_per_cm3(1.3);
/// ```
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd, Deserialize, Serialize)]
pub struct VolumeDensity(f64); // Base unit: g/cm³

impl VolumeDensity {
    /// Creates a new `VolumeDensity` from a value in grams per cubic centimeter.
    ///
    /// This is the most direct constructor since g/cm³ is the base unit and
    /// the standard in planetary science literature.
    ///
    /// # Arguments
    ///
    /// * `value` - The volume density in g/cm³
    ///
    /// # Examples
    ///
    /// ```rust
    /// use units::volume_density::VolumeDensity;
    ///
    /// let earth = VolumeDensity::from_grams_per_cm3(5.5);       // Earth
    /// let iron = VolumeDensity::from_grams_per_cm3(7.9);        // Pure iron
    /// let rock = VolumeDensity::from_grams_per_cm3(3.3);        // Silicate rock
    /// let ice = VolumeDensity::from_grams_per_cm3(0.9);         // Water ice
    /// ```
    pub fn from_grams_per_cm3(value: f64) -> Self {
        Self(value)
    }

    /// Creates a new `VolumeDensity` from a value in kilograms per cubic meter.
    ///
    /// Converts SI units to the internal CGS representation.
    /// 1 kg/m³ = 0.001 g/cm³
    ///
    /// # Arguments
    ///
    /// * `value` - The volume density in kg/m³
    ///
    /// # Examples
    ///
    /// ```rust
    /// use units::volume_density::VolumeDensity;
    ///
    /// let water = VolumeDensity::from_kg_per_m3(1000.0);  // 1.0 g/cm³
    /// ```
    pub fn from_kg_per_m3(value: f64) -> Self {
        // 1 kg/m³ = 1000 g / 1000000 cm³ = 0.001 g/cm³
        Self(value * 0.001)
    }

    /// Returns the volume density value in grams per cubic centimeter.
    ///
    /// # Returns
    ///
    /// The volume density in g/cm³
    ///
    /// # Examples
    ///
    /// ```rust
    /// use units::volume_density::VolumeDensity;
    ///
    /// let density = VolumeDensity::from_kg_per_m3(5500.0);
    /// assert_eq!(density.to_grams_per_cm3(), 5.5);  // Earth density
    /// ```
    pub fn to_grams_per_cm3(&self) -> f64 {
        self.0
    }

    /// Converts the volume density to kilograms per cubic meter.
    ///
    /// # Returns
    ///
    /// The volume density in kg/m³
    ///
    /// # Examples
    ///
    /// ```rust
    /// use units::volume_density::VolumeDensity;
    ///
    /// let earth = VolumeDensity::from_grams_per_cm3(5.5);
    /// assert_eq!(earth.to_kg_per_m3(), 5500.0);
    /// ```
    pub fn to_kg_per_m3(&self) -> f64 {
        // 1 g/cm³ = 1000 kg/m³
        self.0 * 1000.0
    }

    /// Common material density constants for reference.
    ///
    /// Returns density in g/cm³ for common planetary materials.
    pub fn iron() -> Self {
        Self::from_grams_per_cm3(7.9)
    }

    pub fn silicate_rock() -> Self {
        Self::from_grams_per_cm3(3.3)
    }

    pub fn water_ice() -> Self {
        Self::from_grams_per_cm3(0.92)
    }

    pub fn methane_ice() -> Self {
        Self::from_grams_per_cm3(0.47)
    }

    /// Calculates a weighted average density from multiple components.
    ///
    /// This is useful for computing bulk density from compositional models.
    ///
    /// # Arguments
    ///
    /// * `components` - Slice of (density, mass_fraction) tuples
    ///
    /// # Returns
    ///
    /// The weighted average density, or None if mass fractions don't sum to ~1
    ///
    /// # Examples
    ///
    /// ```rust
    /// use units::volume_density::VolumeDensity;
    ///
    /// // Earth-like planet: 32% iron core, 68% silicate mantle
    /// let components = vec![
    ///     (VolumeDensity::iron(), 0.32),
    ///     (VolumeDensity::silicate_rock(), 0.68),
    /// ];
    ///
    /// let bulk_density = VolumeDensity::weighted_average(&components).unwrap();
    /// // Should be around 5.1 g/cm³ (uncompressed)
    /// ```
    pub fn weighted_average(components: &[(VolumeDensity, f64)]) -> Option<Self> {
        let total_fraction: f64 = components.iter().map(|(_, frac)| frac).sum();

        // Check that fractions sum to approximately 1
        if (total_fraction - 1.0).abs() > 0.01 {
            return None;
        }

        let weighted_sum: f64 = components
            .iter()
            .map(|(density, frac)| density.to_grams_per_cm3() * frac)
            .sum();

        Some(Self::from_grams_per_cm3(weighted_sum))
    }
}

impl Add for VolumeDensity {
    type Output = VolumeDensity;

    fn add(self, rhs: VolumeDensity) -> VolumeDensity {
        VolumeDensity(self.0 + rhs.0)
    }
}

impl Sub for VolumeDensity {
    type Output = VolumeDensity;

    fn sub(self, rhs: VolumeDensity) -> VolumeDensity {
        VolumeDensity(self.0 - rhs.0)
    }
}

impl Mul<f64> for VolumeDensity {
    type Output = VolumeDensity;

    fn mul(self, rhs: f64) -> VolumeDensity {
        VolumeDensity(self.0 * rhs)
    }
}

impl Div<f64> for VolumeDensity {
    type Output = VolumeDensity;

    fn div(self, rhs: f64) -> VolumeDensity {
        VolumeDensity(self.0 / rhs)
    }
}

/// Allow f64 * VolumeDensity (commutative multiplication)
impl Mul<VolumeDensity> for f64 {
    type Output = VolumeDensity;

    fn mul(self, rhs: VolumeDensity) -> VolumeDensity {
        rhs * self
    }
}
