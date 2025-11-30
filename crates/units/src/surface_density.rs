use serde::{Deserialize, Serialize};
use std::ops::{Add, Div, Mul, Sub};

/// One astronomical unit in centimeters
const AU_CM: f64 = 1.496e13;

/// A physical surface density quantity using f64 precision.
///
/// The `SurfaceDensity` struct represents surface mass density (mass per area)
/// with grams per square centimeter as the base unit. This is the standard unit
/// for protoplanetary disk surface density in astrophysical literature.
///
/// Surface density is a key parameter in planet formation, following power-law
/// profiles like Σ(r) = Σ₀ (r/r₀)^(-p) where typical values are:
/// - Σ₀ ~ 1700 g/cm² at 1 AU (Minimum Mass Solar Nebula)
/// - p ~ 1.5 (power law index)
///
/// # Examples
///
/// ```rust
/// use units::SurfaceDensity;
///
/// // Typical disk surface density at 1 AU
/// let sigma_1au = SurfaceDensity::from_grams_per_cm2(1700.0);
///
/// // Convert to kg/m²
/// let si_units = sigma_1au.to_kg_per_m2();
/// ```
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd, Deserialize, Serialize)]
pub struct SurfaceDensity(f64); // Base unit: g/cm²

impl SurfaceDensity {
    /// Creates a zero surface density value
    pub fn zero() -> Self {
        Self(0.0)
    }

    /// Creates a new `SurfaceDensity` from a value in grams per square centimeter.
    ///
    /// This is the most direct constructor since g/cm² is the base unit and
    /// the standard in astrophysical disk modeling literature.
    ///
    /// # Arguments
    ///
    /// * `value` - The surface density in g/cm²
    ///
    /// # Examples
    ///
    /// ```rust
    /// use units::SurfaceDensity;
    ///
    /// // Minimum Mass Solar Nebula at 1 AU
    /// let mmsn = SurfaceDensity::from_grams_per_cm2(1700.0);
    ///
    /// // Depleted disk
    /// let depleted = SurfaceDensity::from_grams_per_cm2(100.0);
    /// ```
    pub fn from_grams_per_cm2(value: f64) -> Self {
        Self(value)
    }

    /// Creates a new `SurfaceDensity` from a value in kilograms per square meter.
    ///
    /// Converts SI units to the internal CGS representation.
    /// 1 kg/m² = 0.1 g/cm²
    ///
    /// # Arguments
    ///
    /// * `value` - The surface density in kg/m²
    ///
    /// # Examples
    ///
    /// ```rust
    /// use units::SurfaceDensity;
    ///
    /// let sigma = SurfaceDensity::from_kg_per_m2(1000.0);  // 100 g/cm²
    /// ```
    pub fn from_kg_per_m2(value: f64) -> Self {
        // 1 kg/m² = 1000 g / 10000 cm² = 0.1 g/cm²
        Self(value * 0.1)
    }

    /// Creates a new `SurfaceDensity` from solar masses per square AU.
    ///
    /// This can be useful for very massive disks or when working with
    /// stellar-scale mass distributions.
    ///
    /// # Arguments
    ///
    /// * `value` - The surface density in M☉/AU²
    ///
    /// # Examples
    ///
    /// ```rust
    /// use units::SurfaceDensity;
    ///
    /// let massive_disk = SurfaceDensity::from_solar_masses_per_au2(1e-6);
    /// ```
    pub fn from_solar_masses_per_au2(value: f64) -> Self {
        // M☉ = 1.98847e33 g
        // AU² = (1.496e13)² cm²
        let solar_mass_g = 1.98847e33;
        let au2_cm2 = AU_CM * AU_CM;
        Self(value * solar_mass_g / au2_cm2)
    }

    /// Returns the surface density value in grams per square centimeter.
    ///
    /// # Returns
    ///
    /// The surface density in g/cm²
    ///
    /// # Examples
    ///
    /// ```rust
    /// use units::SurfaceDensity;
    ///
    /// let sigma = SurfaceDensity::from_kg_per_m2(1000.0);
    /// assert_eq!(sigma.to_grams_per_cm2(), 100.0);
    /// ```
    pub fn to_grams_per_cm2(&self) -> f64 {
        self.0
    }

    /// Converts the surface density to kilograms per square meter.
    ///
    /// # Returns
    ///
    /// The surface density in kg/m²
    ///
    /// # Examples
    ///
    /// ```rust
    /// use units::SurfaceDensity;
    ///
    /// let sigma = SurfaceDensity::from_grams_per_cm2(100.0);
    /// assert_eq!(sigma.to_kg_per_m2(), 1000.0);
    /// ```
    pub fn to_kg_per_m2(&self) -> f64 {
        // 1 g/cm² = 10 kg/m²
        self.0 * 10.0
    }

    /// Converts the surface density to solar masses per square AU.
    ///
    /// # Returns
    ///
    /// The surface density in M☉/AU²
    ///
    /// # Examples
    ///
    /// ```rust
    /// use units::SurfaceDensity;
    ///
    /// let sigma = SurfaceDensity::from_grams_per_cm2(1700.0);
    /// let solar_units = sigma.to_solar_masses_per_au2();
    /// ```
    pub fn to_solar_masses_per_au2(&self) -> f64 {
        let solar_mass_g = 1.98847e33;
        let au2_cm2 = AU_CM * AU_CM;
        self.0 * au2_cm2 / solar_mass_g
    }

    /// Raise to integer power
    pub fn powi(&self, n: i32) -> f64 {
        self.0.powi(n)
    }

    /// Power function
    pub fn powf(&self, n: f64) -> f64 {
        self.0.powf(n)
    }

    /// Applies a power-law scaling to the surface density.
    ///
    /// This is useful for modeling radial disk profiles: Σ(r) = Σ₀ (r/r₀)^(-p)
    ///
    /// # Arguments
    ///
    /// * `ratio` - The radius ratio (r/r₀)
    /// * `power` - The power law index (typically 0.5 to 2.0)
    ///
    /// # Returns
    ///
    /// The scaled surface density
    ///
    /// # Examples
    ///
    /// ```rust
    /// use units::SurfaceDensity;
    ///
    /// let sigma_1au = SurfaceDensity::from_grams_per_cm2(1700.0);
    ///
    /// // At 5 AU with p=1.5 power law
    /// let sigma_5au = sigma_1au.power_law_scaling(5.0, 1.5);
    /// // sigma_5au ≈ 1700 * 5^(-1.5) ≈ 152 g/cm²
    /// ```
    pub fn power_law_scaling(&self, ratio: f64, power: f64) -> Self {
        let scale = ratio.powf(-power);
        Self(self.0 * scale)
    }
}

impl Add for SurfaceDensity {
    type Output = SurfaceDensity;

    fn add(self, rhs: SurfaceDensity) -> SurfaceDensity {
        SurfaceDensity(self.0 + rhs.0)
    }
}

impl Sub for SurfaceDensity {
    type Output = SurfaceDensity;

    fn sub(self, rhs: SurfaceDensity) -> SurfaceDensity {
        SurfaceDensity(self.0 - rhs.0)
    }
}

impl Mul<f64> for SurfaceDensity {
    type Output = SurfaceDensity;

    fn mul(self, rhs: f64) -> SurfaceDensity {
        SurfaceDensity(self.0 * rhs)
    }
}

impl Div<f64> for SurfaceDensity {
    type Output = SurfaceDensity;

    fn div(self, rhs: f64) -> SurfaceDensity {
        SurfaceDensity(self.0 / rhs)
    }
}

/// Division of SurfaceDensity by SurfaceDensity returns a ratio
impl Div for SurfaceDensity {
    type Output = f64;

    fn div(self, rhs: SurfaceDensity) -> f64 {
        self.0 / rhs.0
    }
}

/// Allow f64 * SurfaceDensity (commutative multiplication)
impl Mul<SurfaceDensity> for f64 {
    type Output = SurfaceDensity;

    fn mul(self, rhs: SurfaceDensity) -> SurfaceDensity {
        rhs * self
    }
}
