use serde::{Deserialize, Serialize};
use std::ops::{Add, Div, Mul, Sub};

/// Mass of the Sun in grams (1.98847 × 10³³ g)
pub const SOLAR_MASS_G: f64 = 1.98847e33;

/// Mass of the Earth in grams (5.972 × 10²⁷ g)
pub const EARTH_MASS_G: f64 = 5.972e27;

/// Mass of Jupiter in grams (1.898 × 10³⁰ g)
const JUPITER_MASS_G: f64 = 1.898e30;

/// A physical mass quantity using f64 precision.
///
/// The `Mass` struct represents mass values with solar masses as the base unit.
/// This choice allows for convenient astronomical calculations while maintaining
/// precision for planetary and stellar mass comparisons.
///
/// # Examples
///
/// ```rust
/// use units::Mass;
///
/// // Create masses using different units
/// let sun_mass = Mass::from_solar_masses(1.0);
/// let earth_mass = Mass::from_earth_masses(1.0);
/// let gram_mass = Mass::from_grams(1000.0);
///
/// // Convert between units
/// let earth_in_solar = earth_mass.to_solar_masses();
/// let sun_in_grams = sun_mass.to_grams();
/// ```
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd, Deserialize, Serialize)]
pub struct Mass(f64); // Base unit: Solar Masses

impl Mass {
    /// Creates a new `Mass` from a value in solar masses.
    ///
    /// This is the most direct constructor since solar masses are the base unit.
    ///
    /// # Arguments
    ///
    /// * `value` - The mass value in solar masses
    ///
    /// # Examples
    ///
    /// ```rust
    /// use units::Mass;
    ///
    /// let sun = Mass::from_solar_masses(1.0);
    /// let massive_star = Mass::from_solar_masses(25.0);
    /// ```
    pub fn from_solar_masses(value: f64) -> Self {
        Self(value)
    }

    /// Creates a new `Mass` from a value in Earth masses.
    ///
    /// Converts Earth masses to the internal solar mass representation.
    /// One solar mass is approximately 332,946 Earth masses.
    ///
    /// # Arguments
    ///
    /// * `value` - The mass value in Earth masses
    ///
    /// # Examples
    ///
    /// ```rust
    /// use units::Mass;
    ///
    /// let earth = Mass::from_earth_masses(1.0);
    /// let jupiter = Mass::from_earth_masses(317.8);
    /// let super_earth = Mass::from_earth_masses(5.0);
    /// ```
    pub fn from_earth_masses(value: f64) -> Self {
        Self(value * EARTH_MASS_G / SOLAR_MASS_G)
    }

    /// Creates a new `Mass` from a value in Jupiter masses.
    ///
    /// Jupiter mass is approximately 317.8 Earth masses or 0.000954 solar masses.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use units::Mass;
    ///
    /// let jupiter = Mass::from_jupiter_masses(1.0);
    /// let hot_jupiter = Mass::from_jupiter_masses(2.5);
    /// ```
    pub fn from_jupiter_masses(value: f64) -> Self {
        Self(value * JUPITER_MASS_G / SOLAR_MASS_G)
    }

    /// Creates a new `Mass` from a value in grams.
    ///
    /// Converts grams to the internal solar mass representation.
    /// This is useful for laboratory-scale or small astronomical objects.
    ///
    /// # Arguments
    ///
    /// * `value` - The mass value in grams
    ///
    /// # Examples
    ///
    /// ```rust
    /// use units::Mass;
    ///
    /// let kilogram = Mass::from_grams(1000.0);
    /// let moon = Mass::from_grams(7.342e25);
    /// let asteroid = Mass::from_grams(1.0e18);
    /// ```
    pub fn from_grams(value: f64) -> Self {
        Self(value / SOLAR_MASS_G)
    }

    pub fn from_kg(value: f64) -> Self {
        Self::from_grams(value * 1000.0)
    }

    /// Returns the mass value in solar masses.
    ///
    /// Since solar masses are the base unit, this simply returns the stored value.
    ///
    /// # Returns
    ///
    /// The mass in solar masses
    ///
    /// # Examples
    ///
    /// ```rust
    /// use units::Mass;
    ///
    /// let mass = Mass::from_earth_masses(332946.0);
    /// // Approximately 1 solar mass
    /// ```
    pub fn to_solar_masses(&self) -> f64 {
        self.0
    }

    /// Converts the mass to Earth masses.
    ///
    /// Returns the mass value expressed in terms of Earth's mass.
    /// One solar mass equals approximately 332,946 Earth masses.
    ///
    /// # Returns
    ///
    /// The mass in Earth masses
    ///
    /// # Examples
    ///
    /// ```rust
    /// use units::Mass;
    ///
    /// let sun = Mass::from_solar_masses(1.0);
    /// let earth_masses = sun.to_earth_masses();  // ~332,946 Earth masses
    ///
    /// let jupiter = Mass::from_earth_masses(317.8);
    /// assert_eq!(jupiter.to_earth_masses(), 317.8);
    /// ```
    pub fn to_earth_masses(&self) -> f64 {
        self.0 * SOLAR_MASS_G / EARTH_MASS_G
    }

    /// Converts the mass to Jupiter masses.
    ///
    /// # Returns
    ///
    /// The mass in Jupiter masses
    pub fn to_jupiter_masses(&self) -> f64 {
        self.0 * SOLAR_MASS_G / JUPITER_MASS_G
    }

    /// Converts the mass to grams.
    ///
    /// Returns the mass value in the CGS base unit of grams.
    /// This is useful for scientific calculations or when interfacing
    /// with systems that expect gram units.
    ///
    /// # Returns
    ///
    /// The mass in grams
    ///
    /// # Examples
    ///
    /// ```rust
    /// use units::Mass;
    ///
    /// let earth = Mass::from_earth_masses(1.0);
    /// let grams = earth.to_grams();  // ~5.972 × 10²⁷ grams
    ///
    /// let sun = Mass::from_solar_masses(1.0);
    /// let sun_grams = sun.to_grams();  // ~1.98847 × 10³³ grams
    /// ```
    pub fn to_grams(&self) -> f64 {
        self.0 * SOLAR_MASS_G
    }

    pub fn to_kg(&self) -> f64 {
        self.to_grams() / 1000.0
    }

    /// Raise to integer power
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

impl Add for Mass {
    type Output = Mass;

    fn add(self, rhs: Mass) -> Mass {
        Mass(self.0 + rhs.0)
    }
}

impl Sub for Mass {
    type Output = Mass;

    fn sub(self, rhs: Mass) -> Mass {
        Mass(self.0 - rhs.0)
    }
}

impl Mul<f64> for Mass {
    type Output = Mass;

    fn mul(self, rhs: f64) -> Mass {
        Mass(self.0 * rhs)
    }
}

impl Div<f64> for Mass {
    type Output = Mass;

    fn div(self, rhs: f64) -> Mass {
        Mass(self.0 / rhs)
    }
}

/// Division of Mass by Mass returns a dimensionless ratio
impl Div for Mass {
    type Output = f64;

    fn div(self, rhs: Mass) -> f64 {
        self.0 / rhs.0
    }
}

/// Allow f64 * Mass (commutative multiplication)
impl Mul<Mass> for f64 {
    type Output = Mass;

    fn mul(self, rhs: Mass) -> Mass {
        rhs * self
    }
}
