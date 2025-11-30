use serde::{Deserialize, Serialize};
use std::ops::{Add, Div, Mul, Sub};

/// A physical temperature quantity using f64 precision.
///
/// The `Temperature` struct represents temperature with Kelvin as the base unit,
/// following the SI standard and astrophysical conventions. Kelvin is the natural
/// choice for astrophysics as it's an absolute scale starting at zero.
///
/// # Examples
///
/// ```rust
/// use units::Temperature;
///
/// // Create temperatures in different units
/// let disk_temp = Temperature::from_kelvin(280.0);
/// let earth_surface = Temperature::from_celsius(15.0);  // ~288 K
/// let room_temp = Temperature::from_fahrenheit(68.0);   // ~293 K
///
/// // Convert between units
/// let celsius = disk_temp.to_celsius();
/// let fahrenheit = disk_temp.to_fahrenheit();
///
/// // Use physical constants
/// let freezing = Temperature::water_freezing();  // 273.15 K
/// let boiling = Temperature::water_boiling();    // 373.15 K
/// ```
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd, Deserialize, Serialize)]
pub struct Temperature(f64); // Base unit: Kelvin

impl Temperature {
    /// Creates a new `Temperature` from a value in Kelvin.
    ///
    /// This is the most direct constructor since Kelvin is the base unit.
    ///
    /// # Arguments
    ///
    /// * `value` - The temperature in Kelvin
    ///
    /// # Examples
    ///
    /// ```rust
    /// use units::Temperature;
    ///
    /// let freezing = Temperature::from_kelvin(273.15);
    /// let room_temp = Temperature::from_kelvin(293.15);
    /// let disk_inner = Temperature::from_kelvin(1600.0);
    /// ```
    pub fn from_kelvin(value: f64) -> Self {
        Self(value)
    }

    /// Creates a new `Temperature` from a value in Celsius.
    ///
    /// Converts Celsius to Kelvin: K = °C + 273.15
    ///
    /// # Arguments
    ///
    /// * `value` - The temperature in Celsius
    ///
    /// # Examples
    ///
    /// ```rust
    /// use units::Temperature;
    ///
    /// let freezing = Temperature::from_celsius(0.0);     // 273.15 K
    /// let earth_avg = Temperature::from_celsius(15.0);   // 288.15 K
    /// let boiling = Temperature::from_celsius(100.0);    // 373.15 K
    /// ```
    pub fn from_celsius(value: f64) -> Self {
        Self(value + 273.15)
    }

    /// Creates a new `Temperature` from a value in Fahrenheit.
    ///
    /// Converts Fahrenheit to Kelvin: K = (°F - 32) × 5/9 + 273.15
    ///
    /// # Arguments
    ///
    /// * `value` - The temperature in Fahrenheit
    ///
    /// # Examples
    ///
    /// ```rust
    /// use units::Temperature;
    ///
    /// let freezing = Temperature::from_fahrenheit(32.0);  // 273.15 K
    /// let room_temp = Temperature::from_fahrenheit(68.0); // ~293.15 K
    /// ```
    pub fn from_fahrenheit(value: f64) -> Self {
        Self((value - 32.0) * 5.0 / 9.0 + 273.15)
    }

    /// Returns the temperature value in Kelvin.
    ///
    /// # Returns
    ///
    /// The temperature in Kelvin
    ///
    /// # Examples
    ///
    /// ```rust
    /// use units::Temperature;
    ///
    /// let temp = Temperature::from_celsius(0.0);
    /// assert!((temp.to_kelvin() - 273.15).abs() < 0.01);
    /// ```
    pub fn to_kelvin(&self) -> f64 {
        self.0
    }

    /// Converts the temperature to Celsius.
    ///
    /// # Returns
    ///
    /// The temperature in Celsius (°C = K - 273.15)
    ///
    /// # Examples
    ///
    /// ```rust
    /// use units::Temperature;
    ///
    /// let temp = Temperature::from_kelvin(273.15);
    /// assert!((temp.to_celsius() - 0.0).abs() < 0.01);
    /// ```
    pub fn to_celsius(&self) -> f64 {
        self.0 - 273.15
    }

    /// Converts the temperature to Fahrenheit.
    ///
    /// # Returns
    ///
    /// The temperature in Fahrenheit (°F = (K - 273.15) × 9/5 + 32)
    ///
    /// # Examples
    ///
    /// ```rust
    /// use units::Temperature;
    ///
    /// let temp = Temperature::from_kelvin(273.15);
    /// assert!((temp.to_fahrenheit() - 32.0).abs() < 0.01);
    /// ```
    pub fn to_fahrenheit(&self) -> f64 {
        (self.0 - 273.15) * 9.0 / 5.0 + 32.0
    }

    /// Water freezing point at 1 atm (273.15 K / 0°C).
    ///
    /// This is a fundamental physical constant useful across many domains.
    pub fn water_freezing() -> Self {
        Self::from_kelvin(273.15)
    }

    /// Water boiling point at 1 atm (373.15 K / 100°C).
    ///
    /// This is a fundamental physical constant useful across many domains.
    pub fn water_boiling() -> Self {
        Self::from_kelvin(373.15)
    }

    /// Raise to integer power
    pub fn powi(&self, n: i32) -> f64 {
        self.0.powi(n)
    }

    /// Power function
    pub fn powf(&self, n: f64) -> f64 {
        self.0.powf(n)
    }
}

impl Add for Temperature {
    type Output = Temperature;

    fn add(self, rhs: Temperature) -> Temperature {
        Temperature(self.0 + rhs.0)
    }
}

impl Sub for Temperature {
    type Output = Temperature;

    fn sub(self, rhs: Temperature) -> Temperature {
        Temperature(self.0 - rhs.0)
    }
}

impl Mul<f64> for Temperature {
    type Output = Temperature;

    fn mul(self, rhs: f64) -> Temperature {
        Temperature(self.0 * rhs)
    }
}

impl Div<f64> for Temperature {
    type Output = Temperature;

    fn div(self, rhs: f64) -> Temperature {
        Temperature(self.0 / rhs)
    }
}

/// Allow f64 * Temperature (commutative multiplication)
impl Mul<Temperature> for f64 {
    type Output = Temperature;

    fn mul(self, rhs: Temperature) -> Temperature {
        rhs * self
    }
}
