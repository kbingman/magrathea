use serde::{Deserialize, Serialize};

/// RGB color representation for stellar objects
///
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct StellarColor {
    pub r: u8,
    pub g: u8,
    pub b: u8,
}

impl StellarColor {
    // Tanner Helland blackbody-to-RGB algorithm coefficients
    // Reference: https://tannerhelland.com/2012/09/18/convert-temperature-rgb-algorithm-code.html
    //
    // These are curve-fit coefficients for approximating the Planckian locus
    // (the path of blackbody colors through CIE color space)

    /// Red channel coefficient for hot stars (temp > 6600K)
    const RED_COEFF: f64 = 329.698727446;
    const RED_EXP: f64 = -0.1332047592;

    /// Green channel coefficients
    const GREEN_COOL_COEFF: f64 = 99.4708025861;
    const GREEN_COOL_OFFSET: f64 = -161.1195681661;
    const GREEN_HOT_COEFF: f64 = 288.1221695283;
    const GREEN_HOT_EXP: f64 = -0.0755148492;

    /// Blue channel coefficients for warm stars (1900K < temp < 6600K)
    const BLUE_COEFF: f64 = 138.5177312231;
    const BLUE_OFFSET: f64 = -305.0447927307;

    /// Temperature thresholds (in units of temp/100)
    const TEMP_HOT_THRESHOLD: f64 = 66.0; // 6600K - transition point
    const TEMP_BLUE_CUTOFF: f64 = 19.0; // 1900K - below this, no blue

    /// Desaturation blend factor for realistic "whitish" stellar appearance
    const DESATURATION_BLEND: f64 = 0.3;

    /// Valid temperature range for stellar objects (Kelvin)
    const MIN_TEMP: f64 = 1000.0;
    const MAX_TEMP: f64 = 40000.0;

    pub fn new(r: u8, g: u8, b: u8) -> Self {
        Self { r, g, b }
    }

    /// Convert blackbody temperature to RGB color
    ///
    /// Uses an approximation of the Planckian locus for stellar temperatures.
    /// Based on the Tanner Helland algorithm, adjusted for stellar ranges.
    ///
    /// # Arguments
    /// * `temperature` - Temperature in Kelvin (clamped to 1000K-40000K)
    ///
    /// # References
    /// - Tanner Helland (2012) - "How to Convert Temperature to RGB"
    pub fn from_temperature(temperature: f64) -> Self {
        let temp = temperature.clamp(Self::MIN_TEMP, Self::MAX_TEMP) / 100.0;

        let r = match temp {
            t if t <= Self::TEMP_HOT_THRESHOLD => 255.0,
            t => (Self::RED_COEFF * (t - 60.0).powf(Self::RED_EXP)).clamp(0.0, 255.0),
        };

        let g = match temp {
            t if t <= Self::TEMP_HOT_THRESHOLD => {
                (Self::GREEN_COOL_COEFF * t.ln() + Self::GREEN_COOL_OFFSET).clamp(0.0, 255.0)
            }
            t => (Self::GREEN_HOT_COEFF * (t - 60.0).powf(Self::GREEN_HOT_EXP)).clamp(0.0, 255.0),
        };

        let b = match temp {
            t if t >= Self::TEMP_HOT_THRESHOLD => 255.0,
            t if t <= Self::TEMP_BLUE_CUTOFF => 0.0,
            t => (Self::BLUE_COEFF * (t - 10.0).ln() + Self::BLUE_OFFSET).clamp(0.0, 255.0),
        };

        // Apply desaturation for realistic "whitish" stellar colors
        // Stars appear more white to the eye than pure blackbody would suggest
        let avg = (r + g + b) / 3.0;
        let r = r + (avg - r) * Self::DESATURATION_BLEND;
        let g = g + (avg - g) * Self::DESATURATION_BLEND;
        let b = b + (avg - b) * Self::DESATURATION_BLEND;

        Self {
            r: r.round() as u8,
            g: g.round() as u8,
            b: b.round() as u8,
        }
    }

    /// Returns the color as a hex string (e.g., "#FF9944")
    pub fn to_hex(&self) -> String {
        format!("#{:02X}{:02X}{:02X}", self.r, self.g, self.b)
    }

    /// Parse a hex color string (e.g., "#FF9944" or "FF9944")
    pub fn from_hex(s: &str) -> Result<Self, String> {
        let s = s.strip_prefix('#').unwrap_or(s);

        if s.len() != 6 {
            return Err(format!("Invalid hex color length: {}", s));
        }

        let r = u8::from_str_radix(&s[0..2], 16)
            .map_err(|_| format!("Invalid red component: {}", &s[0..2]))?;
        let g = u8::from_str_radix(&s[2..4], 16)
            .map_err(|_| format!("Invalid green component: {}", &s[2..4]))?;
        let b = u8::from_str_radix(&s[4..6], 16)
            .map_err(|_| format!("Invalid blue component: {}", &s[4..6]))?;

        Ok(Self { r, g, b })
    }
}
