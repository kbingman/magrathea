/// Volume density in g/cm³
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd, serde::Serialize, serde::Deserialize)]
pub struct Density(pub f64);

impl Density {
    pub fn from_grams_per_cm3(value: f64) -> Self {
        Self(value)
    }

    pub fn to_grams_per_cm3(&self) -> f64 {
        self.0
    }
}
