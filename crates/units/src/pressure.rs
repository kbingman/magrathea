/// Pressure in dyn/cm² (g/(cm·s²))
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd, serde::Serialize, serde::Deserialize)]
pub struct Pressure(pub f64);

impl Pressure {
    pub fn from_dyn_per_cm2(value: f64) -> Self {
        Self(value)
    }

    pub fn to_dyn_per_cm2(&self) -> f64 {
        self.0
    }
}
