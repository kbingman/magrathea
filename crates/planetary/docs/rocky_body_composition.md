# Rocky Body Composition System

*Design Document for Stellar System Generator*

---

## Overview

This document describes the composition generation system for rocky bodies (planets, asteroids, moons) based on observational data from polluted white dwarfs. White dwarf atmospheric pollution provides direct chemical analysis of exoplanetary material—the same material that formed during main sequence stellar evolution.

**Key insight**: White dwarfs act as mass spectrometers for planetary material that formed billions of years earlier. The statistics we observe (mostly chondritic, some differentiated, occasional volatile-rich) directly inform main sequence rocky body generation.

---

## Observational Basis

Between 25-50% of white dwarfs show metal pollution in their atmospheres, indicating accretion of rocky planetary material. Studies have detected up to 21 different heavy elements in single systems (GD 362, WD 1145+017).

### Composition Categories from Observations

| Category | Occurrence | Description |
|----------|------------|-------------|
| **Chondritic/Bulk Earth** | ~70% | Primitive material matching CI chondrite or bulk Earth. Dominated by O, Mg, Si, Fe. |
| **Differentiated** | ~15% | Fragments from bodies with iron cores. Core-rich (high Fe/Mg) or mantle-rich (low Fe/Mg). |
| **Volatile-rich** | ~10% | Material from beyond snow line. Contains water (5-50%), elevated C, N. |
| **Exotic** | ~5% | Compositions without solar system analogs. Unusual formation or processing. |

### Key Diagnostic Ratios

| Ratio | Indicates | Typical Values |
|-------|-----------|----------------|
| **Fe/Mg** | Core vs mantle | Chondritic: ~1.9, Core-rich: >5, Mantle-rich: <0.5 |
| **Ca/Mg** | Refractory enrichment | High = high-T formation or crust |
| **Mn/Na** | Volatile loss history | High = post-nebula volatilization |

---

## Major Rock-Forming Elements

| Element | Mass Fraction (Chondrite) | Notes |
|---------|---------------------------|-------|
| Oxygen | ~46% | Bound in silicates/oxides |
| Iron | ~18% | Variable with differentiation |
| Silicon | ~17% | Silicate backbone |
| Magnesium | ~15% | Silicate component |
| Sulfur | ~2% | Volatile, depleted inner system |
| Calcium | ~1.3% | Refractory, enriched in crust |
| Aluminum | ~1.2% | Refractory, enriched in crust |
| Nickel | ~1.1% | Siderophile, follows iron |

---

## Reference Compositions

```rust
pub fn ci_chondrite() -> HashMap<Element, f64> {
    hashmap! {
        Element::Oxygen => 0.46,
        Element::Silicon => 0.107,
        Element::Magnesium => 0.097,
        Element::Iron => 0.182,
        Element::Calcium => 0.013,
        Element::Aluminum => 0.012,
        Element::Nickel => 0.011,
        Element::Sulfur => 0.054,
        Element::Sodium => 0.006,
    }
}

pub fn bulk_earth() -> HashMap<Element, f64> {
    hashmap! {
        Element::Oxygen => 0.44,
        Element::Silicon => 0.15,
        Element::Magnesium => 0.14,
        Element::Iron => 0.32,  // Higher due to core
        Element::Calcium => 0.017,
        Element::Aluminum => 0.016,
        Element::Nickel => 0.018,
        Element::Sulfur => 0.006,  // Depleted
    }
}

pub fn core_composition() -> HashMap<Element, f64> {
    hashmap! {
        Element::Iron => 0.85,
        Element::Nickel => 0.10,
        Element::Sulfur => 0.04,
        Element::Chromium => 0.005,
    }
}

pub fn mantle_composition() -> HashMap<Element, f64> {
    hashmap! {
        Element::Oxygen => 0.44,
        Element::Silicon => 0.21,
        Element::Magnesium => 0.23,
        Element::Iron => 0.06,
        Element::Calcium => 0.025,
        Element::Aluminum => 0.023,
    }
}
```

---

## Formation Context

### Condensation Temperature Effects

| Temperature | Disk Region | Volatile Retention |
|-------------|-------------|-------------------|
| >1200 K | Inner (<0.5 AU) | 10% |
| 900-1200 K | Terrestrial zone | 30% |
| 600-900 K | Outer terrestrial | 60% |
| 400-600 K | Main belt | 85% |
| <400 K | Beyond snow line | 100% + ices |

```rust
fn volatile_retention(condensation_temp: Kelvin) -> f64 {
    match condensation_temp.0 {
        t if t > 1200.0 => 0.1,
        t if t > 900.0 => 0.3,
        t if t > 600.0 => 0.6,
        t if t > 400.0 => 0.85,
        _ => 1.0,
    }
}
```

### Snow Line and Water Content

```rust
fn water_content(formation_distance: Au, snow_line: Au, rng: &mut impl Rng) -> f64 {
    if formation_distance < snow_line {
        0.0
    } else {
        let ratio = formation_distance.0 / snow_line.0;
        match ratio {
            r if r < 1.5 => rng.gen_range(0.05..0.15),
            r if r < 3.0 => rng.gen_range(0.10..0.30),
            _ => rng.gen_range(0.20..0.50),
        }
    }
}
```

---

## Differentiation

### When Bodies Differentiate

Bodies >~300 km radius typically differentiate due to radiogenic heating (²⁶Al).

```rust
fn is_differentiated(radius_km: f64, rng: &mut impl Rng) -> bool {
    match radius_km {
        r if r < 100.0 => false,
        r if r > 500.0 => rng.gen_bool(0.98),
        r => rng.gen_bool((r - 100.0) / 500.0),
    }
}
```

### Core Mass Fractions

| Body | Core Fraction | Notes |
|------|---------------|-------|
| Mercury | ~0.75 | Mantle stripping or inner formation |
| Earth | ~0.325 | Typical for large terrestrial |
| Mars | ~0.21 | Smaller body |
| Moon | ~0.02 | Formed from Earth's mantle |

```rust
fn sample_core_fraction(mass: EarthMass, distance: Au, rng: &mut impl Rng) -> f64 {
    let base = match mass.0 {
        m if m > 0.5 => rng.gen_range(0.28..0.40),
        m if m > 0.1 => rng.gen_range(0.18..0.35),
        _ => rng.gen_range(0.15..0.50),
    };
    
    // Mercury effect: close formation concentrates iron
    let proximity_boost = if distance.0 < 0.5 {
        rng.gen_range(0.0..0.3)
    } else { 0.0 };
    
    (base + proximity_boost).min(0.80)
}
```

---

## Collisional Fragments

```rust
pub enum PlanetaryLayer { Core, LowerMantle, UpperMantle, Crust }

fn sample_fragment_origin(rng: &mut impl Rng) -> PlanetaryLayer {
    match rng.gen::<f64>() {
        x if x < 0.25 => PlanetaryLayer::Core,       // M-type
        x if x < 0.70 => PlanetaryLayer::UpperMantle, // S-type
        x if x < 0.90 => PlanetaryLayer::LowerMantle,
        _ => PlanetaryLayer::Crust,                   // V-type
    }
}
```

---

## Asteroid Belt Population

| Class | Fraction | Composition |
|-------|----------|-------------|
| Carbonaceous (C) | ~75% | Primitive, water-bearing |
| Siliceous (S) | ~17% | Differentiated silicates |
| Metallic (M) | ~5% | Core fragments |
| Basaltic (V) | ~3% | Crustal material |

---

## Stellar Metallicity Scaling

Refractory elements scale with stellar [Fe/H]:

```rust
fn scale_for_metallicity(base: &HashMap<Element, f64>, feh: f64) -> HashMap<Element, f64> {
    let scale = 10.0_f64.powf(feh);
    base.iter().map(|(elem, frac)| {
        let scaled = match elem {
            Element::Iron | Element::Magnesium | Element::Silicon |
            Element::Calcium | Element::Aluminum | Element::Nickel => frac * scale,
            Element::Oxygen => frac * scale.powf(0.5),
            _ => *frac,
        };
        (*elem, scaled)
    }).collect()
}
```

---

## Complete Generation Pipeline

```rust
pub struct RockyBodyComposition {
    pub elements: HashMap<Element, f64>,
    pub formation_context: FormationContext,
    pub processing_history: ProcessingHistory,
}

impl RockyBodyComposition {
    pub fn generate_for_planet(
        mass: EarthMass,
        formation_distance: Au,
        star: &Star,
        rng: &mut impl Rng,
    ) -> Self {
        let snow_line = star.snow_line_distance();
        let condensation_temp = disk_temperature(formation_distance, star);
        
        let ctx = FormationContext {
            formation_distance,
            condensation_temp,
            inside_snow_line: formation_distance < snow_line,
            stellar_metallicity: star.metallicity(),
        };
        
        let radius_km = estimate_radius_km(mass);
        let differentiated = radius_km > 300.0 && rng.gen_bool(0.95);
        let core_fraction = differentiated.then(|| 
            sample_core_fraction(mass, formation_distance, rng));
        
        let history = ProcessingHistory {
            differentiated,
            core_fraction,
            volatile_depleted: ctx.inside_snow_line,
            is_collisional_fragment: false,
            fragment_origin: None,
        };
        
        let elements = Self::compute_elements(&ctx, &history, star, rng);
        Self { elements, formation_context: ctx, processing_history: history }
    }
}
```

---

## References

- Trierweiler et al. (2023) - "A Chondritic Solar Neighborhood"
- Putirka & Xu (2021) - "Polluted white dwarfs reveal exotic mantle rock types"
- Harrison et al. (2021) - Population-level differentiation analysis
- Xu et al. (2014, 2019) - Detailed polluted white dwarf abundances
