# Remnant Stellar Systems

*Design Document for Stellar System Generator*

---

## Overview

This document describes planetary system generation around stellar remnants: white dwarfs, neutron stars, and pulsars. These systems require separate generation pathways from main sequence systems due to fundamentally different formation mechanisms.

**Key insight**: The first confirmed exoplanets were found around pulsar PSR B1257+12 in 1992—predating 51 Pegasi b. Remnant systems are rare but real, and their "weirdness" provides excellent variety for simulation.

---

## Remnant Types

### White Dwarfs

End state of stars <8 M☉ (>97% of all stars). Earth-sized but Sun-mass.

| Property | Typical Value |
|----------|---------------|
| Mass | 0.4–1.4 M☉ |
| Radius | ~Earth-sized |
| Temperature | 4,000–100,000 K (cooling) |
| Progenitor | 0.8–8 M☉ main sequence |

### Neutron Stars / Pulsars

End state of stars 8–25 M☉. City-sized but >1 M☉.

| Property | Typical Value |
|----------|---------------|
| Mass | 1.1–2.5 M☉ |
| Radius | ~10–20 km |
| Spin Period | 1.4 ms – 8.5 s |
| Magnetic Field | 10⁸–10¹⁵ Gauss |

---

## White Dwarf Planetary Systems

### Observational Statistics

| Feature | Occurrence Rate | Notes |
|---------|-----------------|-------|
| Metal pollution | 25–50% | Heavy elements indicating accretion |
| Dusty debris disks | 1–4% | Infrared excess from disrupted material |
| Gaseous debris disks | ~0.06–0.1% | Rare subset with Ca II emission |
| Confirmed giant planets | ~1 known | WD 1856+534 (transiting Jupiter) |

### System Evolution Timeline

```
Main Sequence → Red Giant → AGB → Planetary Nebula → White Dwarf
     │              │         │
     │              │         └── Outer planets survive, orbits expand
     │              └── Inner planets engulfed
     └── Original planetary system
```

**Critical transitions:**

1. **RGB/AGB expansion**: Inner planets (a < few AU) are engulfed
2. **Mass loss**: Star loses 50-80% of mass, surviving planet orbits expand
3. **Cooling phase**: Disk formation from tidally disrupted asteroids/planets

### Orbit Expansion

When the star loses mass, surviving planets move outward:

```rust
fn expanded_orbit(original_a: Au, m_initial: SolarMass, m_final: SolarMass) -> Au {
    Au(original_a.0 * m_initial.0 / m_final.0)
}
// Example: 2 M☉ → 0.6 M☉ white dwarf
// Planet at 10 AU → moves to ~33 AU
```

### Generation Logic

```rust
#[derive(Debug, Clone)]
pub enum WhiteDwarfPlanetOrigin {
    Survivor { original_orbit: Au, expansion_factor: f64 },
    SecondGeneration { disk_mass: EarthMass },
}

impl WhiteDwarfSystem {
    pub fn generate(wd: &WhiteDwarf, rng: &mut impl Rng) -> Option<Self> {
        if !rng.gen_bool(0.35) { return None; }
        
        let debris = Self::generate_debris(wd, rng);
        let surviving_planets = Self::generate_survivors(wd, rng);
        
        Some(Self { debris, planets: surviving_planets })
    }
    
    fn generate_debris(wd: &WhiteDwarf, rng: &mut impl Rng) -> DebrisSignature {
        let cooling_age_gyr = wd.cooling_age.as_gyr();
        let has_dusty_disk = rng.gen_bool(if cooling_age_gyr < 0.7 { 0.08 } else { 0.02 });
        let has_gaseous_disk = has_dusty_disk && rng.gen_bool(0.03);
        
        DebrisSignature {
            dusty_disk: has_dusty_disk,
            gaseous_disk: has_gaseous_disk,
            metal_pollution: true,
            accretion_rate: Self::sample_accretion_rate(rng),
        }
    }
    
    fn generate_survivors(wd: &WhiteDwarf, rng: &mut impl Rng) -> Vec<Planet> {
        if wd.progenitor_mass.0 < 1.5 || !rng.gen_bool(0.1) { return vec![]; }
        
        let expansion = wd.progenitor_mass.0 / wd.mass.0;
        let n_planets = rng.gen_range(0..=2);
        
        (0..n_planets).map(|_| {
            let original_a = Au(rng.gen_range(5.0..30.0));
            Planet {
                orbit: Au(original_a.0 * expansion),
                mass: EarthMass(rng.gen_range(50.0..500.0)),
                origin: WhiteDwarfPlanetOrigin::Survivor {
                    original_orbit: original_a,
                    expansion_factor: expansion,
                },
            }
        }).collect()
    }
}
```

---

## Neutron Star / Pulsar Planetary Systems

### Observational Statistics

| Feature | Occurrence | Notes |
|---------|------------|-------|
| Planetary systems | ~0.5% | ~1 per 200 pulsars |
| Formation potential | ~1% | Of NS progenitors |
| Known planets | ~6 | All around millisecond pulsars |

### Planet Formation Pathways

```rust
#[derive(Debug, Clone)]
pub enum PulsarPlanetOrigin {
    Fallback { disk_mass: EarthMass },
    DisruptedCompanion { companion_type: CompanionType },
    EvaporatedCompanion,
    Captured,
}
```

### Generation Logic

```rust
impl PulsarSystem {
    pub fn generate(ns: &NeutronStar, rng: &mut impl Rng) -> Option<Self> {
        let occurrence = if ns.is_millisecond { 0.01 } else { 0.001 };
        if !rng.gen_bool(occurrence) { return None; }
        
        let origin = Self::select_origin(ns, rng);
        let planets = Self::generate_planets(&origin, rng);
        Some(Self { planets, origin })
    }
    
    fn select_origin(ns: &NeutronStar, rng: &mut impl Rng) -> PulsarPlanetOrigin {
        if ns.is_millisecond {
            match rng.gen::<f64>() {
                x if x < 0.6 => PulsarPlanetOrigin::DisruptedCompanion {
                    companion_type: CompanionType::RedDwarf,
                },
                x if x < 0.85 => PulsarPlanetOrigin::EvaporatedCompanion,
                _ => PulsarPlanetOrigin::Fallback {
                    disk_mass: EarthMass(rng.gen_range(0.01..1.0)),
                },
            }
        } else {
            match rng.gen::<f64>() {
                x if x < 0.7 => PulsarPlanetOrigin::Fallback {
                    disk_mass: EarthMass(rng.gen_range(0.001..0.1)),
                },
                _ => PulsarPlanetOrigin::Captured,
            }
        }
    }
    
    fn generate_planets(origin: &PulsarPlanetOrigin, rng: &mut impl Rng) -> Vec<PulsarPlanet> {
        match origin {
            PulsarPlanetOrigin::DisruptedCompanion { .. } => {
                let n = rng.gen_range(1..=4);
                (0..n).map(|i| PulsarPlanet {
                    mass: EarthMass(rng.gen_range(0.02..10.0)),
                    orbit: Au(rng.gen_range(0.2..2.0) * (i as f64 + 1.0)),
                    composition: PulsarPlanetComposition::Rocky,
                }).collect()
            }
            PulsarPlanetOrigin::EvaporatedCompanion => {
                vec![PulsarPlanet {
                    mass: EarthMass(rng.gen_range(100.0..400.0)),
                    orbit: Au(rng.gen_range(0.01..0.1)),
                    composition: PulsarPlanetComposition::Diamond,
                }]
            }
            PulsarPlanetOrigin::Fallback { disk_mass } => {
                if disk_mass.0 > 0.05 {
                    vec![PulsarPlanet {
                        mass: EarthMass(rng.gen_range(0.1..5.0)),
                        orbit: Au(rng.gen_range(0.3..1.5)),
                        composition: PulsarPlanetComposition::MetalRich,
                    }]
                } else { vec![] }
            }
            PulsarPlanetOrigin::Captured => {
                vec![PulsarPlanet {
                    mass: EarthMass(rng.gen_range(50.0..500.0)),
                    orbit: Au(rng.gen_range(5.0..50.0)),
                    composition: PulsarPlanetComposition::Rocky,
                }]
            }
        }
    }
}
```

### Exotic Planet Types

```rust
#[derive(Debug, Clone)]
pub enum PulsarPlanetComposition {
    Rocky,
    Diamond,
    MetalRich { iron_fraction: f64, radioactive: bool },
    DegenerateMatter,
}
```

---

## Known Pulsar Planet Systems

| System | Planets | Notes |
|--------|---------|-------|
| PSR B1257+12 | 3 | First exoplanets. 0.02, 4.3, 3.9 M⊕ |
| PSR B1620-26 | 1 | 2.5 Mⱼ circumbinary, ~12.6 Gyr old |
| PSR J1719-1438 | 1 | Diamond planet from evaporated WD |
| PSR J2322-2650 | 1 | Jupiter-mass, JWST detected C₂/C₃ atmosphere |

---

## Weird Things to Generate

All observationally supported:

1. **Diamond planets**: Carbon remnants from evaporated WD companions
2. **Ancient planets**: PSR B1620-26 b is ~12.6 billion years old
3. **Metal-rich second generation**: From supernova debris
4. **Expanded-orbit survivors**: Giants at 30-100 AU
5. **Circumbinary debris**: Around tight WD+brown dwarf binaries
6. **Pollution without visible disk**: Late-stage accretion

---

## References

- Wolszczan & Frail (1992) - PSR B1257+12 planets
- Sigurdsson et al. (2003) - PSR B1620-26 ancient planet
- Bailes et al. (2011) - PSR J1719-1438 diamond planet
- Veras (2016) - Post-main-sequence planetary evolution
- Farihi (2016) - White dwarf debris disks review
