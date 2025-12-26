# Moon System Generation Plan

## Overview

This document outlines the design for adding moon (satellite) generation to the statistical planetary system generator. Moons are treated as planets orbiting planets - the same `PlanetClass`/`PlanetType` classification system applies.

## Design Philosophy

**Moons are planets too**: Under the geophysical definition, a moon is simply a planet that orbits another planet rather than a star. We reuse the existing classification system:
- A moon has a `PlanetClass` (Compact, Transitional, etc.)
- A moon has a `PlanetType` (Frozen, Oceanic, Lava, etc.)
- Tidal heating from the host planet affects moon classification

**Statistical sampling, not simulation**: Like the planet generator, we sample from occurrence rates rather than simulating formation physics. The formation mechanism (impact, capture, co-accretion) is recorded but doesn't drive a physics simulation.

## Formation Mechanisms

Three primary mechanisms create moons:

### 1. Giant Impact (Earth-Moon style)
- Large impactor strikes proto-planet
- Debris disk coalesces into one large moon
- Typical for terrestrial planets
- Produces single, relatively large moon (1-5% host mass)

### 2. Gravitational Capture (Triton style)
- Passing body captured by planet's gravity
- Often results in irregular/retrograde orbit
- Common for ice/gas giants
- Produces smaller, more distant moons

### 3. Co-accretion (Galilean style)
- Moons form in circumplanetary disk alongside planet
- Produces multiple moons in resonant chains
- Dominant for gas giants
- Regular prograde orbits, often in MMR

## Planet-Type Formation Rules

### Terrestrial Planets (< 2 M⊕)
```
Moon count: [70% none, 25% one, 5% two, 0% three+]
Mechanisms: 60% impact, 30% capture, 10% co-accretion
Mass ratio: 0.1% - 5% of host mass
Orbit range: 0.8% - 2.5% of Hill sphere
```
Example: Earth-Moon system

### Super-Earths (2-10 M⊕)
```
Moon count: [60% none, 28% one, 10% two, 2% three+]
Mechanisms: 30% impact, 15% capture, 25% co-accretion
Mass ratio: 0.5% - 2% of host mass
Orbit range: 0.6% - 2.0% of Hill sphere
```
Note: Reduced formation rate due to higher escape velocity and potential for envelope retention.

### Ice Giants (10-50 M⊕)
```
Moon count: [15% none, 25% one, 35% two, 25% three+]
Mechanisms: 10% impact, 40% capture, 50% co-accretion
Mass ratio: 0.01% - 0.5% of host mass
Orbit range: 0.4% - 1.5% of Hill sphere
```
Examples: Uranus (27 moons), Neptune (14 moons, including captured Triton)

### Gas Giants (> 50 M⊕)
```
Moon count: [3% none, 12% one, 35% two, 50% three+]
Mechanisms: 5% impact, 35% capture, 60% co-accretion
Mass ratio: 0.001% - 0.1% of host mass
Orbit range: 0.2% - 1.2% of Hill sphere
```
Examples: Jupiter (95 moons), Saturn (146 moons)

## Physical Constraints

### Hill Sphere
The Hill sphere defines the region where the planet's gravity dominates:

```
R_Hill = a × (M_planet / (3 × M_star))^(1/3)
```

Moons must orbit within ~30% of Hill sphere for long-term stability.

### Roche Limit
Inside the Roche limit, tidal forces disrupt moons into rings:

```
R_Roche ≈ 2.44 × R_planet × (ρ_planet / ρ_moon)^(1/3)
```

### Tidal Locking
Close-in moons become tidally locked:

```
τ_lock ∝ a^6 / (M_planet × R_moon^2)
```

Most major moons are tidally locked to their host planet.

## Tidal Heating

Critical for moon classification. Io is the most volcanically active body in the Solar System due to tidal heating from Jupiter.

### Heat flux calculation
```
Q_tidal ∝ (M_planet × R_moon^5 × e^2) / (a^6 × Q)
```

Where Q is the tidal dissipation factor.

### Tidal heating effects by distance:
| Orbital Distance | Heating Level | Example |
|------------------|---------------|---------|
| < 5 R_planet | Extreme (Io-like) | Lava world, continuous volcanism |
| 5-15 R_planet | Moderate (Europa-like) | Subsurface ocean, cryovolcanism |
| 15-40 R_planet | Mild (Ganymede-like) | Possible subsurface ocean |
| > 40 R_planet | Negligible | Frozen, Callisto-like |

### Resonance amplification
Moons in mean-motion resonance (like Io-Europa-Ganymede) have forced eccentricities that maintain tidal heating over billions of years.

## Ring Systems

Gas giants may have ring systems. Rings form from:
- Disrupted moons (inside Roche limit)
- Captured debris
- Cometary impacts

### Ring probability by planet type:
| Planet Type | Ring Probability | Ring Extent |
|-------------|------------------|-------------|
| Terrestrial | ~0% | - |
| Super-Earth | ~1% | Minimal |
| Ice Giant | ~50% | Modest (Uranus-like) |
| Gas Giant | ~80% | Extensive (Saturn-like) |

Ring presence should be tracked but detailed structure is out of scope.

## Data Structures

### Moon (Satellite)
```rust
pub struct Moon {
    /// Unique identifier within system
    pub id: String,
    /// Display name (e.g., "KV-4729 b I" for first moon of planet b)
    pub name: String,
    /// Moon mass
    pub mass: Mass,
    /// Moon radius
    pub radius: Length,
    /// Orbital semi-major axis (from host planet center)
    pub semi_major_axis: Length,
    /// Orbital eccentricity
    pub eccentricity: f64,
    /// Physical classification (same as planets)
    pub class: PlanetClass,
    /// Observable expression
    pub moon_type: PlanetType,
    /// Formation mechanism
    pub formation: MoonFormation,
    /// Equilibrium temperature (includes tidal + stellar)
    pub surface_temp: f64,
    /// Tidal heating flux (W/m²)
    pub tidal_heat_flux: f64,
    /// Is tidally locked to host?
    pub tidally_locked: bool,
}
```

### MoonFormation (enum)
```rust
pub enum MoonFormation {
    GiantImpact,
    GravitationalCapture { retrograde: bool },
    CoAccretion { resonance_order: Option<(u8, u8)> },
}
```

### MoonSystem (attached to Planet)
```rust
pub struct MoonSystem {
    /// Major moons (large enough for classification, ~100+ km)
    pub moons: Vec<Moon>,
    /// Has ring system?
    pub has_rings: bool,
}
```

## Naming Convention

Following IAU-style convention for exomoons:
- Planet: `KV-4729 b`
- First moon: `KV-4729 b I` (Roman numerals, ordered by discovery/distance)
- Second moon: `KV-4729 b II`

For generated systems, order by semi-major axis (innermost = I).

## Integration with Planet

Add `moon_system: Option<MoonSystem>` field to `Planet`. Moons are part of the planetary system and belong with their host planet - there's no use case for separating them.

## Generation Pipeline

```
1. After planet generation, for each planet:
   a. Check if planet qualifies (mass > 0.05 M⊕)
   b. Determine planet type (Terrestrial/SuperEarth/IceGiant/GasGiant)
   c. Sample moon count from type-specific distribution
   d. For each moon:
      - Sample formation mechanism
      - Sample mass ratio from type-specific range
      - Sample orbital distance (fraction of Hill sphere)
      - Apply tidal heating based on distance
      - Classify moon using PlanetType::from_environment()
   e. Check for ring system (gas/ice giants)
   f. Sort moons by semi-major axis, assign names

2. Return planet with attached MoonSystem
```

## Implementation Phases

### Phase 1: Core Structure
- [ ] Add `Moon` struct to planetary crate
- [ ] Add `MoonFormation` enum
- [ ] Add `MoonSystem` struct
- [ ] Add `moon_system` field to `Planet`
- [ ] Add moon naming helper (Roman numerals)

### Phase 2: Generation Logic
- [ ] Add planet-type-specific formation rules
- [ ] Implement moon count sampling
- [ ] Implement moon mass/orbit sampling
- [ ] Integrate with Hill sphere calculation
- [ ] Add tidal heating calculation

### Phase 3: Classification
- [ ] Adapt `PlanetType::from_environment()` for moons
- [ ] Add tidal heating to temperature calculation
- [ ] Handle tidally-locked classification (Eyeball moons)
- [ ] Add Io-like lava moon classification

### Phase 4: Ring Systems
- [ ] Add ring probability by planet type
- [ ] Track ring presence in `MoonSystem`
- [ ] Update Giant/IceGiant variants to use MoonSystem rings

### Phase 5: Validation
- [ ] Test moon count distributions match Solar System
- [ ] Test tidal heating produces Io-like moons for close-in orbits
- [ ] Test Galilean-style systems for gas giants

## References

- Canup (2004) - "Simulations of a late lunar-forming impact"
- Agnor & Hamilton (2006) - "Neptune's capture of its moon Triton"
- Sasaki et al. (2010) - "Origin of the regular satellites of giant planets"
- Heller & Barnes (2013) - "Exomoon habitability constrained by illumination and tidal heating"
- Dobos et al. (2017) - "Tidal heating and the habitability of the TRAPPIST-1 exoplanets"

## Design Decisions

1. **Irregular/minor moons**: Don't track at all. Only major moons (large enough to be spherical, ~100+ km) are generated and classified. The swarm of tiny captured rocks isn't useful for world-building.

## Open Questions

1. How to handle binary moons (Pluto-Charon)?
   - Rare edge case, probably out of scope for v1
