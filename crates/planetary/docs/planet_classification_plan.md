# Exoplanet Classification System Design

A physics-based taxonomy for procedural planet generation, moving beyond problematic nomenclature like "Super-Earth" and "Ice Giant."

## Overview

This system uses a two-tier classification approach:

1. **PlanetClass** — The physical regime determined by what pressure support mechanism dominates (mass-based)
2. **PlanetType** — The observable expression of that regime (driven by composition, temperature, stellar environment)

This separation allows the same physical class to manifest as very different observable worlds depending on formation history and stellar context.

---

## Primary Classification: Physical Regimes

| Class | Mass Range | Physical Regime | Key Behavior |
|-------|------------|-----------------|--------------|
| **Terrestrial** | M < 2 M⊕ | Self-compression | ρ increases with M, R ∝ M^0.29 |
| **Transitional** | 2–5 M⊕ | Envelope retention threshold | Stellar flux determines rocky vs. gaseous |
| **Volatile** | 5–160 M⊕ | Envelope accumulation | ρ *decreases* with M, R ∝ M^0.56 |
| **Degenerate** | > 160 M⊕ | Electron degeneracy | R nearly constant, ρ ∝ M |

### Key Physics

- **Terrestrial**: Condensate-dominated. Rock and iron under self-gravity. Thin or no atmospheres.
- **Transitional**: The critical regime where photoevaporation determines fate. Above the flux threshold → stripped to rocky core. Below → retains H/He envelope.
- **Volatile**: The density inversion regime. Adding mass means adding low-density envelope, so planets get puffier faster than heavier.
- **Degenerate**: Electron degeneracy pressure supports against gravity. Radius stays ~1 R_Jupiter regardless of mass until you approach brown dwarf territory.

---

## Secondary Classification: Observable Types

### Terrestrial Expressions (M < 2 M⊕)

| Type | Temperature | Characteristics |
|------|-------------|-----------------|
| **SubEarth** | Any | M < 0.5 M⊕, Mars/Mercury class, cannot retain atmosphere long-term |
| **Barren** | Any | Airless, cratered, Mercury/Moon-like |
| **Lava** | > 1500 K | Magma oceans, silicate/metal vapor atmospheres |
| **Desert** | 300–500 K | Thin CO₂ or N₂ atmosphere, Mars/Venus-like |
| **Frozen** | < 200 K | Ice-covered, possible subsurface oceans |
| **Oceanic** | 250–350 K | Surface liquid water, high water fraction |
| **Terran** | 250–320 K | Earth-like, requires specific conditions |
| **Eyeball** | Tidally locked | Terminator habitability, permanent day/night hemispheres |
| **Carbon** | Any | C/O > 0.8, graphite/carbide surface, diamond interior |
| **Iron** | Any | Exposed metallic core, > 70% iron, mantle-stripped |

### Transitional Expressions (2–5 M⊕)

| Type | Condition | Characteristics |
|------|-----------|-----------------|
| **SuperTerran** | High flux or stripped | Massive rocky world, thick secondary atmosphere |
| **MiniNeptune** | Low flux, retained envelope | H/He envelope 1–10% by mass, no solid surface |
| **WaterWorld** | High volatile fraction | Deep water layer, possibly supercritical at depth |
| **Hycean** | H₂ envelope + liquid water | Astrobiologically interesting, ocean under H₂ atmosphere |
| **Chthonian** | Stripped gas giant core | M > 10 M⊕ rocky, remnant of evaporated giant |

### Volatile Expressions (5–160 M⊕)

| Type | Temperature | Characteristics |
|------|-------------|-----------------|
| **SubNeptune** | Warm (close-in) | Puffy, hydrogen-dominated, cloud layers |
| **WarmNeptune** | Intermediate | Neptune-mass at moderate separation, transitional clouds |
| **IceGiant** | Cold (outer system) | Uranus/Neptune-like, methane clouds, complex moons |

### Degenerate Expressions (> 160 M⊕)

| Type | Temperature | Characteristics |
|------|-------------|-----------------|
| **GasGiant** | Cold | Jupiter-like, banded clouds, ring systems, internal heat |
| **HotJupiter** | > 1000 K | Inflated, tidally locked, atmospheric escape |
| **PuffySaturn** | Warm | Inflated sub-Jovian, low density anomaly, 0.1–0.5 M_J |

---

## Stellar Context

The star determines which types are possible at a given orbital distance.

### Key Stellar Parameters

```rust
pub struct StellarContext {
    pub luminosity: f64,      // L☉
    pub temperature: f64,     // K (determines UV flux)
    pub age: f64,             // Gyr (time for photoevaporation)
}
```

### Derived Quantities

| Quantity | Formula | Effect |
|----------|---------|--------|
| **Incident Flux** | F = L / a² | Drives temperature, evaporation |
| **Equilibrium Temperature** | T_eq = 278 × F^0.25 | Surface/atmosphere temperature |
| **UV Multiplier** | (T_eff / 5778)^4 | Photoevaporation rate |
| **Evaporation Threshold** | F_thresh ∝ M^2.5 | Determines envelope retention |

### Habitable Zone

- **Inner edge**: a = √(L / 1.1) AU — runaway greenhouse limit
- **Outer edge**: a = √(L / 0.53) AU — maximum greenhouse limit

---

## Atmospheric Classification

### Retention Categories

| Category | Mass Fraction | Examples |
|----------|---------------|----------|
| **None** | 0 | Airless moons, Mercury |
| **Thin** | 10⁻⁷ – 10⁻⁴ | Mars, Earth |
| **Thick** | 10⁻⁴ – 0.1 | Venus, mini-Neptunes |
| **Massive** | > 0.1 | Gas giants |
| **Stripped** | ~10⁻⁵ | Former mini-Neptunes |

### Composition Categories

| Category | Dominant Species | Temperature Range |
|----------|------------------|-------------------|
| **None** | — | Any (too small to retain) |
| **Reducing** | H₂, CH₄, NH₃ | Cold or massive |
| **Oxidizing** | CO₂, N₂, O₂ | Moderate |
| **Water Vapor** | H₂O | Runaway greenhouse |
| **Hydrogen** | H₂/He | Mini-Neptunes, giants |
| **Exotic** | SiO, Fe, Na, K | > 2000 K lava worlds |

---

## Water Worlds: Phase Structure

Water-rich planets have complex interior structure depending on pressure and temperature.

### Water Phases

| Phase | Pressure | Temperature | Notes |
|-------|----------|-------------|-------|
| Ice Ih | < 2 kbar | < 273 K | Normal ice |
| Ice III–VI | 2–11 kbar | Various | High-pressure ices |
| Ice VII | > 11 kbar | High | Very high pressure |
| Liquid | < 221 bar | 273–647 K | Normal ocean |
| Supercritical | > 221 bar | > 647 K | No liquid/gas distinction |

### Surface States

| State | Conditions | Description |
|-------|------------|-------------|
| **Frozen** | T_eq < 273 K | Ice shell, possible subsurface ocean |
| **Liquid Ocean** | 273–373 K, moderate P | Surface ocean, possibly over high-P ice |
| **Steam** | > 373 K | Dense steam atmosphere |
| **Supercritical** | > 647 K, > 221 bar | No distinct surface |

---

## Surface Hydrology

For terrestrial worlds with partial or complete liquid coverage (as opposed to water worlds with deep global oceans).

### Rust Types

```rust
struct SurfaceHydrology {
    ocean_fraction: f64,        // 0.0–1.0 (Earth ~0.71)
    mean_depth_km: f64,         // Average ocean depth (Earth ~3.7 km)
    max_depth_km: f64,          // Deepest point (Earth ~11 km)
    ice_cap_fraction: f64,      // Polar/glacial ice coverage
    salinity_ppt: f64,          // Parts per thousand (Earth ~35)
    ocean_type: OceanType,
    distribution: OceanDistribution,
}

enum OceanType {
    Water,                      // H₂O, 273–373 K at 1 bar
    Ammonia,                    // NH₃-H₂O mix, cold outer system (~195 K)
    Hydrocarbon,                // CH₄/C₂H₆, Titan-style (~90–180 K)
    Supercritical,              // Hot, high-pressure, no surface
    Brine,                      // High-salinity, lower freezing point
}

enum OceanDistribution {
    Global,                     // No exposed land
    Hemispheric,                // Ocean/continent dichotomy
    Scattered,                  // Multiple ocean basins (Earth-like)
    Polar,                      // Liquid only at poles (warm world)
    Equatorial,                 // Liquid only at equator (cold world)
    Subsurface,                 // Beneath ice shell
    Terminator,                 // Tidally locked, ring of liquid at day/night boundary
}
```

### Ocean Depth Limits

Maximum ocean depth before high-pressure ice forms at the seafloor:

| Ocean Type | Max Depth | Reason |
|------------|-----------|--------|
| Water | ~100–200 km | Ice VII forms at ~2 GPa |
| Ammonia | ~300 km | Higher pressure tolerance |
| Hydrocarbon | ~200 km | Depends on composition |

Beyond these depths, the ocean floor is high-pressure ice rather than rock — cutting off silicate-water interaction important for habitability.

### Integration with World Types

| World Type | SurfaceHydrology | Notes |
|------------|------------------|-------|
| **Terran** | Required | Earth-like, ocean_fraction 0.3–0.9 |
| **Oceanic** | Required | ocean_fraction > 0.95, deep oceans |
| **Desert** | Optional | Remnant seas, subsurface aquifers |
| **Frozen** | Subsurface | Ice shell over liquid layer |
| **Eyeball** | Terminator | Liquid ring at day/night boundary |
| **Hycean** | Required | Ocean under H₂ atmosphere |

### Habitability Considerations

Surface liquid water requires:

```
273 K < T_surface < 373 K  (at 1 bar)
```

But pressure extends the range:
- Higher pressure → liquid stable to higher T (up to 647 K critical point)
- Dissolved salts → liquid stable to lower T (brines down to ~250 K)

Key factors for ocean stability:
- Stellar flux (not too hot/cold)
- Atmospheric pressure (prevents boiling)
- Magnetic field (prevents atmospheric stripping)
- Carbonate-silicate cycle (long-term climate regulation)

---

## Subsurface Oceans

Low-mass icy bodies can maintain liquid water beneath ice shells regardless of surface conditions. This applies across the full mass range — from Enceladus (0.00002 M⊕) to Mars-sized bodies.

### Key Insight

Surface oceans require atmospheric pressure to prevent boiling/sublimation, which requires M > ~0.1 M⊕.

Subsurface oceans have no minimum mass — the ice shell acts as a pressure vessel. You only need sufficient heat to maintain a liquid layer.

### Examples

| Body | Mass (M⊕) | Ice Shell | Ocean Depth | Heat Source |
|------|-----------|-----------|-------------|-------------|
| Enceladus | 0.00002 | ~20–25 km | ~10 km | Tidal |
| Europa | 0.008 | ~10–30 km | ~60–150 km | Tidal |
| Ganymede | 0.025 | ~150 km | ~100 km | Radiogenic |
| Titan | 0.023 | ~50–100 km | ~50–100 km | Mixed |
| Pluto | 0.002 | ~150–200 km | ~100 km? | Radiogenic |
| Callisto | 0.018 | >100 km | Possible | Minimal |

### Rust Types

```rust
struct SubsurfaceOcean {
    ice_shell_km: f64,              // Thickness of ice shell
    ocean_depth_km: f64,            // Liquid layer thickness
    ocean_type: SubsurfaceOceanType,
    heat_source: HeatSource,
    silicate_contact: bool,         // Ocean touches rock? (habitability key)
    cryovolcanism: bool,            // Active venting through ice?
}

enum SubsurfaceOceanType {
    Water,                          // Pure H₂O
    Brine {                         // Salts lower freezing point
        salinity_ppt: f64,
        dominant_salt: Salt,
    },
    AmmoniaWater {                  // NH₃-H₂O eutectic (~176 K)
        ammonia_fraction: f64,
    },
}

enum Salt {
    NaCl,           // Sodium chloride
    MgSO4,          // Magnesium sulfate (Europa?)
    NH4Cl,          // Ammonium chloride
}

enum HeatSource {
    Tidal {
        power_tw: f64,
        eccentricity: f64,
        resonance: Option<String>,  // e.g., "Laplace" for Io-Europa-Ganymede
    },
    Radiogenic {
        power_tw: f64,
    },
    Combined {
        tidal_tw: f64,
        radiogenic_tw: f64,
    },
    Residual,       // Primordial heat only, slowly freezing
}
```

### Ice Shell Dynamics

The ice shell thickness reflects thermal equilibrium:

```
Q_heat = Q_conduction + Q_convection
```

Thin shells (~10 km) indicate high heat flux — convection dominates, possible ice tectonics.

Thick shells (>100 km) indicate low heat flux — conductive lid, geologically dead surface.

### Pluto: A Special Case

Pluto demonstrates that subsurface oceans can exist far from any tidal heat source:

- No significant tidal heating (Charon is tidally locked)
- Radiogenic heating from rocky core (~3–4 TW)
- Ammonia antifreeze likely (lowers freezing point to ~176 K)
- Evidence: Sputnik Planitia's location suggests positive mass anomaly (dense subsurface liquid)

This implies subsurface oceans may be common even on isolated Kuiper Belt objects with sufficient rock fraction.

### Habitability Implications

The key variable is `silicate_contact`:

| Configuration | Habitability | Example |
|---------------|--------------|---------|
| Ocean → Rock | High potential | Europa, Enceladus |
| Ocean → High-P Ice → Rock | Lower potential | Ganymede, Titan |
| Ocean → Ice (no rock) | Minimal | Pure ice bodies |

Silicate-water interaction provides:
- Chemical energy (serpentinization)
- Nutrient cycling
- Redox gradients

Enceladus is particularly interesting: small enough that its ocean likely contacts rock directly, with confirmed hydrothermal activity.

### Integration with World Types

`SubsurfaceOcean` can attach to:

| World Type | Typical Configuration |
|------------|----------------------|
| Frozen | Ice shell over liquid, any mass |
| SubEarth | Small icy bodies (Europa-class) |
| Eyeball | Subsurface ocean beneath frozen nightside |

This replaces the simple `Subsurface` variant in `OceanDistribution` with a full structural model.

---

## Exotic World Types

### Lava Worlds (T > 1500 K)

- Permanent magma oceans (tidally locked: dayside only)
- Silicate vapor atmospheres (SiO, Na, K)
- Metal vapor at T > 2500 K (Fe, Mg)
- Example: CoRoT-7b

### Carbon Worlds

- Form around carbon-rich stars (C/O > 0.8)
- Graphite crust, silicon carbide mantle
- Diamond layers at depth under high pressure
- Tar-like hydrocarbon deposits possible

### Iron Worlds

- Exposed metallic cores (stripped mantles via giant impacts)
- Mercury-like but more extreme
- May have metal vapor atmospheres if hot enough

---

## Mass-Radius Relations

Each class has its own power law with intrinsic scatter (from Wolfgang et al. 2016):

```
M/M⊕ ~ Normal(μ = C × (R/R⊕)^γ, σ = σ_M)
```

| Class | Coefficient (C) | Exponent (γ) | Scatter (σ_M) |
|-------|-----------------|--------------|---------------|
| Terrestrial | 0.55 | 0.29 | 0.5 M⊕ |
| Transitional | 2.7 | 1.3 | 1.9 M⊕ |
| Volatile | 7×10⁻⁸ | 0.56 | 1.9 M⊕ |
| Degenerate | 4×10⁸ | -0.024 | 0.3 M⊕ |

Note: Coefficients assume SI units (kg, m). Adjust for Earth units in implementation.

---

## Generation Flow

```
1. Generate mass (from formation model)
         ↓
2. Classify by mass → PlanetClass
         ↓
3. Sample radius from M-R relation (with scatter)
         ↓
4. Compute bulk density
         ↓
5. Generate composition fractions
         ↓
6. Get stellar context (L, T_eff, age)
         ↓
7. Compute incident flux, T_eq, UV flux
         ↓
8. Apply photoevaporation (strip envelope if F > F_thresh)
         ↓
9. Determine PlanetType from class + temperature + composition
         ↓
10. Generate atmospheric profile
         ↓
11. Generate surface/interior details specific to type
```

---

## Rust Type Hierarchy

```rust
// Primary classification
enum PlanetClass {
    Terrestrial,
    Transitional,
    Volatile,
    Degenerate,
}

// Secondary classification with type-specific data
enum PlanetType {
    // Terrestrial (M < 2 M⊕)
    SubEarth(SubEarthWorld),
    Barren(BarrenWorld),
    Lava(LavaWorld),
    Desert(DesertWorld),
    Oceanic(OceanicWorld),
    Frozen(FrozenWorld),
    Terran(TerranWorld),
    Eyeball(EyeballWorld),
    Carbon(CarbonWorld),
    Iron(IronWorld),
    
    // Transitional (2–5 M⊕)
    SuperTerran(SuperTerranWorld),
    MiniNeptune(MiniNeptuneWorld),
    WaterWorld(WaterWorldData),
    Hycean(HyceanWorld),
    Chthonian(ChthonianWorld),
    
    // Volatile (5–160 M⊕)
    SubNeptune(SubNeptuneWorld),
    WarmNeptune(WarmNeptuneWorld),
    IceGiant(IceGiantWorld),
    
    // Degenerate (> 160 M⊕)
    GasGiant(GasGiantWorld),
    HotJupiter(HotJupiterWorld),
    PuffySaturn(PuffySaturnWorld),
}

// Unified planet struct
struct Planet {
    // Physical
    mass: Mass,
    radius: Length,
    density: Density,
    semi_major_axis: Length,
    
    // Classification
    class: PlanetClass,
    planet_type: PlanetType,
    
    // Composition
    composition: Composition,
    
    // Derived
    surface_gravity: Acceleration,
    escape_velocity: Velocity,
    equilibrium_temp: Temperature,
}
```

---

## Moon Classification

The same physical regimes apply to moons, but with different energy sources and additional constraints.

### Energy Sources

Unlike planets, moons may derive significant heating from tidal forces rather than stellar flux:

```rust
enum MoonEnergySource {
    Tidal {
        eccentricity: f64,
        orbital_period: Duration,
        parent_mass: Mass,
    },
    Radiogenic,
    Stellar,
    Combined {
        tidal_fraction: f64,
        stellar_fraction: f64,
        radiogenic_fraction: f64,
    },
}
```

### Tidal Heating

Tidal dissipation power scales as:

```
P_tidal ∝ (M_parent² × R_moon⁵ × e²) / (a⁶ × Q)
```

Where Q is the tidal dissipation factor (~100 for rocky, ~10⁴ for icy bodies).

### Moon-Specific Constraints

| Constraint | Effect |
|------------|--------|
| **Hill Sphere** | Parent's Hill radius caps maximum moon mass |
| **Roche Limit** | Destroys large moons that orbit too close |
| **Atmospheric Retention** | No H/He retention (too small + parent's Hill sphere) |
| **Tidal Locking** | Common for close-in moons |
| **Formation Mode** | Capture vs accretion affects composition |

### Moon Types

Moons use the same PlanetType variants where applicable, but some types are moon-specific:

| Type | Energy Source | Characteristics |
|------|---------------|-----------------|
| **Volcanic (Io-like)** | Tidal | Silicate volcanism, SO₂ atmosphere |
| **Cryovolcanic (Europa-like)** | Tidal | Ice shell, subsurface ocean, water geysers |
| **Frozen (Callisto-like)** | Radiogenic only | Ancient surface, minimal activity |
| **Captured (Triton-like)** | Mixed | Retrograde orbit, exotic composition |
| **Titan-like** | Radiogenic + Solar | Dense N₂ atmosphere, hydrocarbon lakes |

### Subsurface Ocean Conditions

For icy moons, liquid water requires:

```
T_interior > 273 K
```

Heat balance: `Q_tidal + Q_radiogenic = Q_conductive_loss`

Ice shell thickness depends on heat flux:
- High tidal heating → thin ice (Europa: ~10–30 km)
- Low heating → thick ice or frozen solid (Callisto: >100 km)

---

## References

- **Durand-Manterola (2011)** — Three-class system based on pressure support physics
- **Wolfgang & Lopez (2015)** — Probabilistic M-R relation with intrinsic scatter
- **Wolfgang, Rogers & Ford (2016)** — Bayesian M-R relation framework
- **Lopez & Fortney (2013, 2014)** — Photoevaporation thresholds, thermal evolution
- **Fulton et al. (2017)** — Radius valley / Fulton gap observational confirmation

---

## Future Extensions

- [ ] Biosignature potential assessment
- [x] Moon system classification
- [ ] Moon system generation (orbital mechanics, resonances)
- [ ] Ring system generation  
- [ ] Magnetic field modeling
- [ ] Tectonic state from internal heat
- [ ] Detailed cloud layer modeling
- [ ] Spectral signatures for each type
- [ ] UltraHotJupiter type (T > 2000 K, molecular dissociation, metal weather)
- [ ] Rogue planets (StellarContext::None, radiogenic heat only)
- [ ] Binary star support (S-type, P-type orbits)
