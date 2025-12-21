# Binary Systems for Magrathea

## Why Binaries Matter

Roughly 50% of Sun-like stars and 25-30% of M-dwarfs are in binary or higher multiple systems. Currently modeling only single stars means modeling the minority case.

Binaries don't just add variety—they add *constrained* variety. The dynamics dictate what's architecturally possible, producing naturally distinct catalog signatures.

## Binary Types and Planetary Constraints

### S-Type (Circumstellar)
Planets orbit one star, companion truncates the outer disk.

- Stability limit roughly 1/3 to 1/2 of binary separation (depends on eccentricity, mass ratio)
- Produces compact systems with a hard outer edge
- Wide binaries (>100 AU) allow nearly normal inner systems
- Medium binaries (1-50 AU) create visibly truncated architectures

### P-Type (Circumbinary)
Planets orbit both stars, inner exclusion zone where orbits are unstable.

- Inner stability boundary roughly 2-4× binary separation (more for eccentric binaries)
- No close-in planets possible—systems start at some minimum distance
- Visually distinct: empty near the star pair, then planets

### Dynamical Effects

- **Kozai-Lidov cycles**: Inclined companions pump eccentricity over long timescales, can destabilize or eject planets
- **Disk truncation during formation**: Binary shapes what *could* form, not just what survives
- **Misaligned disks**: Circumbinary disks can be tilted relative to binary plane

## As a Variety Engine

Binary forcing compresses dynamically interesting timescales. In your n-body experiments, 20-minute runs produced visible ejections and stable survivors—single-star systems need geological time for comparable evolution.

This suggests:
- Single-star systems: Generate statistically, accept as "end states"
- Binary systems: Worth simulating briefly, dynamics produce emergent architecture

Different binary configurations naturally produce different outcomes:
- Close equal-mass → tight circumbinary only
- Wide unequal → nearly independent S-type system
- Eccentric medium → chaotic, survivors have sculpted orbits

## Implementation Path

### Binary Generation
Sample from observed distributions:
- Period: log-normal (Raghavan+ 2010 for FGK stars)
- Mass ratio: roughly uniform, slight twin excess
- Eccentricity: thermal distribution, increases with period

### Stability Zones
Holman & Wiegert 1999 provides empirical fits:
- S-type critical semi-major axis as function of (separation, eccentricity, mass ratio)
- P-type inner boundary similarly parameterized

### Integration with Existing Archetypes
Feed current archetypes the constrained orbital range:
- CompactMulti within allowed S-type zone
- Mixed with truncated outer system
- GiantDominated only if stability zone permits giant placement

### Optional N-Body Validation
For interesting configurations:
- Generate initial placement statistically
- Run short integration (few thousand binary periods)
- Keep survivors with their dynamically-settled elements

## Traveller's Approach

The Scout rules handle companions simply:
- Roll for system nature (solo/binary/trinary) early
- Companion orbit affects "maximum orbits available"
- Close companion adds luminosity, use hotter star's temperature for zones

Crude but produces playable variety. Your physics-based stability calculations would refine this.

## References

- Raghavan+ 2010: Binary statistics for solar-type stars
- Holman & Wiegert 1999: Stability criteria for planetary orbits in binaries
- Quarles+ 2020: Updated circumbinary stability limits
