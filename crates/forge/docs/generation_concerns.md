# Planet Generation: Initial Diagnostic Concerns

## Giant Planet Formation Around M-Dwarfs

The current output shows gas giants and ice giants forming readily around low-mass M-dwarfs. For example:

- `0.242 M3V gasGiant 0.00242M☉` — an 800 M⊕ giant around a 0.24 M☉ star (1:100 mass ratio)
- Multiple ice giants around stars in the 0.1–0.2 M☉ range

Empirical occurrence rates suggest this is too generous:

| Stellar Type | Gas Giant Occurrence | Ice Giant Occurrence |
|--------------|---------------------|----------------------|
| G/K stars    | 10–15%              | 15–20%               |
| Early M (M0–M3) | 3–5%             | ~10%                 |
| Late M (M5+) | <1–2%               | Suppressed           |

The core accretion timescale problem is severe for M-dwarf disks: less disk mass, longer orbital periods at equivalent temperatures, and disk dissipation before giants can complete formation.

**Questions to investigate:**
- Is the suppression factor being applied at the correct point in the pipeline?
- Is the factor strong enough?
- Is it applied to gas giants only, or also ice giants?

## System Sparseness

Many systems are empty or contain only distant frozen bodies. This may be physically reasonable (debris disk systems exist), but for a worldbuilding tool, consider whether a tuneable "minimum interesting content" parameter would be useful.

## Terrestrial Planet Rarity

Across both diagnostic outputs, terrestrial-class planets appear infrequently. Only one `terran` was observed in the IMF sample. This could be:

- Correct (small planets are harder to form/retain, especially around M-dwarfs)
- A detection bias in the display (are sub-Earth masses being filtered?)
- A genuine gap in the formation logic

Worth validating against Kepler occurrence rates, which show small planets are actually very common.
