import type {
  MainSequenceStar,
  PlanetarySystem,
  SpectralType,
} from "@magrathea/planetary-wasm";

export function filterBySpectralType(
  systems: PlanetarySystem[],
  spectralType: SpectralType
) {
  return systems.filter(
    ({ stars }) => (stars[0] as MainSequenceStar)?.spectralType === spectralType
  );
}
