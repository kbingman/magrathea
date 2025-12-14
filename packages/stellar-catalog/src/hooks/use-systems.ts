import { useEffect, useState } from "react";
import init, {
  generateSystemsBatch,
  type MainSequenceStar,
  type PlanetarySystem,
} from "@magrathea/planetary-wasm";
import { filterBySpectralType } from "../utils/stars";

export function useGenerateSystems(seed: string, count: number) {
  const [systems, setSystems] = useState<PlanetarySystem[]>([]);
  const [initialized, setInitialized] = useState(false);

  useEffect(() => {
    init().then(() => setInitialized(true));
  }, [setInitialized]);

  useEffect(() => {
    if (initialized) {
      Promise.resolve(generateSystemsBatch(seed, count)).then((systems) => {
        setSystems(filterBySpectralType(systems, "K"));
      });
    }
  }, [seed, count, initialized, setSystems]);

  return {
    systems,
    initialized,
  };
}
