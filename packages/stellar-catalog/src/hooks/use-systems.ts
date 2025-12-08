import { useEffect, useState } from "react";
import init, {
  generateSystemsBatch,
  type PlanetarySystem,
} from "@magrathea/planetary-wasm";

export function useGenerateSystems(seed: string, count: number) {
  const [systems, setSystems] = useState<PlanetarySystem[]>([]);
  const [initialized, setInitialized] = useState(false);

  useEffect(() => {
    init().then(() => setInitialized(true));
  }, [setInitialized]);

  useEffect(() => {
    if (initialized) {
      Promise.resolve(generateSystemsBatch(seed, count)).then(setSystems);
    }
  }, [seed, count, initialized, setSystems]);

  return {
    systems,
    initialized,
  };
}
