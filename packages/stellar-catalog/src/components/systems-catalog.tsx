import type { PlanetarySystem } from "@magrathea/planetary-wasm";
import { SystemCard } from "./system-card";

type Props = {
  systems: PlanetarySystem[];
  maxWidth: number;
};

export function SystemCatalog({ systems, maxWidth }: Props) {
  return (
    <div className="grid gap-2 fixed left-4 w-[1060px] h-[calc(100vh-80px)] overflow-scroll">
      {systems
        // .filter(({ stars }) => stars[0]?.spectralType === "G")
        .map(({ metadata, planets, stars }) => (
          <SystemCard
            key={`system-${metadata.id}`}
            id={metadata.id}
            stars={stars}
            planets={planets}
            name={metadata.catalogName || ""}
            maxWidth={maxWidth}
          />
        ))}
    </div>
  );
}
