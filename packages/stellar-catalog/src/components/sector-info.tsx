import type { PlanetarySystem } from "@magrathea/planetary-wasm";
import { useMemo } from "react";
import { useStars } from "../hooks/use-stars";
import { collectPlanets, filterByClass } from "../utils/planets";
import { DataList } from "./data-list";

type Props = {
  systems: PlanetarySystem[];
};

export function SectorInfo({ systems }: Props) {
  const stars = useStars(systems);

  const planets = useMemo(() => collectPlanets(systems), [systems]);

  const averagePlanets = (planets.length / systems.length).toFixed(3);

  const { compact, transitional, volatile, giant } = useMemo(() => {
    return {
      compact: filterByClass(planets, "Compact"),
      transitional: filterByClass(planets, "Transitional"),
      volatile: filterByClass(planets, "Volatile"),
      giant: filterByClass(planets, "Giant"),
    };
  }, [planets]);

  return (
    <>
      <DataList
        title="System Overview"
        data={[
          ["Total Systems", systems.length, "systems"],
          ["Total Planets", planets.length, "planets"],
          ["Average Planets / System", averagePlanets, "average"],
        ]}
      />

      <DataList
        title="Stellar Distribution"
        data={[
          ["M-Type", stars.M.length, "m-type"],
          ["K-Type", stars.K.length, "k-type"],
          ["G-Type", stars.G.length, "g-type"],
          ["F-Type", stars.F.length, "f-type"],
          ["A-Type", stars.A.length, "a-type"],
        ]}
      />

      <DataList
        title="Planet Classes"
        data={[
          [
            <div key="compact" className="flex gap-1 items-center">
              <div className="bg-slate-300 size-4 rounded-full"></div>Compact
            </div>,
            compact.length,
            "compact",
          ],
          [
            <div key="transitional" className="flex gap-1 items-center">
              <div className="bg-sky-100 size-4 rounded-full"></div>
              Transitional
            </div>,
            transitional.length,
            "transitional",
          ],
          [
            <div key="volatile" className="flex gap-1 items-center">
              <div className="bg-blue-300 size-4 rounded-full"></div>
              Volatile
            </div>,
            volatile.length,
            "volatile",
          ],
          [
            <div key="giant" className="flex gap-1 items-center">
              <div className="bg-orange-300 size-4 rounded-full"></div>
              <div>Giant</div>
            </div>,
            giant.length,
            "giant",
          ],
        ]}
      />
    </>
  );
}
