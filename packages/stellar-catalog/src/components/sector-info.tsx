import type { PlanetarySystem } from "@magrathea/planetary-wasm";
import { useMemo } from "react";
import { collectPlanets, filterByClass } from "../utils/planets";
import { DataList } from "./data-list";
import { useStars } from "../hooks/use-stars";

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
          ["Total Systems", systems.length],
          ["Total Planets", planets.length],
          ["Average Planets / System", averagePlanets],
        ]}
      />

      <DataList
        title="Stellar Distribution"
        data={[
          ["M-Type", stars.M.length],
          ["K-Type", stars.K.length],
          ["G-Type", stars.G.length],
          ["F-Type", stars.F.length],
          ["A-Type", stars.A.length],
        ]}
      />

      <DataList
        title="Planet Classes"
        data={[
          [
            <div className="flex gap-1 items-center">
              <div className="bg-slate-300 size-4 rounded-full"></div>Compact
            </div>,
            compact.length,
          ],
          [
            <div className="flex gap-1 items-center">
              <div className="bg-sky-100 size-4 rounded-full"></div>
              Transitional
            </div>,
            transitional.length,
          ],
          [
            <div className="flex gap-1 items-center">
              <div className="bg-blue-300 size-4 rounded-full"></div>
              Volatile
            </div>,
            volatile.length,
          ],
          [
            <div className="flex gap-1 items-center">
              <div className="bg-orange-300 size-4 rounded-full"></div>
              <div>Giant</div>
            </div>,
            giant.length,
          ],
        ]}
      />
    </>
  );
}
