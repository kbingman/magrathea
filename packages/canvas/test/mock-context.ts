import { vi } from "vitest";

export function createMockContext(
  width: number = 800,
  height: number = 600
): CanvasRenderingContext2D {
  return {
    canvas: {
      width,
      height,
    },
    fillStyle: "",
    strokeStyle: "",
    lineWidth: 1,
    beginPath: vi.fn(),
    arc: vi.fn(),
    fill: vi.fn(),
    stroke: vi.fn(),
    fillRect: vi.fn(),
    strokeRect: vi.fn(),
  } as unknown as CanvasRenderingContext2D;
}
