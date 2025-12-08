import { describe, expect, test } from "vitest";
import { translateX, translateY } from "../src/coords";
import { createMockContext } from "./mock-context";

describe("translateX", () => {
  test("returns a translated center", () => {
    const ctx = createMockContext(128, 128);
    expect(translateX(ctx, 1, 0)).toBe(64);
  });

  test("returns a translated left side", () => {
    const ctx = createMockContext(128, 128);
    expect(translateX(ctx, 1, -64)).toBe(0);
  });

  test("returns a translated right side", () => {
    const ctx = createMockContext(128, 128);
    expect(translateX(ctx, 1, 64)).toBe(128);
  });

  test("returns a translated center", () => {
    const ctx = createMockContext(128, 128);
    expect(translateX(ctx, 1.5, 64)).toBe(160);
  });
});

describe("translateY", () => {
  test("returns a translated center", () => {
    const ctx = createMockContext(128, 128);
    expect(translateY(ctx, 1, 0)).toBe(64);
  });

  test("returns a translated top", () => {
    const ctx = createMockContext(128, 128);
    expect(translateY(ctx, 1, 64)).toBe(0);
  });

  test("returns a translated bottom", () => {
    const ctx = createMockContext(128, 128);
    expect(translateY(ctx, 1, -64)).toBe(128);
  });

  test("returns a translated center", () => {
    const ctx = createMockContext(128, 128);
    expect(translateY(ctx, 1.5, 64)).toBe(-32);
  });
});
