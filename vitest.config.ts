import { defineConfig } from "vitest/config";
// import { playwright } from "@vitest/browser-playwright";

export default defineConfig({
  test: {
    browser: {
      enabled: true,
      provider: playwright(),
      instances: [
        {
          browser: "chromium",
        },
      ],
      headless: true,
    },
    // Global test setup
    setupFiles: ["./crates/render/test/setup.ts"],
    // Include test files
    include: [
      "crates/**/test/*.{test,spec}.{js,ts,jsx,tsx}",
      "packages/**/*.{test,spec}.{js,ts,jsx,tsx}",
    ],
    // Global variables for WASM
    globals: true,
    // Timeout for WASM loading
    testTimeout: 10000,
  },
  // Vite configuration for WASM
  server: {
    fs: {
      allow: [".."],
    },
  },
  // Build configuration for tests
  build: {
    target: "esnext",
    rollupOptions: {
      external: ["fs", "path"],
    },
  },
  // Resolve configuration
  resolve: {
    alias: {
      "@": "./src",
      "@/wasm": "./pkg",
    },
  },
  // Define global constants
  define: {
    global: "globalThis",
  },
});
