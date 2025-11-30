import tailwindcss from "@tailwindcss/vite";
import { defineConfig } from "vite";

export default defineConfig({
  build: {
    sourcemap: "hidden",
  },
  plugins: [tailwindcss()],
  server: {
    port: 3003,
    strictPort: true,
  },
});
