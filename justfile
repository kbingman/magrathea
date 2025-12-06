# List available commands
default:
  @just --list

# Install all dependencies
install: install-rust install-ts

install-rust:
  cargo fetch

install-ts:
  npm install

build-wasm:
  wasm-pack build crates/magrathea-wasm --target web --scope magrathea
  wasm-pack build crates/planetary-wasm --target web --scope magrathea
  # Patch TypeScript types to include 'name' and 'description' fields in PlanetType variants
  sed -i '' 's/{ type: "\([^"]*\)"/{ type: "\1"; name: string; description: string/g' crates/planetary-wasm/pkg/planetary_wasm.d.ts

test-rust:
  cargo test --quiet
