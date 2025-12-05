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

test-rust:
  cargo test --quiet
