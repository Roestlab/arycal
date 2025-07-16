#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/.."

# 1) Build the crate for wasm:
RUSTFLAGS='--cfg getrandom_backend="wasm_js"' cargo build --release --target wasm32-unknown-unknown

# 2) Run wasm-bindgen:
WASM=target/wasm32-unknown-unknown/release/arycal_gui.wasm
OUT=static/pkg
rm -rf $OUT && mkdir -p $OUT
wasm-bindgen "$WASM" \
  --out-dir $OUT \
  --target web

# 3) Generate index.html if you like:
cat > static/index.html <<'EOF'
<!DOCTYPE html>
<html>
  <head><meta charset="utf-8"/><title>Arycal GUI</title></head>
  <body>
    <div id="root"></div>
    <script type="module">
      import init from './pkg/arycal_gui.js';
      async function run() {
        await init('./pkg/arycal_gui_bg.wasm');
      }
      run();
    </script>
  </body>
</html>
EOF
