# syntax=docker/dockerfile:1.4

# ─── Stage 1: build both binaries ────────────────────────────────────────
FROM rust:1.85-slim AS builder
WORKDIR /app

# 1) Install any system deps your crate needs (MPI, GTK, etc.) plus a C++ toolchain
RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
      build-essential \
      pkg-config \
      libssl-dev \
      libopenmpi-dev openmpi-bin \
      libx11-dev libxrandr-dev libxi-dev \
      libgl1-mesa-dev libegl1-mesa libgtk-3-dev && \
    rm -rf /var/lib/apt/lists/*

# 2) Copy just the Cargo manifests & crates folder so `cargo fetch` can cache
COPY Cargo.toml .
COPY crates/ ./crates/

# 3) Fetch all Rust dependencies (glibc target)
RUN cargo fetch

# 4) Copy the rest of your source code
COPY . .

# 5) Build optimized, stripped release binaries into target/
ENV RUSTFLAGS="-C lto=thin -C strip=symbols"

# mount /app/target as a cache so it lives outside the build container
RUN --mount=type=cache,target=/app/target \
    cargo build --release --bin arycal

RUN --mount=type=cache,target=/app/target \
    cargo build --release --bin arycal-gui

# ─── Stage 2: runtime ────────────────────────────────────────────────────
FROM debian:bullseye-slim
RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
      ca-certificates && \
    update-ca-certificates && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /app
COPY --from=builder /app/target/release/arycal    /usr/local/bin/arycal
COPY --from=builder /app/target/release/arycal-gui /usr/local/bin/arycal-gui

ENTRYPOINT ["arycal"]
CMD ["--help"]
