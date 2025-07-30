# Stage 1: Build binaries with musl toolchain
FROM rust:1.85-slim AS builder

WORKDIR /app

# 1. Install musl headers & tools, clang for C++, pkg-config & OpenSSL dev
# 2. Add the x86_64‐unknown‐linux‐musl Rust target
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
      build-essential \
      musl-tools \
      musl-dev \
      clang \
      lld \
      pkg-config \
      libssl-dev && \
    rustup target add x86_64-unknown-linux-musl && \
    rm -rf /var/lib/apt/lists/*

# 3. Tell cargo/cc-rs to use musl-gcc for C, clang++ for C++ (with the musl sysroot)
ENV CC=musl-gcc
ENV CXX=clang++
ENV CXXFLAGS="--target=x86_64-linux-musl --sysroot=/usr/x86_64-linux-musl"

# Copy your Cargo files first (for caching)
COPY Cargo.toml Cargo.lock ./
COPY arycal/Cargo.toml arycal/
COPY arycal-gui/Cargo.toml arycal-gui/

# Fetch deps
RUN mkdir arycal/src && \
    echo "fn main() {}" > arycal/src/main.rs && \
    mkdir arycal-gui/src && \
    echo "fn main() {}" > arycal-gui/src/main.rs && \
    cargo build --release --target x86_64-unknown-linux-musl && \
    rm -rf target/x86_64-unknown-linux-musl/release/deps/*

# Now copy your real source
COPY . .

# Build optimized & stripped binaries
ENV RUSTFLAGS="-C strip=symbols -C lto=yes"
RUN cargo build --release --target x86_64-unknown-linux-musl --bin arycal
RUN cargo build --release --target x86_64-unknown-linux-musl --bin arycal-gui

# Stage 2: Minimal runtime image
FROM debian:bullseye-slim

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
      procps \
      ca-certificates && \
    update-ca-certificates && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY --from=builder /app/target/x86_64-unknown-linux-musl/release/arycal    /app/arycal
COPY --from=builder /app/target/x86_64-unknown-linux-musl/release/arycal-gui /app/arycal-gui

ENV PATH="/app:$PATH"
