# Stage 1: Build binaries with musl toolchain
FROM rust:1.85-slim AS builder

WORKDIR /app

# 1. Install the MUSL toolchain, Clang/LLD for C++ support, pkg-config & OpenSSL headers
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

# 2. Tell cc-rs / cargo what compilers to use for the MUSL target
ENV CC=musl-gcc
ENV CXX=clang++
ENV CXXFLAGS="--target=x86_64-linux-musl --sysroot=/usr/x86_64-linux-musl"

# 3. Copy only your workspace manifests first (to get a cached layer of fetched dependencies)
COPY Cargo.toml Cargo.lock ./
COPY crates/arycal/Cargo.toml          crates/arycal/Cargo.toml
COPY crates/arycal-cli/Cargo.toml      crates/arycal-cli/Cargo.toml
COPY crates/arycal-cloudpath/Cargo.toml crates/arycal-cloudpath/Cargo.toml
COPY crates/arycal-common/Cargo.toml   crates/arycal-common/Cargo.toml
COPY crates/arycal-gui/Cargo.toml      crates/arycal-gui/Cargo.toml

# 4. Fetch all dependency crates for the MUSL target (caches on changes to any Cargo.toml/Cargo.lock above)
RUN cargo fetch --target x86_64-unknown-linux-musl

# 5. Now copy the rest of your source code
COPY . .

# 6. Build your two binaries with optimizations & stripping
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

# Copy in the statically-built binaries
COPY --from=builder /app/target/x86_64-unknown-linux-musl/release/arycal    /app/arycal
COPY --from=builder /app/target/x86_64-unknown-linux-musl/release/arycal-gui /app/arycal-gui

ENV PATH="/app:$PATH"
