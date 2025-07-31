# Stage 1: Build both binaries with the MUSL toolchain
FROM clux/muslrust:1.85.1-stable AS builder
WORKDIR /app

# 1) Install only what we need for a static build
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
      pkg-config \
      libssl-dev \
      musl-tools && \
    rm -rf /var/lib/apt/lists/*

# 2) Tell cargo/cc-rs to use the musl wrappers
ENV CC=musl-gcc
ENV CXX=musl-g++

# 3) Copy your workspace manifest & crates folder and pre-fetch deps
COPY Cargo.toml ./
COPY crates/ ./crates/
RUN cargo fetch --target x86_64-unknown-linux-musl

# 4) Copy the rest of your source and build
COPY . .
ENV RUSTFLAGS="-C strip=symbols -C lto=yes"
RUN cargo build --release --target x86_64-unknown-linux-musl --bin arycal
RUN cargo build --release --target x86_64-unknown-linux-musl --bin arycal-gui

# Stage 2: Minimal runtime
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
