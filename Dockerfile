# Stage 1: Build both binaries with the MUSL toolchain
FROM clux/muslrust:1.85.1-stable AS builder

WORKDIR /app

# Install only pkg-config & OpenSSL headers (musl-g++ is already bundled)
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
      pkg-config \
      libssl-dev && \
    rm -rf /var/lib/apt/lists/*

# **Ensure cargo’s cc-rs finds a C++ compiler**
ENV CC=musl-gcc
ENV CXX=musl-g++

# (Optional fallback symlink in case cc-rs still probes by target triple)
RUN ln -sf "$(which musl-g++)" /usr/local/bin/x86_64-unknown-linux-musl-g++

# Copy just the workspace manifest and crate directories, so `cargo fetch` can use the cache
COPY Cargo.toml ./
COPY crates/ ./crates/

# Pre-fetch deps for the MUSL target
RUN cargo fetch --target x86_64-unknown-linux-musl

# Now copy in the rest of the source (including each crate’s src/)
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

# Copy in the two static binaries
COPY --from=builder /app/target/x86_64-unknown-linux-musl/release/arycal    /app/arycal
COPY --from=builder /app/target/x86_64-unknown-linux-musl/release/arycal-gui /app/arycal-gui

ENV PATH="/app:$PATH"
