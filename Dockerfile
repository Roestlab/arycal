# Stage 1: Build both binaries with the MUSL toolchain
FROM clux/muslrust:1.85.1-stable AS builder

WORKDIR /app

# 1. Only install pkg-config & OpenSSL headers (musl-g++ is already present)
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
      pkg-config \
      libssl-dev && \
    rm -rf /var/lib/apt/lists/*

# 2. Copy workspace manifests and the crates/ folder for caching
COPY Cargo.toml ./
COPY crates/ ./crates/

# 3. Pre-fetch everything for the MUSL target
RUN cargo fetch --target x86_64-unknown-linux-musl

# 4. Copy the rest of your source (including each crate’s src/)
COPY . .

# 5. Build optimized & stripped binaries
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
