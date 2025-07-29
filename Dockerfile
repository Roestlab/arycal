# Stage 1: Build binaries with musl toolchain
FROM clux/muslrust:1.85.0 AS builder

WORKDIR /app

# Install extra build deps (if needed for crates like openssl)
RUN apt-get update && \
    apt-get install -y pkg-config libssl-dev && \
    rustup target add x86_64-unknown-linux-musl

# Copy source code
COPY . .

# Build optimized & stripped binaries
ENV RUSTFLAGS="-C strip=symbols -C lto=yes"
RUN cargo build --release --target x86_64-unknown-linux-musl --bin arycal
RUN cargo build --release --target x86_64-unknown-linux-musl --bin arycal-gui

# Stage 2: Minimal runtime image
FROM debian:bullseye-slim

RUN apt-get update && \
    apt-get install -y --no-install-recommends procps ca-certificates && \
    update-ca-certificates && \
    rm -rf /var/lib/apt /var/lib/dpkg /var/lib/cache /var/lib/log

WORKDIR /app
COPY --from=builder /app/target/x86_64-unknown-linux-musl/release/arycal /app/arycal
COPY --from=builder /app/target/x86_64-unknown-linux-musl/release/arycal-gui /app/arycal-gui

ENV PATH="/app:$PATH"
