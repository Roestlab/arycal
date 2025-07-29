# Stage 1: Build binaries
FROM rust:1.85-slim AS builder

WORKDIR /app

# Install musl toolchain + C++ support
RUN apt-get update && \
    apt-get install -y musl-tools g++-musl pkg-config libssl-dev && \
    rustup target add x86_64-unknown-linux-musl

# Copy source
COPY . .

# Build optimized + stripped
ENV RUSTFLAGS="-C strip=symbols -C lto=yes"
RUN cargo build --release --target x86_64-unknown-linux-musl --bin arycal
RUN cargo build --release --target x86_64-unknown-linux-musl --bin arycal-gui

# Stage 2: Runtime image
FROM debian:bullseye-slim

RUN apt-get update && \
    apt-get install -y --no-install-recommends procps ca-certificates && \
    update-ca-certificates && \
    rm -rf /var/lib/apt /var/lib/dpkg /var/lib/cache /var/lib/log

WORKDIR /app
COPY --from=builder /app/target/x86_64-unknown-linux-musl/release/arycal /app/arycal
COPY --from=builder /app/target/x86_64-unknown-linux-musl/release/arycal-gui /app/arycal-gui

ENV PATH="/app:$PATH"
