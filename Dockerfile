# Stage 1: Build binaries with cross
FROM rust:1.85-slim AS builder

WORKDIR /app

# Install cross (prebuilt Docker images handle musl + C++ toolchains)
RUN cargo install cross --git https://github.com/cross-rs/cross

# Copy source code
COPY . .

# Build optimized & stripped binaries using cross
ENV RUSTFLAGS="-C strip=symbols -C lto=yes"
RUN cross build --release --target x86_64-unknown-linux-musl --bin arycal
RUN cross build --release --target x86_64-unknown-linux-musl --bin arycal-gui

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
