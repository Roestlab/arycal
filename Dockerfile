# Stage 1: compile inside a Rust container (uses glibc)
FROM rust:1.85-slim AS builder
WORKDIR /app

# Install any system deps your crate needs (e.g. for MPI, GTK, etc.)
RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
      build-essential \
      pkg-config \
      libssl-dev \
      libopenmpi-dev openmpi-bin \
      libx11-dev libxrandr-dev libxi-dev \
      libgl1-mesa-dev libegl1-mesa libgtk-3-dev && \
    rm -rf /var/lib/apt/lists/*

# Bring in your Cargo.toml / lock first for caching
COPY Cargo.toml ./
COPY crates/ ./crates/

# Pre-fetch deps for glibc target
RUN cargo fetch

# Copy the rest of your code
COPY . .

# Build optimized release binaries
RUN cargo build --release --bin arycal
RUN cargo build --release --bin arycal-gui

# Stage 2: slim runtime
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
