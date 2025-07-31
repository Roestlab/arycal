# Stage 1: Build both binaries with MUSL + Debian/Ubuntu toolchain
FROM ubuntu:22.04 AS builder
WORKDIR /app

# 1) Install system deps *including* ca-certificates
RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
      curl \
      ca-certificates \
      build-essential \
      pkg-config \
      libssl-dev \
      musl-tools && \
    update-ca-certificates && \
    rm -rf /var/lib/apt/lists/*

# 2) Install rustup (now that SSL works)
RUN curl https://sh.rustup.rs -sSf | sh -s -- -y --default-toolchain 1.85.0

# 3) Make sure `~/.cargo/bin` is on PATH for subsequent steps
ENV PATH="/root/.cargo/bin:${PATH}"

# 4) Add the MUSL target
RUN rustup target add x86_64-unknown-linux-musl

# Copy just Cargo manifests & crates dir to leverage Docker cache
COPY Cargo.toml ./
COPY crates/ ./crates/

# Pre-fetch dependencies for the MUSL target
RUN cargo fetch --target x86_64-unknown-linux-musl

# Copy the rest of your source code
COPY . .

# Build optimized, stripped static binaries
ENV RUSTFLAGS="-C lto=thin -C strip=symbols"
RUN cargo build --release --target x86_64-unknown-linux-musl --bin arycal
RUN cargo build --release --target x86_64-unknown-linux-musl --bin arycal-gui

# Stage 2: Minimal runtime image
FROM debian:bullseye-slim
RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
      ca-certificates && \
    update-ca-certificates && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /app
COPY --from=builder /app/target/x86_64-unknown-linux-musl/release/arycal    /usr/local/bin/arycal
COPY --from=builder /app/target/x86_64-unknown-linux-musl/release/arycal-gui /usr/local/bin/arycal-gui

ENTRYPOINT ["arycal"]
CMD ["--help"]
