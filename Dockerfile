# Stage 1: Build both binaries with MUSL + Debian/Ubuntu toolchain
FROM ubuntu:22.04 AS builder
WORKDIR /app

# Install system deps: curl for rustup, pkg-config, OpenSSL dev,
# build-essential (gcc/g++), and musl-tools for static linking
RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
      curl \
      build-essential \
      pkg-config \
      libssl-dev \
      musl-tools && \
    rm -rf /var/lib/apt/lists/*

# Install Rust 1.85.0 via rustup
RUN curl https://sh.rustup.rs -sSf | sh -s -- -y --default-toolchain 1.85.0
ENV PATH="/root/.cargo/bin:${PATH}"

# Add the MUSL target and symlink compilers so cc-rs will find them
RUN rustup target add x86_64-unknown-linux-musl && \
    ln -sf "$(which musl-gcc)" /usr/local/bin/x86_64-unknown-linux-musl-gcc && \
    ln -sf "$(which musl-g++)" /usr/local/bin/x86_64-unknown-linux-musl-g++

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
