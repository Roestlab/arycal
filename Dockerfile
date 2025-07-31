# 1) Builder stage: plain Ubuntu + Rust via rustup
FROM ubuntu:22.04 AS builder
WORKDIR /app

# Install system deps: curl for rustup, pkg-config, OpenSSL dev, C/C++ toolchains,
# clang/clang++ (needed by DuckDB), and musl-gcc/g++
RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
      curl \
      build-essential \
      pkg-config \
      libssl-dev \
      clang \
      clang++ \
      musl-tools && \
    rm -rf /var/lib/apt/lists/*

# Install Rust toolchain
RUN curl https://sh.rustup.rs -sSf | sh -s -- -y --default-toolchain 1.85.0
ENV PATH="/root/.cargo/bin:${PATH}"

# Add the MUSL target
RUN rustup target add x86_64-unknown-linux-musl

# Create compiler symlinks so that cc-rs probes (e.g. libduckdb-sys) will find the
# right compiler for the MUSL target triple.
RUN ln -sf "$(which musl-gcc)" /usr/local/bin/x86_64-unknown-linux-musl-gcc && \
    ln -sf "$(which musl-g++)"  /usr/local/bin/x86_64-unknown-linux-musl-g++

# Copy just the manifests and your crates folder so that Docker cache
# can re-use `cargo fetch` when the code hasn’t changed.
COPY Cargo.toml ./
COPY crates/ ./crates/

# Pre-fetch all crates for the MUSL target
RUN cargo fetch --target x86_64-unknown-linux-musl

# Now copy the rest of the code and build your binaries
COPY . .

# Enable LTO and strip symbols for small static binaries
ENV RUSTFLAGS="-C lto=thin -C strip=symbols"

RUN cargo build --release --target x86_64-unknown-linux-musl --bin arycal
RUN cargo build --release --target x86_64-unknown-linux-musl --bin arycal-gui

# 2) Final “runtime” stage: slim Debian image with only the two statically-linked binaries
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
