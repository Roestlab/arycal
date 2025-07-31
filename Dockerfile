# Use Cross's pre-built MUSL image that has both C and C++ toolchains
FROM rustembedded/cross:x86_64-unknown-linux-musl AS builder
WORKDIR /app

# Copy your workspace, fetch deps then build
COPY Cargo.toml ./
COPY crates/ ./crates/
RUN cargo fetch --target x86_64-unknown-linux-musl

COPY . .
ENV RUSTFLAGS="-C strip=symbols -C lto=yes"
RUN cargo build --release --target x86_64-unknown-linux-musl --bin arycal
RUN cargo build --release --target x86_64-unknown-linux-musl --bin arycal-gui

# Now pack into a tiny Debian runtime
FROM debian:bullseye-slim
RUN apt-get update && apt-get install -y procps ca-certificates \
    && update-ca-certificates && rm -rf /var/lib/apt/lists/*

WORKDIR /app
COPY --from=builder /app/target/x86_64-unknown-linux-musl/release/arycal    /app/arycal
COPY --from=builder /app/target/x86_64-unknown-linux-musl/release/arycal-gui /app/arycal-gui

ENV PATH="/app:$PATH"
ENTRYPOINT ["/app/arycal"]
