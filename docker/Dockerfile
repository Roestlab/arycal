FROM rust:latest AS builder

# Set working directory
WORKDIR /app

# Clone the repository
RUN git clone https://github.com/singjc/arycal.git .

# Build the binaries
RUN cargo build --release --bin arycal
RUN cargo build --release --bin arycal-gui

# Create a minimal runtime image
FROM debian:bookworm-slim

# Install necessary dependencies for GUI (if needed)
RUN apt-get update && apt-get install -y \
    libx11-dev libgl1-mesa-glx libgtk-3-dev \
    && rm -rf /var/lib/apt/lists/*

# Copy binaries from the builder stage
COPY --from=builder /app/target/release/arycal /usr/local/bin/arycal
COPY --from=builder /app/target/release/arycal-gui /usr/local/bin/arycal-gui

# Allow the user to specify which binary to run
ENTRYPOINT []
CMD ["/usr/local/bin/arycal"]

