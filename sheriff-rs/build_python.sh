#!/bin/bash
# Build script for Sheriff-rs Python bindings
#
# This script builds and installs the Python bindings for sheriff-rs
# using maturin.

set -e  # Exit on error

echo "=========================================="
echo "  Sheriff-rs Python Bindings Build"
echo "=========================================="
echo ""

# Check if we're in the right directory
if [ ! -f "Cargo.toml" ]; then
    echo "Error: Cargo.toml not found!"
    echo "Please run this script from the sheriff-rs directory."
    exit 1
fi

# Check for Rust
if ! command -v cargo &> /dev/null; then
    echo "Error: Rust is not installed!"
    echo "Install it from: https://rustup.rs/"
    echo "  curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh"
    exit 1
fi

echo "✓ Rust compiler found: $(rustc --version)"

# Check for Python
if ! command -v python3 &> /dev/null; then
    echo "Error: Python 3 is not installed!"
    exit 1
fi

echo "✓ Python found: $(python3 --version)"

# Check for maturin
if ! command -v maturin &> /dev/null; then
    echo ""
    echo "Maturin not found. Installing..."
    pip install maturin

    if [ $? -ne 0 ]; then
        echo "Error: Failed to install maturin"
        echo "Try: pip install --user maturin"
        exit 1
    fi
fi

echo "✓ Maturin found: $(maturin --version)"
echo ""

# Determine build mode
BUILD_MODE="${1:-develop}"

case "$BUILD_MODE" in
    develop|dev)
        echo "Building in DEVELOPMENT mode..."
        echo "This will install an editable version of the module."
        echo ""
        maturin develop --release --features python

        if [ $? -eq 0 ]; then
            echo ""
            echo "=========================================="
            echo "  Build successful!"
            echo "=========================================="
            echo ""
            echo "The module is now installed and ready to use."
            echo "Test it with:"
            echo "  python3 -c 'import sheriff_rs; print(sheriff_rs.__version__)'"
            echo ""
            echo "Run the demo:"
            echo "  python3 examples/python_demo.py"
        else
            echo ""
            echo "Build failed!"
            exit 1
        fi
        ;;

    wheel|release)
        echo "Building WHEEL for distribution..."
        echo ""
        maturin build --release --features python

        if [ $? -eq 0 ]; then
            echo ""
            echo "=========================================="
            echo "  Wheel built successfully!"
            echo "=========================================="
            echo ""
            echo "Wheel file created in: target/wheels/"
            ls -lh target/wheels/*.whl | tail -1
            echo ""
            echo "To install the wheel:"
            echo "  pip install target/wheels/sheriff_rs-*.whl"
        else
            echo ""
            echo "Build failed!"
            exit 1
        fi
        ;;

    *)
        echo "Usage: $0 [develop|wheel]"
        echo ""
        echo "Modes:"
        echo "  develop  - Build and install in development mode (default)"
        echo "  wheel    - Build a redistributable wheel package"
        echo ""
        echo "Examples:"
        echo "  $0              # Development build"
        echo "  $0 develop      # Development build"
        echo "  $0 wheel        # Production wheel"
        exit 1
        ;;
esac
