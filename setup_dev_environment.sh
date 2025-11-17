#!/bin/bash
# Sheriff Development Environment Setup
# Sets up complete dev environment including Rust acceleration

set -e

echo "Sheriff Development Environment Setup"
echo "======================================"
echo ""

# Check Python version
PYTHON_VERSION=$(python3 --version 2>&1 | awk '{print $2}')
PYTHON_MAJOR=$(echo $PYTHON_VERSION | cut -d. -f1)
PYTHON_MINOR=$(echo $PYTHON_VERSION | cut -d. -f2)

echo "[1/6] Checking Python version..."
if [ "$PYTHON_MAJOR" -lt 3 ] || ([ "$PYTHON_MAJOR" -eq 3 ] && [ "$PYTHON_MINOR" -lt 10 ]); then
    echo "ERROR: Python 3.10+ required, found $PYTHON_VERSION"
    exit 1
fi
echo "  [OK] Python $PYTHON_VERSION"

# Check for required system tools
echo ""
echo "[2/6] Checking system dependencies..."
MISSING_TOOLS=""

if ! command -v samtools >/dev/null 2>&1; then
    MISSING_TOOLS="$MISSING_TOOLS samtools"
fi

if ! command -v cargo >/dev/null 2>&1; then
    MISSING_TOOLS="$MISSING_TOOLS cargo/rustc"
fi

if [ -n "$MISSING_TOOLS" ]; then
    echo "  WARNING: Missing optional tools:$MISSING_TOOLS"
    echo "  These are needed for full functionality but not required for basic setup."
else
    echo "  [OK] All system tools found"
fi

# Create virtual environment
echo ""
echo "[3/6] Setting up Python virtual environment..."
if [ ! -d "venv" ]; then
    python3 -m venv venv
    echo "  [OK] Created virtual environment"
else
    echo "  [OK] Virtual environment already exists"
fi

# Activate venv
source venv/bin/activate
echo "  [OK] Activated virtual environment"

# Install Python dependencies
echo ""
echo "[4/6] Installing Python dependencies..."
pip install --upgrade pip setuptools wheel >/dev/null

# Core dependencies
pip install pysam numpy pandas anndata scanpy >/dev/null 2>&1 || {
    echo "  Installing core packages one by one..."
    pip install pysam >/dev/null
    pip install numpy >/dev/null
    pip install pandas >/dev/null
    pip install anndata >/dev/null
    pip install scanpy >/dev/null
}
echo "  [OK] Core dependencies installed"

# Install maturin for Rust builds
pip install maturin >/dev/null
echo "  [OK] Maturin installed (for Rust builds)"

# Install Sheriff in development mode
echo ""
echo "[5/6] Installing Sheriff in development mode..."
pip install -e . >/dev/null
echo "  [OK] Sheriff installed"

# Build Rust acceleration
echo ""
echo "[6/6] Building Rust acceleration..."
if command -v cargo >/dev/null 2>&1; then
    cd sheriff-rs
    maturin develop --release 2>&1 | head -20
    cd ..
    echo "  [OK] Rust acceleration built"
else
    echo "  SKIPPED: cargo not found (install Rust: https://rustup.rs/)"
fi

# Verify installation
echo ""
echo "======================================"
echo "Verifying Installation"
echo "======================================"

python3 -c "
import sys
sys.path.insert(0, '.')

print('Testing imports...')
try:
    from sheriff import count_t7
    print('  [OK] sheriff.count_t7')
except Exception as e:
    print(f'  [FAIL] sheriff.count_t7: {e}')

try:
    from sheriff import helpers
    print('  [OK] sheriff.helpers')
except Exception as e:
    print(f'  [FAIL] sheriff.helpers: {e}')

try:
    import sheriff_rs
    print('  [OK] sheriff_rs (Rust acceleration)')
except ImportError:
    print('  [SKIP] sheriff_rs (Rust not built)')
"

echo ""
echo "======================================"
echo "Setup Complete!"
echo "======================================"
echo ""
echo "Next steps:"
echo ""
echo "1. Activate the environment:"
echo "   source venv/bin/activate"
echo ""
echo "2. Download test data (if not already):"
echo "   cd test_data && bash download_reference.sh"
echo ""
echo "3. Run validation:"
echo "   python test_data/run_validation.py"
echo ""
echo "4. Run a quick test:"
echo "   python -m sheriff.count_t7 --help"
echo ""
echo "For Claude Cloud deployment:"
echo "  - Package test_data directory (without reference.fa/gtf)"
echo "  - Use GitHub releases for distribution"
echo "  - See test_data/README.md for details"
