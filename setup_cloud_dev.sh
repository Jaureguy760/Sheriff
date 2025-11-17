#!/bin/bash
# Sheriff Cloud Development Environment Setup
# One-command setup for Claude Cloud or similar environments
#
# For rust-optimization-upgrade branch

set -e

BRANCH="${1:-rust-optimization-upgrade}"

echo "================================================"
echo "Sheriff Cloud Development Setup"
echo "Branch: $BRANCH"
echo "================================================"
echo ""

# Check Python
if ! command -v python3 &> /dev/null; then
    echo "ERROR: Python 3 not found"
    exit 1
fi
PYTHON_VERSION=$(python3 --version)
echo "[OK] $PYTHON_VERSION"

# Download test data if not present
if [ ! -f "test_data/test_chr19.bam" ]; then
    echo ""
    echo "[1/5] Downloading test data (23MB)..."
    # Try branch-specific release first, then latest
    wget -q --show-progress -O test_data_chr19.tar.gz \
        "https://github.com/BradBalderson/Sheriff/releases/download/rust-opt-v1/test_data_chr19.tar.gz" \
        2>&1 | head -20 || \
    wget -q --show-progress -O test_data_chr19.tar.gz \
        "https://github.com/BradBalderson/Sheriff/releases/latest/download/test_data_chr19.tar.gz" \
        2>&1 | head -20 || {
            echo "Download failed. You may need to upload test_data_chr19.tar.gz manually."
            echo "Or create a GitHub release with the test data attached."
            echo ""
            echo "To create release:"
            echo "  1. Push branch: git push origin rust-optimization-upgrade"
            echo "  2. Go to: https://github.com/BradBalderson/Sheriff/releases/new"
            echo "  3. Tag: rust-opt-v1, Target: rust-optimization-upgrade branch"
            echo "  4. Attach: test_data_chr19.tar.gz"
        }

    if [ -f "test_data_chr19.tar.gz" ]; then
        tar xzf test_data_chr19.tar.gz
        rm test_data_chr19.tar.gz
        echo "  [OK] Test data extracted"
    fi
else
    echo "[1/5] Test data already present"
fi

# Decompress reference if needed
if [ -f "test_data/chr19.fa.gz" ] && [ ! -f "test_data/chr19.fa" ]; then
    echo ""
    echo "[2/5] Decompressing reference genome..."
    gunzip -k test_data/chr19.fa.gz
    gunzip -k test_data/chr19.gtf.gz 2>/dev/null || true
    echo "  [OK] Reference decompressed (58MB)"
else
    echo "[2/5] Reference already decompressed"
fi

# Create venv if needed
if [ ! -d "venv" ]; then
    echo ""
    echo "[3/5] Creating Python virtual environment..."
    python3 -m venv venv
    echo "  [OK] Virtual environment created"
else
    echo "[3/5] Virtual environment exists"
fi

# Activate and install
source venv/bin/activate
echo ""
echo "[4/5] Installing dependencies..."
pip install --upgrade pip -q
pip install -e . -q 2>&1 | tail -5
pip install maturin pytest -q
echo "  [OK] Sheriff and dependencies installed"

# Build Rust if cargo available
echo ""
if command -v cargo &> /dev/null; then
    echo "[5/5] Building Rust acceleration..."
    cd sheriff-rs
    maturin develop --release 2>&1 | grep -E "Built|Finished|warning" | head -10
    cd ..
    echo "  [OK] Rust acceleration built"
else
    echo "[5/5] Cargo not found - skipping Rust build"
    echo "  To install Rust: curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh"
    echo "  Rust acceleration provides 1.5-50x speedup on key operations"
fi

# Validate
echo ""
echo "================================================"
echo "Validating Installation"
echo "================================================"

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
    print('  [SKIP] sheriff_rs not available (Rust not built)')

try:
    import pysam
    print('  [OK] pysam')
except ImportError:
    print('  [FAIL] pysam not installed')
"

# Run validation tests
echo ""
echo "Running validation tests..."
if [ -f "test_data/ci_validation.py" ]; then
    python test_data/ci_validation.py 2>&1 | tail -20
fi

echo ""
echo "================================================"
echo "Setup Complete!"
echo "================================================"
echo ""
echo "To activate environment:"
echo "  source venv/bin/activate"
echo ""
echo "Quick tests:"
echo "  python test_data/ci_validation.py         # Rust function tests (1 sec)"
echo "  python test_data/validate_chr19_test.py   # Data integrity (2 sec)"
echo ""
echo "Development workflow:"
echo "  1. Edit Python: vim sheriff/count_t7.py"
echo "  2. Test: python test_data/ci_validation.py"
echo "  3. Edit Rust: vim sheriff-rs/src/edit_clustering.rs"
echo "  4. Rebuild: cd sheriff-rs && maturin develop --release && cd .."
echo "  5. Commit: git add -A && git commit -m 'changes'"
echo ""
echo "Test data: 42k reads on chr19 with reference genome"
echo "Ready for development!"
