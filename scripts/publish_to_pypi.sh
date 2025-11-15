#!/bin/bash
# Sheriff PyPI Publication Script
# Run this script to publish Sheriff to PyPI

set -e  # Exit on error

echo "====================================="
echo "Sheriff v1.2.0 PyPI Publication"
echo "====================================="
echo

# Color codes
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Check if we're in the right directory
if [ ! -f "pyproject.toml" ]; then
    echo -e "${RED}Error: Must run from Sheriff root directory${NC}"
    exit 1
fi

# Check if build tools are installed
echo "Checking dependencies..."
python -m pip install --upgrade pip build twine -q

echo -e "${GREEN}✓${NC} Build tools installed"
echo

# Clean previous builds
echo "Cleaning previous builds..."
rm -rf dist/ build/ *.egg-info
echo -e "${GREEN}✓${NC} Clean complete"
echo

# Build package
echo "Building package..."
python -m build
echo -e "${GREEN}✓${NC} Build complete"
echo

# Show built files
echo "Built packages:"
ls -lh dist/
echo

# Check package
echo "Validating package..."
twine check dist/* || echo -e "${YELLOW}Warning: Metadata validation showed warnings (safe to ignore)${NC}"
echo

# Ask user which repository to upload to
echo "Select upload target:"
echo "  1) TestPyPI (https://test.pypi.org) - For testing"
echo "  2) Production PyPI (https://pypi.org) - PERMANENT"
echo "  3) Exit without uploading"
echo

read -p "Enter choice [1-3]: " choice

case $choice in
    1)
        echo
        echo -e "${YELLOW}Uploading to TestPyPI...${NC}"
        echo "Make sure you have TEST_PYPI_API_TOKEN set or ~/.pypirc configured"
        echo
        read -p "Continue? [y/N]: " confirm
        if [[ $confirm == [yY] ]]; then
            twine upload --repository testpypi dist/*
            echo
            echo -e "${GREEN}✓${NC} Upload complete!"
            echo
            echo "View at: https://test.pypi.org/project/sheriff/1.2.0/"
            echo
            echo "To test install:"
            echo "  pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/ sheriff"
        fi
        ;;
    2)
        echo
        echo -e "${RED}WARNING: This will publish to PRODUCTION PyPI!${NC}"
        echo -e "${RED}This action is PERMANENT and cannot be undone!${NC}"
        echo
        echo "Make sure:"
        echo "  - Version 1.2.0 has not been published before"
        echo "  - All tests pass"
        echo "  - CHANGELOG.md is updated"
        echo "  - You have PYPI_API_TOKEN set or ~/.pypirc configured"
        echo
        read -p "Type 'publish' to confirm: " confirm
        if [[ $confirm == "publish" ]]; then
            twine upload dist/*
            echo
            echo -e "${GREEN}✓${NC} Upload complete!"
            echo
            echo "View at: https://pypi.org/project/sheriff/1.2.0/"
            echo
            echo "Users can now install with:"
            echo "  pip install sheriff"
            echo
            echo -e "${GREEN}Next steps:${NC}"
            echo "  1. Update README with PyPI badge"
            echo "  2. Create GitHub release: git tag v1.2.0 && git push --tags"
            echo "  3. Announce on social media"
            echo "  4. Update documentation site"
        else
            echo -e "${YELLOW}Upload cancelled${NC}"
        fi
        ;;
    3)
        echo -e "${YELLOW}Upload cancelled${NC}"
        ;;
    *)
        echo -e "${RED}Invalid choice${NC}"
        exit 1
        ;;
esac

echo
echo "Done!"
