#!/usr/bin/env bash
# Exit immediately if a command exits with a non-zero status
set -e 

# 1. Define Paths (Configurable via Environment Variables)
# Users can export these before running the script, otherwise defaults are used.
GAMSDIST="${GAMSDIST:-$HOME/gamsdist}"
WORKSPACE="${WORKSPACE:-$(pwd)}"
CUOPT_VERSION="${CUOPT_VERSION:-26.06}"
CUOPT_HASH="${CUOPT_HASH:-12345}"

# 2. Dynamically find the site-packages directory to avoid hardcoding the Python version
SITE_PACKAGES=$(find "$WORKSPACE/venv/lib" -maxdepth 2 -type d -name "site-packages" | head -n 1)

if [ -z "$SITE_PACKAGES" ]; then
    echo "Error: Could not find site-packages in $WORKSPACE/venv/lib."
    echo "Please ensure you have created the virtual environment and installed dependencies."
    exit 1
fi

# 3. Export specific build paths
export GAMSCAPI="$GAMSDIST/apifiles/C/api"
export CUOPT="$SITE_PACKAGES/libcuopt"
export JITLINK="$SITE_PACKAGES/nvidia/cu13/lib"

echo "Building GAMS cuOpt solver link..."
echo "Using GAMSDIST: $GAMSDIST"
echo "Using SITE_PACKAGES: $SITE_PACKAGES"

# 4. Compile the solver link
gcc -g3 -O0 -Wall gmscuopt.c -o gmscuopt-cu13.out \
    -DCUOPT_VERSION=\"$CUOPT_VERSION\" \
    -DCUOPT_HASH=\"$CUOPT_HASH\" \
    -I "$GAMSCAPI" "$GAMSCAPI/gmomcc.c" "$GAMSCAPI/optcc.c" "$GAMSCAPI/gevmcc.c" \
    -I "$CUOPT/include" "$JITLINK/libnvJitLink.so.13" \
    -L "$CUOPT/lib64" -lcuopt

# 5. Patch RPATH
patchelf --set-rpath \$ORIGIN gmscuopt-cu13.out

# 6. Move output and copy dependencies
echo "Copying binaries and dependencies to $GAMSDIST..."

# Create target directory if it doesn't exist
mkdir -p "$GAMSDIST"

mv gmscuopt-cu13.out "$GAMSDIST/gmscuopt.out"

cp "$CUOPT/lib64/libcuopt.so" "$GAMSDIST/"
cp "$SITE_PACKAGES"/libcuopt_cu13.libs/libgomp-*.so* "$GAMSDIST/"
cp "$SITE_PACKAGES"/libcuopt_cu13.libs/libtbb-*.so* "$GAMSDIST/"
cp "$SITE_PACKAGES"/libcuopt_cu13.libs/libtbbmalloc-*.so* "$GAMSDIST/"
cp "$SITE_PACKAGES"/rapids_logger/lib64/librapids_logger.so "$GAMSDIST/"
cp "$SITE_PACKAGES"/librmm/lib64/librmm.so "$GAMSDIST/"
cp "$SITE_PACKAGES"/nvidia/cu13/lib/libcudss.so* "$GAMSDIST/"
cp "$SITE_PACKAGES"/nvidia/cu13/lib/libcudss_mtlayer_gomp.so* "$GAMSDIST/"
cp "$SITE_PACKAGES"/nvidia/cu13/lib/libnvJitLink.so* "$GAMSDIST/"
cp "$SITE_PACKAGES"/nvidia/cu13/lib/libcublas.so* "$GAMSDIST/"
cp "$SITE_PACKAGES"/nvidia/cu13/lib/libcublasLt.so* "$GAMSDIST/"
cp "$SITE_PACKAGES"/nvidia/cu13/lib/libcurand.so* "$GAMSDIST/"
cp "$SITE_PACKAGES"/nvidia/cu13/lib/libcusolver.so* "$GAMSDIST/"
cp "$SITE_PACKAGES"/nvidia/cu13/lib/libcusparse.so* "$GAMSDIST/"
cp "$SITE_PACKAGES"/libcuopt_cu13.libs/libcares-*.so* "$GAMSDIST/"

# Copy assets (suppress errors if directory is empty or missing)
cp -r assets/* "$GAMSDIST/" 2>/dev/null || true

echo "Build and copy completed successfully."