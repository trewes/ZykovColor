#!/usr/bin/env bash
set -euo pipefail

# Parallelism
JOBS="${JOBS:-$(nproc)}"

# Layout
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEPS="${ROOT}/deps"
mkdir -p "${DEPS}"

# Helpers
msg() { printf "\033[1;36m[build]\033[0m %s\n" "$*"; }
die() { printf "\033[1;31m[error]\033[0m %s\n" "$*" >&2; exit 1; }

need() { command -v "$1" >/dev/null 2>&1 || die "Missing tool: $1"; }
need git
need make
need cmake
need c++


# CaDiCaL
CADICAL_TAG="rel-2.1.2"
if [[ ! -d "${DEPS}/cadical" ]]; then
  msg "Cloning CaDiCaL ${CADICAL_TAG}…"
  git clone --branch "${CADICAL_TAG}" --depth=1 https://github.com/arminbiere/cadical "${DEPS}/cadical"
fi
if [[ ! -f "${DEPS}/cadical/build/libcadical.a" ]]; then
  msg "Building CaDiCaL…"
  pushd "${DEPS}/cadical" >/dev/null
    ./configure
    make -j"${JOBS}"
  popd >/dev/null
fi


# install a GMP version (needed for open-wbo) or use from GMP_ROOT
GMP_VERSION="6.3.0"
GMP_ROOT_DEFAULT="${DEPS}/gmp"
if [[ "${GMP_ROOT:-}" == "" ]]; then
  msg "GMP_ROOT not set; will build a private GMP ${GMP_VERSION} under ${GMP_ROOT_DEFAULT}"
  GMP_ROOT="${GMP_ROOT_DEFAULT}"
  if [[ ! -e "${GMP_ROOT}/include/gmpxx.h" ]]; then
    need curl
    pushd "${DEPS}" >/dev/null
      if [[ ! -d "gmp-${GMP_VERSION}" ]]; then
        msg "Downloading GMP ${GMP_VERSION}…"
        curl -fsSLo "gmp-${GMP_VERSION}.tar.xz" "https://gmplib.org/download/gmp/gmp-${GMP_VERSION}.tar.xz"
        tar xf "gmp-${GMP_VERSION}.tar.xz"
      fi
      cd "gmp-${GMP_VERSION}"
      msg "Configuring GMP (with C++)…"
      ./configure --prefix="${GMP_ROOT}" --enable-cxx --libdir="${GMP_ROOT}/lib"
      msg "Building GMP…"
      make -j"${JOBS}"
      msg "Installing GMP to ${GMP_ROOT}…"
      make install
    popd >/dev/null
  fi
else
  msg "Using existing GMP at GMP_ROOT=${GMP_ROOT}"
fi

[[ -f "${GMP_ROOT}/include/gmpxx.h" ]] || die "gmpxx.h not found in ${GMP_ROOT}/include (build/point GMP correctly)"
GMP_LIBDIR="${GMP_ROOT}/lib"


# OpenWBO
if [[ ! -d "${DEPS}/open-wbo" ]]; then
  msg "Cloning OpenWBO"
  git clone --depth=1 https://github.com/sat-group/open-wbo.git "${DEPS}/open-wbo"
  # some variants use git submodules; safe to try:
  (cd "${DEPS}/open-wbo" && git submodule update --init --recursive || true)
fi

# Build with GMP include/lib. We pass paths via CPATH/LFLAGS (the Makefile honors LFLAGS; CPATH is picked by gcc).
if [[ ! -f "${DEPS}/open-wbo/lib/lib.a" && ! -f "${DEPS}/open-wbo/lib.a" ]]; then
  OPENWBO_SOLVER="glucose4.1"
  msg "Building OpenWBO (solver: ${OPENWBO_SOLVER})…"
  pushd "${DEPS}/open-wbo" >/dev/null
    # ensure solver headers exist
    [[ -f "solvers/${OPENWBO_SOLVER}/core/Solver.h" ]] || die "Missing solver headers for ${OPENWBO_SOLVER}"
    CPATH="${GMP_ROOT}/include:${CPATH:-}" \
    LFLAGS="-L${GMP_LIBDIR} -Wl,-rpath,${GMP_LIBDIR} ${LFLAGS:-}" \
    PWD="${DEPS}/open-wbo" \
    make -j"${JOBS}" SOLVER="${OPENWBO_SOLVER}" r libr
  popd >/dev/null
fi

if [[ ! -f "${DEPS}/open-wbo/lib.a" ]]; then
  die "OpenWBO 'libr' did not produce lib.a"
fi

# ExactColors (needs GUROBI_HOME)
if [[ -z "${GUROBI_HOME:-}" ]]; then
  die "GUROBI_HOME is not set; required to build ExactColors."
fi
if [[ ! -d "${DEPS}/exactcolors" ]]; then
  msg "Cloning exactcolors"
  git clone --depth=1 https://github.com/heldstephan/exactcolors.git "${DEPS}/exactcolors"
fi
# Build the libs and the 'color' binary
if [[ ! -f "${DEPS}/exactcolors/libexactcolor.a" || ! -x "${DEPS}/exactcolors/color" ]]; then
  msg "Building exactcolors…"
  pushd "${DEPS}/exactcolors" >/dev/null
    make clean || true
    export GUROBI_HOME
    make -j"${JOBS}" libexactcolor.a
    (cd mwis_sewell && make -j"${JOBS}" libsewell.a)
    # ensure 'color' is built
    make -j"${JOBS}" color || true
  popd >/dev/null
fi

# CliSAT binary
mkdir -p "${DEPS}/clisat/bin"
if [[ ! -x "${DEPS}/clisat/bin/CliSAT_2024_07_06" ]]; then
  msg "Fetching CliSAT binary…"
  curl -fsSLo "${DEPS}/clisat/bin/CliSAT_2024_07_06" \
    "https://raw.githubusercontent.com/psanse/CliSAT/29d8e68cd29affe55cc9d81fd2ab226fe8652cfa/bin/CliSAT_2024_07_06"
  chmod +x "${DEPS}/clisat/bin/CliSAT_2024_07_06"
fi

# Configure & build ZykovColor (using its cmake)
BUILD_DIR="${ROOT}/cmake-build"
mkdir -p "${BUILD_DIR}"
msg "Configuring ZykovColor…"

cmake -S "${ROOT}" -B "${BUILD_DIR}" \
  -DCMAKE_BUILD_TYPE=Release \
  -DCADICAL_ROOT_DIR="${DEPS}/cadical" \
  -DOPENWBO_ROOT_DIR="${DEPS}/open-wbo" \
  -DEXACTCOLORS_ROOT_DIR="${DEPS}/exactcolors" \
  -DCLISAT_BINARY_PATH="${DEPS}/clisat/bin/CliSAT_2024_07_06" \
  -DGUROBI_HOME="${GUROBI_HOME}" \
  -DBOOST_ROOT="${BOOST_ROOT}" \
  -DCMAKE_BUILD_TYPE=Release

msg "Building ZykovColor…"
cmake --build "${BUILD_DIR}" -j"${JOBS}"

msg "Done. Artifacts:"
echo "  CaDiCaL:        ${DEPS}/cadical (lib/include)"
echo "  OpenWBO:        ${DEPS}/openwbo (lib + include -> source tree)"
echo "  exactcolors:    ${DEPS}/exactcolors/root (binary 'color' inside)"
echo "  CliSAT:         ${DEPS}/clisat/bin/CliSAT_2024_07_06"
echo "  ZykovColor bin: $(find "${BUILD_DIR}" -maxdepth 2 -type f -name ZykovColor -print 2>/dev/null || true)"
