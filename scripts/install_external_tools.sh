#!/bin/sh
# Installs the msms and discoal command-line binaries into the active
# conda/mamba environment prefix (CONDA_PREFIX). Never uses sudo and never
# writes outside the active environment prefix.
set -eu

MSMS_URL="${MSMS_URL:-https://www.mabs.at/fileadmin/user_upload/p_mabs/msms3.2rc-b163.jar}"
MSMS_SHA256="${MSMS_SHA256-}"
DISCOAL_REF="${DISCOAL_REF-}"
RAISD_AI_REF="${RAISD_AI_REF-}"

if [ -z "${CONDA_PREFIX:-}" ]; then
  echo "CONDA_PREFIX is not set. Activate the popsimwraptor conda environment first:" >&2
  echo "  conda activate popsimwraptor" >&2
  exit 1
fi

BIN_DIR="$CONDA_PREFIX/bin"
SHARE_DIR="$CONDA_PREFIX/share/popsimwraptor"
mkdir -p "$BIN_DIR" "$SHARE_DIR/msms"

command -v java >/dev/null 2>&1 || { echo "java not found in environment; ensure openjdk is installed." >&2; exit 1; }

# ---- msms ------------------------------------------------------------------
jar="$SHARE_DIR/msms/msms.jar"
if [ ! -f "$jar" ]; then
  echo "Downloading msms..."
  tmp="$(mktemp)"
  if command -v curl >/dev/null 2>&1; then
    curl -fsSL "$MSMS_URL" -o "$tmp"
  else
    wget -q -O "$tmp" "$MSMS_URL"
  fi
  head -c4 "$tmp" | od -An -t x1 | tr -d " \n" | grep -qi 504b0304 || { echo "Downloaded msms.jar failed signature check" >&2; rm -f "$tmp"; exit 1; }
  if [ -n "$MSMS_SHA256" ]; then
    if command -v sha256sum >/dev/null 2>&1; then got="$(sha256sum "$tmp" | awk '{print $1}')"; else got="$(shasum -a 256 "$tmp" | awk '{print $1}')"; fi
    [ "$got" = "$MSMS_SHA256" ] || { echo "msms.jar checksum mismatch" >&2; rm -f "$tmp"; exit 1; }
  fi
  mv "$tmp" "$jar"
  chmod 0644 "$jar"
else
  echo "msms.jar already present."
fi

cat > "$BIN_DIR/msms" <<EOF
#!/bin/sh
set -eu
exec java -jar "$jar" "\$@"
EOF
chmod 0755 "$BIN_DIR/msms"
echo "msms ready."

# ---- discoal -----------------------------------------------------------------
if command -v git >/dev/null 2>&1 && command -v make >/dev/null 2>&1 && (command -v gcc >/dev/null 2>&1 || command -v clang >/dev/null 2>&1); then
  repo="$SHARE_DIR/discoal"
  if [ ! -d "$repo/.git" ]; then
    rm -rf "$repo"
    git clone --depth 1 https://github.com/kern-lab/discoal.git "$repo"
  else
    git -C "$repo" fetch --all --tags
    git -C "$repo" pull --ff-only
  fi
  if [ -n "$DISCOAL_REF" ]; then
    git -C "$repo" checkout --quiet "$DISCOAL_REF"
  fi
  ( cd "$repo" && make -s discoal )
  cp "$repo/discoal" "$BIN_DIR/discoal"
  chmod 0755 "$BIN_DIR/discoal"
  echo "discoal ready."
else
  echo "Compiler toolchain (git/make/gcc) not available; skipping discoal build." >&2
fi

# ---- RAiSD-AI ------------------------------------------------------------
# Non-blocking: never fails the overall install, just skips/warns if the
# toolchain or build scripts are unavailable.
if command -v git >/dev/null 2>&1 && command -v make >/dev/null 2>&1 && command -v gcc >/dev/null 2>&1; then
  repo="$SHARE_DIR/RAiSD-AI"
  if [ ! -d "$repo/.git" ]; then
    rm -rf "$repo"
    git clone --depth 1 https://github.com/alachins/RAiSD-AI.git "$repo" || true
    if [ -n "$RAISD_AI_REF" ] && [ -d "$repo/.git" ]; then
      git -C "$repo" fetch --tags || true
      git -C "$repo" checkout --quiet "$RAISD_AI_REF" || true
    fi
  else
    git -C "$repo" fetch --all --tags || true
    git -C "$repo" pull --ff-only || true
  fi

  if [ -d "$repo/.git" ]; then
    (
      set +e
      cd "$repo"
      [ -x ./compile-RAiSD-AI.sh ] && ./compile-RAiSD-AI.sh
      [ -x ./compile-RAiSD-AI-ZLIB.sh ] && ./compile-RAiSD-AI-ZLIB.sh
      exit 0
    )

    ai="$BIN_DIR/RAiSD-AI"
    zlib="$BIN_DIR/RAiSD-AI-ZLIB"
    [ -e "$repo/bin/release/RAiSD-AI" ] && cp -L "$repo/bin/release/RAiSD-AI" "$ai" && chmod 0755 "$ai"
    [ -e "$repo/bin/release/RAiSD-AI-ZLIB" ] && cp -L "$repo/bin/release/RAiSD-AI-ZLIB" "$zlib" && chmod 0755 "$zlib"
    [ ! -x "$ai" ] && [ -x "$zlib" ] && ln -sf RAiSD-AI-ZLIB "$ai"
    [ ! -x "$zlib" ] && [ -x "$ai" ] && ln -sf RAiSD-AI "$zlib"

    if [ -x "$ai" ] || [ -x "$zlib" ]; then
      echo "RAiSD-AI ready."
    else
      echo "RAiSD-AI clone succeeded but no binary was built; check compile-RAiSD-AI.sh output." >&2
    fi
  else
    echo "Could not clone RAiSD-AI; skipping." >&2
  fi
else
  echo "Compiler toolchain (git/make/gcc) not available; skipping RAiSD-AI build." >&2
fi
