#!/bin/sh
# Installs the msms and discoal command-line binaries into the active
# conda/mamba environment prefix (CONDA_PREFIX). Never uses sudo and never
# writes outside the active environment prefix.
set -eu

MSMS_URL="${MSMS_URL:-https://www.mabs.at/fileadmin/user_upload/p_mabs/msms3.2rc-b163.jar}"
MSMS_SHA256="${MSMS_SHA256-}"
DISCOAL_REF="${DISCOAL_REF-}"

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
