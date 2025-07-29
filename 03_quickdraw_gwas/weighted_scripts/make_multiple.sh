#!/usr/bin/env bash
# ---------------------------------------------------------------------------
#  make_multiple.sh  –  create chr{1..22}_weighted.sh from chr01_weighted.sh
#  Works on both GNU sed (Linux) and BSD sed (macOS)
# ---------------------------------------------------------------------------
set -euo pipefail

TEMPLATE="chr01_weighted.sh"
[[ -s $TEMPLATE ]] || { echo "❌  $TEMPLATE is missing or empty"; exit 1; }

# Helper: portable in‑place sed
SED_INPLACE () {
  local expr=$1 file=$2
  if sed --version >/dev/null 2>&1; then        # GNU sed
    sed -i -e "$expr" "$file"
  else                                          # BSD sed (macOS)
    sed -i '' -e "$expr" "$file"
  fi
}

for CHR in {1..22}; do
  NEW="chr${CHR}_weighted.sh"
  cp "$TEMPLATE" "$NEW"

  SED_INPLACE "s/chr01/chr${CHR}/g"           "$NEW"
  SED_INPLACE "s/step2_chr01_/step2_chr${CHR}_/g" "$NEW"

  chmod +x "$NEW"
  echo "✅  created $NEW"
done
