#!/usr/bin/env bash
# Build the metabinR shaded JAR and place it under inst/java/.
#
# Requires: Java 17+ and Maven on PATH.
# Usage:  ./tools/build-jar.sh
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
POM="${REPO_ROOT}/java/metabinR/pom.xml"
VERSION="$(grep -m1 -oE '<version>[^<]+</version>' "${POM}" | head -n1 | sed -E 's|</?version>||g')"

mvn -f "${POM}" -q clean package

SHADED="${REPO_ROOT}/java/metabinR/target/metabinR-${VERSION}-jar-with-dependencies.jar"
DEST_DIR="${REPO_ROOT}/inst/java"
DEST="${DEST_DIR}/metabinR-${VERSION}-jar-with-dependencies.jar"

mkdir -p "${DEST_DIR}"
# Drop any previous metabinR or MetaTarget jars so only one artefact ships.
find "${DEST_DIR}" -maxdepth 1 -type f \( -name 'metabinR-*-jar-with-dependencies.jar' -o -name 'MetaTarget-*-jar-with-dependencies.jar' \) -delete

cp "${SHADED}" "${DEST}"
echo "Built: ${DEST}"
sha256sum "${DEST}"
