#!/usr/bin/env bash
set -euo pipefail

cd ~/Aetherion-Astrophysics-Suite-git
git pull
flatpak-builder --user --install --force-clean build-flatpak \
    flatpak/io.github.ArandomHitman.AetherionSuite.json
flatpak run io.github.ArandomHitman.AetherionSuite


# For lazy people who don't wanna type everything line by line (me fr)

# also fun fact the "force-clean" option is actually required to get the latest changes to be reflected in the build, otherwise it just reuses the old build artifacts 
# and doesn't actually update anything. so it's not just for cleanliness, it's actually necessary for development. --- IGNORE ---