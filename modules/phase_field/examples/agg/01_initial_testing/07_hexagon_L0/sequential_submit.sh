#!/usr/bin/env bash

# Iso

python ~/scripts/phase_field_run.py -i 07_iso_hexagon_L0.i -n 6 -c='L0val=1'

python ~/scripts/phase_field_run.py -i 07_iso_hexagon_L0.i -n 6 -c='L0val=2'

python ~/scripts/phase_field_run.py -i 07_iso_hexagon_L0.i -n 6 -c='L0val=3'

python ~/scripts/phase_field_run.py -i 07_iso_hexagon_L0.i -n 6 -c='L0val=4'

# Aniso

python ~/scripts/phase_field_run.py -i 07_hexagon_L0.i -n 6 -c='L0val=1'

python ~/scripts/phase_field_run.py -i 07_hexagon_L0.i -n 6 -c='L0val=2'

python ~/scripts/phase_field_run.py -i 07_hexagon_L0.i -n 6 -c='L0val=3'

python ~/scripts/phase_field_run.py -i 07_hexagon_L0.i -n 6 -c='L0val=4'
