#!/bin/bash

rank="${OMPI_COMM_WORLD_RANK:-${PMI_RANK:-${PMIX_RANK:-unknown}}}"
label="${MEM_LABEL:-run}"

exec /usr/bin/time -l "$@" 2>"${label}_rank_${rank}.log"
