#!/bin/sh -f

## FIXME: thermal_problem_17_1.mbdyn disabled due to timeout of 30s

MBDYN_EXEC="${MBDYN_EXEC:-mbdyn}"

exec ${MBDYN_EXEC} $*
