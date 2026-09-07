#!/bin/sh -f

## FIXME: thermal_problem8_2.mbdyn disabled due to timeout of 30s

MBDYN_EXEC="${MBDYN_EXEC:-mbdyn}"

exec ${MBDYN_EXEC} $*
