#!/bin/sh -f

## FIXME: piston_blow_by.mbdyn disabled due to timeout of 30s

MBDYN_EXEC="${MBDYN_EXEC:-mbdyn}"

exec ${MBDYN_EXEC} $*
