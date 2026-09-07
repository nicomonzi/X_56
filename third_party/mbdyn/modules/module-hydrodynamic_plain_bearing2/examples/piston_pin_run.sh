#!/bin/sh -f

## FIXME: piston_pin.mbdyn disabled due to timeout of 10s

MBDYN_EXEC="${MBDYN_EXEC:-mbdyn}"

exec ${MBDYN_EXEC} $*
