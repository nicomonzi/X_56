#!/bin/sh -f

## FIXME: conrod_bearing6.mbdyn disabled due to timeout of 10s

MBDYN_EXEC="${MBDYN_EXEC:-mbdyn}"

exec ${MBDYN_EXEC} $*
