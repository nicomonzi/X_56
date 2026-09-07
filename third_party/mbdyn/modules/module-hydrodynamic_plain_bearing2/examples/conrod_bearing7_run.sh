#!/bin/sh -f

## FIXME: conrod_bearing7.mbdyn disabled due to timeout of 30s

MBDYN_EXEC="${MBDYN_EXEC:-mbdyn}"

exec ${MBDYN_EXEC} $*
