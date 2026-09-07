#!/bin/sh -f

## FIXME: flow_problem12_2.mbdyn disabled due to timeout of 10s

MBDYN_EXEC="${MBDYN_EXEC:-mbdyn}"

exec ${MBDYN_EXEC} $*
