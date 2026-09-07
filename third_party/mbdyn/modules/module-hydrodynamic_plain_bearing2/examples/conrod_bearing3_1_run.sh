#!/bin/sh -f

## FIXME: conrod_bearing3_1.mbdyn does not pass the patched testsuite yet

MBDYN_EXEC="${MBDYN_EXEC:-mbdyn}"

exec ${MBDYN_EXEC} $*
