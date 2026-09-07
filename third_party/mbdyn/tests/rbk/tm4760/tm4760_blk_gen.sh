#!/bin/sh -f

OCTAVE_EXEC=${OCTAVE_EXEC:-octave}

exec ${OCTAVE_EXEC} ${GTEST_OCTAVE_ARGS} --eval `printf 'cd("%s");tm4760_blk;' "$(dirname $0)"`
