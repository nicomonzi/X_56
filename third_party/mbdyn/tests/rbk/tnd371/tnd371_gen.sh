#!/bin/sh -f

OCTAVE_EXEC=${OCTAVE_EXEC:-octave}

exec ${OCTAVE_EXEC} ${GTEST_OCTAVE_ARGS} --eval `printf 'cd("%s");tnd371;' "$(dirname $0)"`
