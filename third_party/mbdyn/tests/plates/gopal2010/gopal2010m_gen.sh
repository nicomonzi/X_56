#!/bin/sh -f

OCTAVE_EXEC=${OCTAVE_EXEC:-octave}

exec ${OCTAVE_EXEC} ${GTEST_OCTAVE_ARGS} --eval `printf 'cd("%s");gopal2010m;' "$(dirname $0)"`
