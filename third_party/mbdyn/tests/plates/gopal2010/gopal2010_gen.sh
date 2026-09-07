#!/bin/sh -f

OCTAVE_EXEC=${OCTAVE_EXEC:-octave}

exec ${OCTAVE_EXEC} ${GTEST_OCTAVE_ARGS} --eval `printf 'cd("%s");gopal2010;' "$(dirname $0)"`
