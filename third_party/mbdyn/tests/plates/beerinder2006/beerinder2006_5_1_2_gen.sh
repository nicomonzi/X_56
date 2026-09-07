#!/bin/sh -f

OCTAVE_EXEC=${OCTAVE_EXEC:-octave}

exec ${OCTAVE_EXEC} ${GTEST_OCTAVE_ARGS} beerinder2006_5_1_2.m $*
