#! /bin/sh

MBDYN_EXEC=${MBDYN_EXEC:-mbdyn}

${MBDYN_EXEC} ${GTEST_MBDYN_ARGS} -N1 -o xxx model 2>&1 \
	| grep expected \
	| sed -e 's/[^:]*:\(.*\) expected: \(.*=\)\?\([-0-9.]\+\), actual: \([-0-9.]\+\)$/\1: \3 \4/'
rm -f xxx.*

