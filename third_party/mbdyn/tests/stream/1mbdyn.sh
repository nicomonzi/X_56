#!/bin/sh

PYTHON_EXEC="${PYTHON_EXEC:-python3}"
MBDYN_EXEC=${MBDYN_EXEC:-mbdyn}

echo "starting body"
nohup ${MBDYN_EXEC} -o 1b body ${GTEST_MBDYN_ARGS} > 1b.txt 2>&1 &
PIDS=$!

echo -n "waiting for body to start... "
sleep 1
echo "done"

echo "starting spring"
${PYTHON_EXEC} spring.py $*

echo -n "waiting for body to complete... "
wait $PIDS
echo "done"

echo "when done: rm 1b.* *.echo"

