#!/bin/sh

MBDYN_EXEC=${MBDYN_EXEC:-mbdyn}

echo "starting body"
nohup ${MBDYN_EXEC} -o 2b body ${GTEST_MBDYN_ARGS} > 2b.txt 2>&1 &
PIDS=$!

echo -n "waiting for body to start... "
sleep 1
echo "done"

echo "starting spring"
nohup ${MBDYN_EXEC} -o 2s spring ${GTEST_MBDYN_ARGS} > 2s.txt 2>&1 &
PIDS="$PIDS $!"

echo -n "waiting for body and spring to complete... "
wait $PIDS
echo "done"

echo "when done: rm 2b.* 2s.* *.echo"
