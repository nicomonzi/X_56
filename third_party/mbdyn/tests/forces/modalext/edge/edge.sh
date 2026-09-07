#!/bin/sh

MBDYN_EXEC=${MBDYN_EXEC:-mbdyn}

which ${MBDYN_EXEC}

RC=$?
if test $RC != 0 ; then
	echo "mbdyn is not in PATH"
	exit 1
fi

which test_modalext_edge
RC=$?
if test $RC != 0 ; then
	echo "test_modalext_edge is not in PATH"
	exit 1
fi

rm -f mflag.dat mdata.dat rflag.dat rdata.dat

echo "executing test_modalext_edge..."
test_modalext_edge -c 4 -s 1 \
	-M 4 -m flag=mflag.dat -m data=mdata.dat \
	-r flag=rflag.dat -r data=rdata.dat \
	> trashme.test1 2>&1 &
T1PID=$!

for i in 1 2 3 4 5 ; do
	if test -f mflag.dat \
		-a -f mdata.dat \
		-a -f rflag.dat \
		-a -f rdata.dat ; \
	then
		break;
	fi
	sleep $i
done


echo "executing mbdyn..."
${MBDYN_EXEC} ${GTEST_MBDYN_ARGS} -f edge -o trashme > trashme.mbdyn 2>&1 &
MBPID=$!

wait $MBPID $T1PID
echo "... end of simulation"

