#! /bin/sh
# $Header: /var/cvs/mbdyn/mbdyn/mbdyn-1.0/tests/plates/agard/agard_mapping.sh,v 1.9 2011/10/19 09:57:30 masarati Exp $
#
# 1) run
#	awk -f agard_solid.awk > agard_solid.dat
#    to create the external grid
#
# 2) edit file "agard_mapping"
#	choose the desired mapping (include: "agard_mapping_*.elm":)
# 	uncomment the "echo" line in the "external structural mapping" element
# 
# 3)  run 
#       mbdyn -f agard_mapping
#
# 5) run with "octave" 
# 	create_mls_interface("agard_mb.dat")    
#
# 6) edit file "agard_mapping_*.elm"
# 	comment the "echo" line
#
# 7) run this script

MBDYN_EXEC=${MBDYN_EXEC:-mbdyn}

which ${MBDYN_EXEC} >/dev/null
RC=$?
if test $RC != 0 ; then
	echo "mbdyn is not in PATH"
	exit 1
fi

which test_strext_socket >/dev/null
RC=$?
if test $RC != 0 ; then
	echo "test_strext_socket is not in PATH"
	exit 1
fi

FORCES=
if test -f mapped_forces.dat ; then
	FORCES="-i mapped_forces.dat"
fi

OUTPUT=trashme

echo "executing mbdyn..."
${MBDYN_EXEC} ${GTEST_MBDYN_ARGS} -f agard_mapping -o $OUTPUT > $OUTPUT.mbdyn 2>&1 &
MBPID=$!

for i in 1 2 3 4 5 ; do
	if test -S /tmp/mbdyn_agard.sock ; then
		break;
	fi
	sleep $i
done

echo "executing test_strext_socket for socket #1..."
test_strext_socket \
	-v \
	-c 4 \
	-N 12025 \
	-n \
	$FORCES \
	-s 0 \
	-r \
	-o $OUTPUT.mapped \
	-x \
	-H local:///tmp/mbdyn_agard.sock \
	>$OUTPUT.test1 2>&1 &
T1PID=$!

wait $MBPID $T1PID
echo "... end of simulation"

# 8) result visualization
#       execute gnuplot
#
#       sp"<awk '$1==\"POS\" {i++;p=1; getline} $1==\"VEL\" {p=0} i==<STEP*ITER> && p==1 {print $0}' $OUTPUT.mapped", \
#         "<awk '$1==1 {i++} i==(STEP + 1) {print $0}' $OUTPUT.mov"u 2:3:4
#
# 9) safely remove all files $OUTPUT.*
