#! /bin/sh

MBDYN_EXEC=${MBDYN_EXEC:-mbdyn}

${MBDYN_EXEC} ${GTEST_MBDYN_ARGS} -f planar_3link_ergo.mbd -o z -s

awk '$1==12 {x1=$7} $1==13 {x2=$7} $1==14 {print x1, x2, $7}' z.jnt > z_c.dat
awk '$1==12 {x1=$19} $1==13 {x2=$19} $1==14 {print x1, x2, $19}' z.jnt > z_t.dat
awk '$1==12 {x1=$25} $1==13 {x2=$25} $1==14 {print x1, x2, $25}' z.jnt > z_w.dat

${MBDYN_EXEC} ${GTEST_MBDYN_ARGS} -f planar_3link_ergo_check.mbd -o c -s
