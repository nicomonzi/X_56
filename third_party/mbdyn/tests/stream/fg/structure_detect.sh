#!/bin/sh

### uncomment desired version
### it assumes that flightgear-$FGVERSION and simgear-$FGVERSION are in the parent dir
#FGVERSION="2.4.0"
FGVERSION="3.2.0"
#FGVERSION="3.4.0"

### uncomment to "explode" all data
#FGVER="0`echo $FGVERSION | sed -e 's/\./0/g'`"

### uncomment to explode listed data only
FGVER=0
CTRLS="version aileron elevator rudder num_engines throttle"
FDM="version longitude latitude altitude agl phi theta psi alpha"

### do not edit below this line

CTRLS_HCC="ctrls.hcc"
echo "std::cout << \"CTRLS:\" << std::endl;" > "${CTRLS_HCC}"
for c in $CTRLS ; do
	echo "MYFG_OFFSETOF(FGNetCtrls, ctrls, ${c}) << std::endl;" >> "${CTRLS_HCC}"
done

FDM_HCC="fdm.hcc"
echo "std::cout << \"FDM:\" << std::endl;" > "${FDM_HCC}"
for f in $FDM ; do
	echo "MYFG_OFFSETOF(FGNetFDM, fdm, ${f}) << std::endl;" >> "${FDM_HCC}"
done

g++ \
	-Wall \
	-D HAVE_CXXABI_H=1 \
	-D FGVERSION="\"${FGVERSION}\"" \
	-D FGVER="${FGVER}" \
	-I ../simgear-"${FGVERSION}"/ \
	-I ../flightgear-"${FGVERSION}"/src/ \
	-o structure_detect \
	structure_detect.cc

./structure_detect

