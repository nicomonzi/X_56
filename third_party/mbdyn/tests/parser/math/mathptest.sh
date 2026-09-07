#!/bin/bash

FS='[[:space:]]*@[[:space:]]*'
OUT="out.txt"

# ee
#MBDYN="/home/masarati/Lavoro/mbdyn/GSoC/2015/Ankit/mbdyn-HEAD-28jun.X/mbdyn/mbdyn"
MBDYN="../../../mbdyn/mbdyn"
SEP=","
QUOTES="\""

# original
#MBDYN="/usr/local/mbdyn/bin/mbdyn"
#SEP=":"
#QUOTES=""

awk -F "$FS" -v MBDYN="$MBDYN" -f mathptest.awk mathptest.dat | bash > $OUT
sed -e 's/^.*=[[:space:]]\+\([^'$SEP']\+\)'$QUOTES$SEP'[[:space:]"]\+\([^"]*\)["]*/\1 \2 \0/' $OUT \
        | awk 'BEGIN {cnt=0; ok=0; err=0} {cnt++; if ($1 != $2) {out = "ERR"; err++} else {out = "ok"; ok++} print $0, out} END {print "ok=" ok "/" cnt, "ERR=" err "/" cnt}'

rm -f $OUT
