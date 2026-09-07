#!/bin/sh -f

cd $(dirname $0)
exec awk -v fname=sze36_loadinc_m -v loadinc=1 -v mesh=40x40 -f sze36.awk
