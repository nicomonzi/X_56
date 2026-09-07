#!/bin/sh -f

cd $(dirname $0)
exec awk -v loadinc=0 -v mesh=40x40 -f sze36.awk
