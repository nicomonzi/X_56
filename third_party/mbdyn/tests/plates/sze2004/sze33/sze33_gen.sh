#!/bin/sh -f

cd $(dirname $0)
exec awk -v mesh=6x30 -f sze33.awk
