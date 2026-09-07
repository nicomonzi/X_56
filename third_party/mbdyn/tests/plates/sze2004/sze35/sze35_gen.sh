#!/bin/sh -f

cd $(dirname $0)
exec awk -v mesh=16x24 -f sze35.awk
