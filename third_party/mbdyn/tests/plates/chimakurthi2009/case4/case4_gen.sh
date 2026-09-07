#!/bin/sh -f

exec awk -v mesh=16x6 -f case4.awk
