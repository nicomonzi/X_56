#!/bin/sh -f

PYTHON_EXEC="${PYTHON_EXEC:-python3}"

exec ${PYTHON_EXEC} gopal2010_static.py $*
