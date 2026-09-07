#!/bin/sh -f

## Allow all our emacs and vim users to update their "TAGS" and "BROWSE" files conveniently by means of etags(1) and ebrowse(1).

## MBDyn (C) is a multibody analysis code.
## http://www.mbdyn.org
##
## Copyright (C) 1996-2023
##
## Pierangelo Masarati     <pierangelo.masarati@polimi.it>
## Paolo Mantegazza        <paolo.mantegazza@polimi.it>
##
## Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
## via La Masa, 34 - 20156 Milano, Italy
## http://www.aero.polimi.it
##
## Changing this copyright notice is forbidden.
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation (version 2 of the License).
##
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
##
## ----------------------------------------------------------------

program_name="$0"

program_dir=$(realpath $(dirname "${program_name}"))

if ! test -f "${program_dir}/emacs-vi-update-tags.sh"; then
    program_dir=$(realpath $(which "${program_name}"))
fi

if ! test -d "${program_dir}"; then
    program_dir="."
fi

mbdyn_locate_all_cxx_source_files()
{
    cd "$1"

    find . '(' -type f -and -not -path './.git/*' -and -not -path './tests/*' -and -not -path './manual/*' ')' -print0 | \
        xargs -0 file | \
        awk -v FPAT='[^:]+' '/[[:space:]]+C\+\+[[:space:]]+(\<source\>|\<header\>)/{printf("%s\0", $1);}'
}

mbdyn_locate_all_cxx_source_files "${program_dir}" | xargs -0 etags
mbdyn_locate_all_cxx_source_files "${program_dir}" | xargs -0 ebrowse
