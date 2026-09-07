# MBDyn (C) is a multibody analysis code.
# http://www.mbdyn.org
#
# Copyright (C) 1996-2023
#
# Pierangelo Masarati	<pierangelo.masarati@polimi.it>
# Paolo Mantegazza	<paolo.mantegazza@polimi.it>
#
# Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
# via La Masa, 34 - 20156 Milano, Italy
# http://www.aero.polimi.it
#
# Changing this copyright notice is forbidden.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation (version 2 of the License).
#
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#

# AUTHOR: Reinhard Resch <mbdyn-user@a1.net>
#        Copyright (C) 2023(-2023) all rights reserved.
#
#        The copyright of this code is transferred
#        to Pierangelo Masarati and Paolo Mantegazza
#        for use in the software MBDyn as described
#        in the GNU Public License version 2.1

# parse_test_suite_status.awk: helper to parse the output from Octave's __run_test_suite__ function

BEGINFILE {
    failed = 0;
    passed = 0;
}

/^  PASS\>/ {
    ## Output from Octave's function "__run_test_suite__"
    passed += strtonum($2);
}

/^  FAIL\>/ {
    ## Output from Octave's function "__run_test_suite__"
    failed += strtonum($2);
}

/^PASSES [0-9]+ out of [0-9]+ test$/ {
    ## Output from Octave's function "test"
    passes = strtonum($2);
    total = strtonum($5);
    passed += passes;
    failed += total - passes;
}

/^FAILED\>/ {
    ## Output from Octave's function "test"
    failed += strtonum($2);
}

$2~/\<PASSED\>/ && $4~/^[0-9]+/ && $5~/\<test\.$/ {
    ## Output from gtest
    passed += strtonum($4);
}

$2~/\<FAILED\>/ && $4~/^[0-9]/ && $5~/\<test\>/ {
    ## Output from gtest
    failed += strtonum($4);
}

/\?\?\?\?\? .+has no tests available$|\?\?\?\?\? .+source code with tests for dynamically linked function not found$/ {
    ## Output from Octave's function "test"
    ## This is not considered as a failure.
    ++passed;
}

/^!!!!! test failed$/ {
    ++failed;
}

$1~/\<testsuites\>/ && $2~/\<tests="[0-9]+"/ && $3~/\<failures="[0-9]+"/ && $5~/\<errors="[0-9]+"/ {
    ntests = split($2, tests, "\"");
    nfailures = split($3, failures, "\"");
    nerrors = split($5, errors, "\"");

    if (ntests == 3 && nfailures == 3 && nerrors == 3) {
        passed += tests[2] - failures[2] - errors[2];
        failed += failures[2] + errors[2];
    }
}

ENDFILE {
    if (failed > 0) {
        exit 1;
    }

    ## Use zero terminated strings.
    ## So, we can use a similar command to clean up all the output which does not indicate a failure:
    ## find ${TMPDIR} -name 'fntests.out' -print0 | awk -f parse_test_suite_status.awk | xargs -0 rm
    printf("%s\0", FILENAME);
}
