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

# AUTHOR: Reinhard Resch <mbdyn-user@a1.net>
# Copyright (C) 2024(-2025) all rights reserved.

# The copyright of this code is transferred
# to Pierangelo Masarati and Paolo Mantegazza
# for use in the software MBDyn as described
# in the GNU Public License version 2.1

function gtest_octave_cli_args = octave_pkg_testsuite_hook(idx, flags, test_data, pid, status)
  gtest_octave_cli_args = {};

  switch (flags)
    case "pre:spawn"
      printf("%d:%s:%s:%s\n", idx, flags, test_data.pkg_name{test_data.pkg_index(idx)}, test_data.pkg_functions{idx});
      output_dir = fullfile(test_data.octave_pkg_test_dir, sprintf('%03d', idx));

      [status, msg] = mkdir(output_dir);

      if (~status)
        error("failed to create directory \"%s\"", output_dir);
      endif

      putenv("TMPDIR", output_dir);

      gtest_octave_cli_args = {"--test-suite-name", ...
                               test_data.pkg_name{test_data.pkg_index(idx)}, ...
                               "--test-name", ...
                               test_data.pkg_functions{idx}};
    case "post:spawn"
      printf("%d:%s:%s:%s:pid=%d:status=%d\n", idx, flags, test_data.pkg_name{test_data.pkg_index(idx)}, test_data.pkg_functions{idx}, pid, status);
  endswitch
endfunction
