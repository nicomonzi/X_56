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
# Copyright (C) 2024(-2024) all rights reserved.

# The copyright of this code is transferred
# to Pierangelo Masarati and Paolo Mantegazza
# for use in the software MBDyn as described
# in the GNU Public License version 2.1

function status = octave_pkg_testsuite_exec(idx, test_data)
  status = struct();

  try
    printf("%d:%s\n", idx, test_data.pkg_functions{idx});

    sigterm_dumps_octave_core(false);

    status.pkg_index = test_data.pkg_index(idx);
    status.pkg_function = test_data.pkg_functions{idx};

    output_dir = fullfile(test_data.octave_pkg_test_dir, sprintf('%d', idx));

    ## NOTE: mboct-mbdyn-pkg will fill use the format "%03d", just in case that multiple instances of MBDyn will be started from the same test function.

    junit_xml_report_file = fullfile(output_dir, "junit_xml_report_octave_%03d.xml");
    mbdyn_solver_command = sprintf("%s %s --gtest_output=xml:%s", test_data.mbdyn_exec, test_data.mbdyn_args_add, junit_xml_report_file);
    mbdyn_solver_command_var_name = "MBOCT_MBDYN_PKG_MBDYN_SOLVER_COMMAND";

    printf("%s=%s\n", mbdyn_solver_command_var_name, mbdyn_solver_command);

    putenv(mbdyn_solver_command_var_name, mbdyn_solver_command);

    if (test_data.enable_profile)
      profile('off');
      profile('clear');
      profile('on');
    endif

    printf("Run test %d: \"%s\" ...\n", idx, test_data.pkg_functions{idx});

    start = tic();

    [status.test.N, status.test.NMAX] = test(test_data.pkg_functions{idx}, test_data.test_args{:});

    status.test.elapsed_time = toc(start);

    printf("Exit status of test %d: \"%s\"\n", idx, test_data.pkg_functions{idx});
    printf("N=%d\n", status.test.N);
    printf("NMAX=%d\n", status.test.NMAX);

    if (test_data.enable_profile)
      profile('off');

      status.prof = profile('info');

      profile('clear');
    endif
  catch
    err = lasterror();
    error_report(err);
    rethrow(err);
  end_try_catch
endfunction
