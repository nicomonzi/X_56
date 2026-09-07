#!/usr/bin/env gtest-octave-cli

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

clear all;

total.N = 0;
total.NMAX = 0;
total.NREGRESSION = 0;
total_pkg = total([]);
test_data.pkg_name = {};
test_data.octave_pkg_test_dir = "";
test_data.mbdyn_exec = "mbdyn";
test_data.mbdyn_args_add ="-CGF";
test_data.octave_pkg_prefix = getenv("OCT_PKG_INSTALL_PREFIX");
thread_data.number_of_threads = 1;

try
  if (~isempty(test_data.octave_pkg_prefix))
    fprintf(stderr, "\npkg(\"local_list\", \"%s\");\n\n", test_data.octave_pkg_prefix);
    pkg("local_list", fullfile(test_data.octave_pkg_prefix, "octave_packages"));
  endif

  pkg load mboct-octave-pkg;

  [prog_dir, prog_name, prog_ext] = fileparts(__FILE__);

  addpath(prog_dir);

  args = argv();

  output_file = "";

  idx = int32(0);

  opts.number_of_processors = int32(1);
  opts.reuse_subprocess = false;
  opts.verbose = false;
  opts.user_hook_func = @octave_pkg_testsuite_hook;

  ## FIXME: We should use semaphores instead of polling!
  ## FIXME: However this appears to be less robust,
  ## FIXME: just in case that one of our jobs terminates without releasing the semaphore.
  opts.waitpid_polling_period = 100e-3; ## Higher values will reduce the CPU time which is wasted inside the main loop!

  while (++idx < numel(args))
    switch(args{idx})
      case {"--octave-pkg-list", "-p"}
        pkg_list = strsplit(args{++idx});
        add_pkg_name = cell(1, numel(pkg_list));
        num_pkg = int32(0);
        for i=1:numel(pkg_list)
          pkg_flags = strsplit(pkg_list{i}, ":");
          switch (pkg_flags{4})
            case "yes"
              add_pkg_name{++num_pkg} = pkg_flags{1};
          endswitch
        endfor
        test_data.pkg_name = {test_data.pkg_name{:}, add_pkg_name{1:num_pkg}};
      case "--octave-exec"
        opts.octave_exec = args{++idx};
      case "--octave-args-append"
        opts.octave_args_append = args{++idx};
      case "--octave-pkg-test-dir"
        test_data.octave_pkg_test_dir = make_absolute_filename(args{++idx});
      case "--verbose"
        opts.verbose = true;
      case "--mbdyn-exec"
        test_data.mbdyn_exec = args{++idx};
      case "--mbdyn-args-add"
        test_data.mbdyn_args_add = args{++idx};
      case {"--tasks", "-t"}
        [opts.number_of_processors, cnt, msg] = sscanf(args{++idx}, "%d", "C");

        if (cnt ~= 1)
          error("invalid argument %s \"%s\": %s", args{idx - 1}, args{idx}, msg);
        endif
      case "--threads"
        [thread_data.number_of_threads, cnt, msg] = sscanf(args{++idx}, "%d", "C");

        if (cnt ~= 1)
          error("invalid argument %s \"%s\": %s", args{idx - 1}, args{idx}, msg);
        endif
      case {"--help", "-h"}
        fprintf(stderr, "octave_pkg_testsuite.m\n");
        fprintf(stderr, "\t\t\t--package-name|-p <test_data.pkg_name>\n");
        return;
      otherwise
        error("unknown argument \"%s\"", args{idx});
    endswitch
  endwhile

  if (isempty(test_data.pkg_name))
    error("missing argument --octave-pkg-list <PKG_NAME>:<BUILD>:<BRANCH>:<RUN>:<LIMIT>");
  endif

  if (isempty(test_data.octave_pkg_test_dir))
    error("missing argument --octave-pkg-test-dir <DIR_NAME>");
  endif

  sigterm_dumps_octave_core(false);

  test_data.pkg_functions = {};
  test_data.pkg_index = zeros(0, 0, "int32");

  for j=1:numel(test_data.pkg_name)
    pkg("load", test_data.pkg_name{j});
    pkg_list = pkg('list', '-verbose', test_data.pkg_name{j});
    pkg_files = {dir(fullfile(pkg_list{1}.dir, '*.m')).name, dir(fullfile(pkg_list{1}.dir, '*.tst')).name};

    test_data.pkg_index(end + (1:numel(pkg_files))) = repmat(int32(j), 1, numel(pkg_files));
    test_data.pkg_functions = {test_data.pkg_functions{:}, pkg_files{:}};
  endfor

  test_data.enable_profile = true;
  test_data.test_args = {"normal"};

  opts.number_of_parameters = numel(test_data.pkg_functions);

  opts.gtest_output_junit_xml = fullfile(test_data.octave_pkg_test_dir, '%03d', 'junit_xml_report_octave_assert.xml');
  opts.redirect_stdout = fullfile(test_data.octave_pkg_test_dir, '%03d', 'fntests.out');

  putenv("MBD_NUM_THREADS", sprintf("%d", thread_data.number_of_threads));

  status = run_parallel(opts, @octave_pkg_testsuite_exec, test_data);

  total_pkg = repmat(total, 1, numel(test_data.pkg_name));

  for i=1:numel(status)
    NREGRESSION = status{i}.test.NMAX - status{i}.test.N;
    total_pkg(status{i}.pkg_index).N += status{i}.test.N;
    total_pkg(status{i}.pkg_index).NMAX += status{i}.test.NMAX;
    total_pkg(status{i}.pkg_index).NREGRESSION += NREGRESSION;

    total.N += status{i}.test.N;
    total.NMAX += status{i}.test.NMAX;
    total.NREGRESSION += NREGRESSION;

    printf("%3d: test(\"%s\":\"%s\"): %2d/%2d passed, %2d/%2d failed, elapsed time: %.2fs\n", i, test_data.pkg_name{status{i}.pkg_index}, status{i}.pkg_function, status{i}.test.N, status{i}.test.NMAX, NREGRESSION, status{i}.test.NMAX, status{i}.test.elapsed_time);
  endfor
catch
  gtest_error = lasterror();
  gtest_fail(gtest_error, __FILE__);
  printf("%s FAILED\n", __FILE__);
  rethrow(gtest_error);
end_try_catch

if (total.NREGRESSION > 0)
  rc = 1;
else
  rc = 0;
endif

printf("\nTest summary per package:\n\n");
for i=1:numel(total_pkg)
  printf("%2d: Test summary \"%s\":\n", i, test_data.pkg_name{i});
  printf("%3d/%3d tests passed\n", total_pkg(i).N, total_pkg(i).NMAX);
  printf("%3d/%3d tests failed\n\n", total_pkg(i).NREGRESSION, total_pkg(i).NMAX);
endfor

printf("\nTest summary for all packages:\n");
printf("%3d/%3d tests passed\n", total.N, total.NMAX);
printf("%3d/%3d tests failed\n\n", total.NREGRESSION, total.NMAX);

exit(rc);
