try
  clear all;
  close all;
  
  pkg_prefix = getenv("OCT_PKG_INSTALL_PREFIX");

  if (~isempty(pkg_prefix))
    pkg("local_list", fullfile(pkg_prefix, "octave_packages"));
  endif
  
  pkg load mboct-mbdyn-pkg;

  ad = [false, true];
  static = [false];

  max_diff = 0;

  for l=1:numel(static)
    fd = -1;
    unwind_protect
      fd = fopen("srscm_static.con", "wt");
      if (fd == -1)
        error("failed to open file");
      endif
      if (static(l))
        fputs(fd, "model: static;\n");
      else
        fputs(fd, "## model: dynamic;\n");
      endif
    unwind_protect_cleanup
      if (fd ~= -1)
        fclose(fd);
      endif
    end_unwind_protect
    for k=1:numel(ad)
      fd = -1;
      unwind_protect
        fd = fopen("srscm_ad.con", "wt");
        if (fd == -1)
          error("failed to open file");
        endif
        if (ad(k))
          fputs(fd, "use automatic differentiation;\n");
        else
          fputs(fd, "## use automatic differentiation;\n");
        endif
      unwind_protect_cleanup
        if (fd ~= -1)
          fclose(fd);
        endif
      end_unwind_protect
      opts.verbose = false;
      opts.output_file = "stage1";
      opts.logfile = "stage1.stdout";
      mbdyn_solver_run("srscm_stage1.mbd", opts);
      opts.output_file = "stage2";
      opts.logfile = "stage2.stdout";
      mbdyn_solver_run("srscm_stage2.mbd", opts);

      [t1, traj1, def1, vel1, acc1, node_id1] = mbdyn_post_load_output_struct("stage1");
      [t2, traj2, def2, vel2, acc2, node_id2] = mbdyn_post_load_output_struct("stage2");

      curr_diff = 0;

      for j=1:numel(traj1)
        x2 = traj2{j};
        x1 = zeros(size(x2));

        for i=1:columns(x1)
          x1(:, i) = interp1(t1, traj1{j}(:, i), t2);
        endfor

        curr_diff = max(curr_diff, max(max(abs(x1 - x2))) / max(1,max(max(abs(x1)))));
      endfor

      max_diff = max(max_diff, curr_diff);
      fprintf(stderr, "[%d, %d]: max difference: %e\n", l, k, curr_diff);
    endfor
  endfor
  tol = 1e-10;
  assert(max_diff < tol);
catch
  gtest_error = lasterror();
  gtest_fail(gtest_error, __FILE__);
  rethrow(gtest_error);
end_try_catch
