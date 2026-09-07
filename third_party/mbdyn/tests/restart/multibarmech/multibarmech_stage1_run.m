% MBDyn (C) is a multibody analysis code.
% http://www.mbdyn.org
%
% Copyright (C) 1996-2017
%
% Pierangelo Masarati	<masarati@aero.polimi.it>
% Paolo Mantegazza	<mantegazza@aero.polimi.it>
%
% Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
% via La Masa, 34 - 20156 Milano, Italy
% http://www.aero.polimi.it
%
% Changing this copyright notice is forbidden.
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation (version 2 of the License).
%
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

try
  clear all;
  close all;

  pkg_prefix = getenv("OCT_PKG_INSTALL_PREFIX");

  if (~isempty(pkg_prefix))
    pkg("local_list", fullfile(pkg_prefix, "octave_packages"));
  endif

  N_Y_LOOPS = 5;
  N_Z_LOOPS = 5;

  Z_OFFSET = 100;
  Y_OFFSET = 10000;

  fd = fopen('multibarmech_m.set', 'w');

  fprintf(fd, '# MBDyn (C) is a multibody analysis code.\n');
  fprintf(fd, '# http://www.mbdyn.org\n');
  fprintf(fd, '# \n');
  fprintf(fd, '# Copyright (C) 1996-2017\n');
  fprintf(fd, '# \n');
  fprintf(fd, '# Pierangelo Masarati	<masarati@aero.polimi.it>\n');
  fprintf(fd, '# Paolo Mantegazza	<mantegazza@aero.polimi.it>\n');
  fprintf(fd, '# \n');
  fprintf(fd, '# Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano\n');
  fprintf(fd, '# via La Masa, 34 - 20156 Milano, Italy\n');
  fprintf(fd, '# http://www.aero.polimi.it\n');
  fprintf(fd, '# \n');
  fprintf(fd, '# Changing this copyright notice is forbidden.\n');
  fprintf(fd, '# \n');
  fprintf(fd, '# This program is free software; you can redistribute it and/or modify\n');
  fprintf(fd, '# it under the terms of the GNU General Public License as published by\n');
  fprintf(fd, '# the Free Software Foundation (version 2 of the License).\n');
  fprintf(fd, '# \n');
  fprintf(fd, '# \n');
  fprintf(fd, '# This program is distributed in the hope that it will be useful,\n');
  fprintf(fd, '# but WITHOUT ANY WARRANTY; without even the implied warranty of\n');
  fprintf(fd, '# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n');
  fprintf(fd, '# GNU General Public License for more details.\n');
  fprintf(fd, '# \n');
  fprintf(fd, '# You should have received a copy of the GNU General Public License\n');
  fprintf(fd, '# along with this program; if not, write to the Free Software\n');
  fprintf(fd, '# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA\n');
  fprintf(fd, '\n');

  fprintf(fd, 'set: const integer N_NODES = %d;\n', (2*N_Y_LOOPS + 1)*N_Z_LOOPS);
  fprintf(fd, 'set: const integer N_BODIES = %d;\n', (2*N_Y_LOOPS + 1)*N_Z_LOOPS);
  fprintf(fd, 'set: const integer N_JOINTS = %d;\n', (3*N_Y_LOOPS + 1)*N_Z_LOOPS);
  fprintf(fd, '\n');
  fprintf(fd, 'set: const real L = 1.; # m\n');
  fprintf(fd, 'set: const real M = 1.; # kg\n');
  fprintf(fd, 'set: real DOTPHI = 2*pi/3; # rad/s\n');
  fprintf(fd, '\n');
  fprintf(fd, 'set: real T1 = 10; # s\n');

  fclose(fd);



  fd = fopen('multibarmech_m.nod', 'w');

  fprintf(fd, '# MBDyn (C) is a multibody analysis code.\n');
  fprintf(fd, '# http://www.mbdyn.org\n');
  fprintf(fd, '# \n');
  fprintf(fd, '# Copyright (C) 1996-2017\n');
  fprintf(fd, '# \n');
  fprintf(fd, '# Pierangelo Masarati	<masarati@aero.polimi.it>\n');
  fprintf(fd, '# Paolo Mantegazza	<mantegazza@aero.polimi.it>\n');
  fprintf(fd, '# \n');
  fprintf(fd, '# Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano\n');
  fprintf(fd, '# via La Masa, 34 - 20156 Milano, Italy\n');
  fprintf(fd, '# http://www.aero.polimi.it\n');
  fprintf(fd, '# \n');
  fprintf(fd, '# Changing this copyright notice is forbidden.\n');
  fprintf(fd, '# \n');
  fprintf(fd, '# This program is free software; you can redistribute it and/or modify\n');
  fprintf(fd, '# it under the terms of the GNU General Public License as published by\n');
  fprintf(fd, '# the Free Software Foundation (version 2 of the License).\n');
  fprintf(fd, '# \n');
  fprintf(fd, '# \n');
  fprintf(fd, '# This program is distributed in the hope that it will be useful,\n');
  fprintf(fd, '# but WITHOUT ANY WARRANTY; without even the implied warranty of\n');
  fprintf(fd, '# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n');
  fprintf(fd, '# GNU General Public License for more details.\n');
  fprintf(fd, '# \n');
  fprintf(fd, '# You should have received a copy of the GNU General Public License\n');
  fprintf(fd, '# along with this program; if not, write to the Free Software\n');
  fprintf(fd, '# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA\n');
  fprintf(fd, '\n');

  for iz = 0:N_Z_LOOPS - 1,
    fprintf(fd, '# loop %d vertical bodies\n', iz + 1);
    for iy = 0:N_Y_LOOPS,
      fprintf(fd, 'structural: %d, dynamic,\n', Z_OFFSET*iz + iy);
      fprintf(fd, '\treference, global, 0., %d*L, -%d*L - L/2,\n', iy, iz);
      fprintf(fd, '\treference, global, eye,\n');
    %if (iz == 0),
    %	fprintf(fd, '\treference, global, 0., DOTPHI*L/2, 0.,\n');
    %	fprintf(fd, '\treference, global, DOTPHI, 0., 0.;\n');
    %else
    %	fprintf(fd, '\treference, global, 0., DOTPHI*L, 0.,\n');
    %	fprintf(fd, '\treference, global, null;\n');
    %end
      fprintf(fd, '\treference, global, 0., DOTPHI*(%d*L + L/2), 0.,\n', iz);
      fprintf(fd, '\treference, global, DOTPHI, 0., 0.;\n');
      fprintf(fd, '\n');
    end

    fprintf(fd, '# loop %d horizontal bodies\n', iz + 1);
    for iy = 1:N_Y_LOOPS,
      fprintf(fd, 'structural: %d, dynamic,\n', Y_OFFSET + Z_OFFSET*iz + iy);
      fprintf(fd, '\treference, global, 0., %d*L + L/2, -%d*L,\n', iy - 1, iz + 1);
      fprintf(fd, '\treference, global, eye,\n');
            % fprintf(fd, '\treference, global, 0., DOTPHI*L, 0.,\n');
      fprintf(fd, '\treference, global, 0., DOTPHI*%d*L, 0.,\n', iz + 1);
      fprintf(fd, '\treference, global, null;\n');
      fprintf(fd, '\n');
    end
  end

  fclose(fd);



  fd = fopen('multibarmech_m.elm', 'w');

  fprintf(fd, '# MBDyn (C) is a multibody analysis code.\n');
  fprintf(fd, '# http://www.mbdyn.org\n');
  fprintf(fd, '# \n');
  fprintf(fd, '# Copyright (C) 1996-2017\n');
  fprintf(fd, '# \n');
  fprintf(fd, '# Pierangelo Masarati	<masarati@aero.polimi.it>\n');
  fprintf(fd, '# Paolo Mantegazza	<mantegazza@aero.polimi.it>\n');
  fprintf(fd, '# \n');
  fprintf(fd, '# Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano\n');
  fprintf(fd, '# via La Masa, 34 - 20156 Milano, Italy\n');
  fprintf(fd, '# http://www.aero.polimi.it\n');
  fprintf(fd, '# \n');
  fprintf(fd, '# Changing this copyright notice is forbidden.\n');
  fprintf(fd, '# \n');
  fprintf(fd, '# This program is free software; you can redistribute it and/or modify\n');
  fprintf(fd, '# it under the terms of the GNU General Public License as published by\n');
  fprintf(fd, '# the Free Software Foundation (version 2 of the License).\n');
  fprintf(fd, '# \n');
  fprintf(fd, '# \n');
  fprintf(fd, '# This program is distributed in the hope that it will be useful,\n');
  fprintf(fd, '# but WITHOUT ANY WARRANTY; without even the implied warranty of\n');
  fprintf(fd, '# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n');
  fprintf(fd, '# GNU General Public License for more details.\n');
  fprintf(fd, '# \n');
  fprintf(fd, '# You should have received a copy of the GNU General Public License\n');
  fprintf(fd, '# along with this program; if not, write to the Free Software\n');
  fprintf(fd, '# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA\n');
  fprintf(fd, '\n');

  for iz = 0:N_Z_LOOPS - 1,
    fprintf(fd, '# loop %d vertical bodies\n', iz + 1);
    for iy = 0:N_Y_LOOPS,
      fprintf(fd, 'body: %d, %d,\n', Z_OFFSET*iz + iy, Z_OFFSET*iz + iy);
      fprintf(fd, '\tM,\n');
      fprintf(fd, '\tnull,\n');
      fprintf(fd, '\tdiag, M*L^2/12, M*L^2/12, 1.;\n');
      fprintf(fd, '\n');
    end

    fprintf(fd, '# loop %d horizontal bodies\n', iz + 1);
    for iy = 1:N_Y_LOOPS,
      fprintf(fd, 'body: %d, %d,\n', Y_OFFSET + Z_OFFSET*iz + iy, Y_OFFSET + Z_OFFSET*iz + iy);
      fprintf(fd, '\tM,\n');
      fprintf(fd, '\tnull,\n');
      fprintf(fd, '\tdiag, M*L^2/12, 1., M*L^2/12;\n');
      fprintf(fd, '\n');
    end
  end

  for iz = 0:N_Z_LOOPS - 1,
    fprintf(fd, '# loop %d of joints\n', iz + 1);

                                %	|
                                %	o
                                %	|
    if (iz == 0),
      fprintf(fd, 'joint: %d, total pin joint,\n', 0);
      fprintf(fd, '\t%d,\n', 0);
      fprintf(fd, '\t\tposition, reference, node, 0., 0., L/2,\n');
      fprintf(fd, '\t\tposition orientation, reference, node, eye,\n');
      fprintf(fd, '\t\trotation orientation, reference, node, eye,\n');
      fprintf(fd, '\t# ground\n');
      fprintf(fd, '\t\tposition, reference, global, 0., %d*L, 0.,\n', 0);
      fprintf(fd, '\t\tposition orientation, reference, global, eye,\n');
      fprintf(fd, '\t\trotation orientation, reference, global, eye,\n');
      fprintf(fd, '\tposition constraint, 1, 1, 1, null,\n');
      fprintf(fd, '\torientation constraint, 0, 1, 1, null;\n');
    else
      fprintf(fd, 'joint: %d, total joint,\n', iz*Z_OFFSET + 0);
      fprintf(fd, '\t%d,\n', iz*Z_OFFSET + 0);
      fprintf(fd, '\t\tposition, reference, node, 0., 0., L/2,\n');
      fprintf(fd, '\t\tposition orientation, reference, node, eye,\n');
      fprintf(fd, '\t\trotation orientation, reference, node, eye,\n');
      fprintf(fd, '\t%d,\n', (iz - 1)*Z_OFFSET + 0);
      fprintf(fd, '\t\tposition, reference, node, 0., 0., -L/2,\n', 0);
      fprintf(fd, '\t\tposition orientation, reference, node, eye,\n');
      fprintf(fd, '\t\trotation orientation, reference, node, eye,\n');
      fprintf(fd, '\tposition constraint, 1, 1, 1, null,\n');
      fprintf(fd, '\torientation constraint, 0, 1, 1, null;\n');
    end
    fprintf(fd, '\n');

                                %	|
                                %	o
                                %	|
    if (iz == 0),
      for iy = 1:N_Y_LOOPS,
        fprintf(fd, 'joint: %d, total pin joint,\n', iy);
        fprintf(fd, '\t%d,\n', iy);
        fprintf(fd, '\t\tposition, reference, node, 0., 0., L/2,\n');
        fprintf(fd, '\t\tposition orientation, reference, node, eye,\n');
        fprintf(fd, '\t\trotation orientation, reference, node, eye,\n');
        fprintf(fd, '\t# ground\n');
        fprintf(fd, '\t\tposition, reference, global, 0., %d*L, 0.,\n', iy);
        fprintf(fd, '\t\tposition orientation, reference, global, eye,\n');
        fprintf(fd, '\t\trotation orientation, reference, global, eye,\n');
        fprintf(fd, '\tposition constraint, 0, 1, 1, null,\n');
        fprintf(fd, '\torientation constraint, 0, 0, 0, null;\n');
        fprintf(fd, '\n');
      end
    else
      for iy = 1:N_Y_LOOPS,
        fprintf(fd, 'joint: %d, total joint,\n', iz*Z_OFFSET + iy);
        fprintf(fd, '\t%d,\n', iz*Z_OFFSET + iy);
        fprintf(fd, '\t\tposition, reference, node, 0., 0., L/2,\n');
        fprintf(fd, '\t\tposition orientation, reference, node, eye,\n');
        fprintf(fd, '\t\trotation orientation, reference, node, eye,\n');
        fprintf(fd, '\t%d,\n', (iz - 1)*Z_OFFSET + iy);
        fprintf(fd, '\t\tposition, reference, node, 0., 0., -L/2,\n');
        fprintf(fd, '\t\tposition orientation, reference, node, eye,\n');
        fprintf(fd, '\t\trotation orientation, reference, node, eye,\n');
        fprintf(fd, '\tposition constraint, 0, 1, 1, null,\n');
        fprintf(fd, '\torientation constraint, 0, 0, 0, null;\n');
        fprintf(fd, '\n');
      end
    end

                                %	|
                                %	o--
    fprintf(fd, 'joint: %d, total joint,\n', Y_OFFSET + iz*Z_OFFSET + 0);
    fprintf(fd, '\t%d,\n', iz*Z_OFFSET + 0);
    fprintf(fd, '\t\tposition, reference, node, 0., 0., -L/2,\n');
    fprintf(fd, '\t\tposition orientation, reference, node, eye,\n');
    fprintf(fd, '\t\trotation orientation, reference, node, eye,\n');
    fprintf(fd, '\t%d,\n', Y_OFFSET + iz*Z_OFFSET + 1);
    fprintf(fd, '\t\tposition, reference, node, 0., -L/2, 0,\n');
    fprintf(fd, '\t\tposition orientation, reference, node, eye,\n');
    fprintf(fd, '\t\trotation orientation, reference, node, eye,\n');
    fprintf(fd, '\tposition constraint, 1, 1, 1, null,\n');
    fprintf(fd, '\torientation constraint, 0, 1, 1, null;\n');
    fprintf(fd, '\n');

    for iy = 1:N_Y_LOOPS - 1,
      fprintf(fd, 'joint: %d, total joint,\n', Y_OFFSET + iz*Z_OFFSET + iy);
      fprintf(fd, '\t%d,\n', iz*Z_OFFSET + iy);
      fprintf(fd, '\t\tposition, reference, node, 0., 0., -L/2,\n');
      fprintf(fd, '\t\tposition orientation, reference, node, eye,\n');
      fprintf(fd, '\t\trotation orientation, reference, node, eye,\n');
      fprintf(fd, '\t%d,\n', Y_OFFSET + iz*Z_OFFSET + iy + 1);
      fprintf(fd, '\t\tposition, reference, node, 0., -L/2, 0,\n');
      fprintf(fd, '\t\tposition orientation, reference, node, eye,\n');
      fprintf(fd, '\t\trotation orientation, reference, node, eye,\n');
      fprintf(fd, '\tposition constraint, 0, 1, 1, null,\n');
      fprintf(fd, '\torientation constraint, 0, 0, 0, null;\n');
      fprintf(fd, '\n');
    end

                                %	  |
                                %	--o
    for iy = 1:N_Y_LOOPS,
      fprintf(fd, 'joint: %d, total joint,\n', 2*Y_OFFSET + iz*Z_OFFSET + iy);
      fprintf(fd, '\t%d,\n', iz*Z_OFFSET + iy);
      fprintf(fd, '\t\tposition, reference, node, 0., 0., -L/2,\n');
      fprintf(fd, '\t\tposition orientation, reference, node, eye,\n');
      fprintf(fd, '\t\trotation orientation, reference, node, eye,\n');
      fprintf(fd, '\t%d,\n', Y_OFFSET + iz*Z_OFFSET + iy);
      fprintf(fd, '\t\tposition, reference, node, 0., L/2, 0,\n');
      fprintf(fd, '\t\tposition orientation, reference, node, eye,\n');
      fprintf(fd, '\t\trotation orientation, reference, node, eye,\n');
      fprintf(fd, '\tposition constraint, 0, 1, 1, null,\n');
      fprintf(fd, '\torientation constraint, 0, 0, 0, null;\n');
      fprintf(fd, '\n');
    end
  end

  fclose(fd);

  pkg load mboct-mbdyn-pkg;

  ad = [false, true];

  max_diff = 0;

  for k=1:numel(ad)
    fd = -1;
    unwind_protect
      fd = fopen("multibarmech_ad.con", "wt");
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
    mbdyn_solver_run("multibarmech_stage1.mbd", opts);
    opts.output_file = "stage2";
    opts.logfile = "stage2.stdout";
    mbdyn_solver_run("multibarmech_stage2.mbd", opts);

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
    fprintf(stderr, "[%d]: max difference: %e\n", k, curr_diff);
  endfor
  tol = 1e-10;
  assert(max_diff < tol);
catch
  gtest_error = lasterror();
  gtest_fail(gtest_error, __FILE__);
  rethrow(gtest_error);
end_try_catch
