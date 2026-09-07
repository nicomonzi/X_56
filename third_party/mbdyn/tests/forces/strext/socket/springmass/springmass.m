% $Header: /var/cvs/mbdyn/mbdyn/mbdyn-1.0/tests/forces/strext/socket/springmass/springmass.m,v 1.1 2011/12/06 20:23:10 masarati Exp $

clear all
clc
close all
clear pymex

script = {
'import matlab',
'import spread',
'import sys',
'sys.path.append(''/usr/local/mbdyn/libexec/mbpy/'')',
'from mbc_py_interface import mbcNodal',
'from mbc_py import mbc_py_nodal_initialize',
'import os',
'import tempfile',
'tmpdir = tempfile.mkdtemp('''', ''.mbdyn_'');',
'path = tmpdir + ''/mbdyn.sock'';',
'os.environ[''MBSOCK''] = path;',
'os.system(''./mbdyn.sh -f springmass -o output > output.txt 2>&1 &'');',
'from numpy import *',
'host = "";',
'port = 0;',
'timeout = -1;	# forever',
'verbose = 1;',
'data_and_next = 1;',
'refnode = 0;	# no reference node',
'nodes = 1;	# number of nodes',
'labels = 0;',
'rot = 0x100;	# orientation vector',
'# rot = 0x200;	# orientation matrix',
'# rot = 0x400;	# Euler 123',
'accels = 0;',
'structure = mbcNodal(path, host, port, timeout, verbose, data_and_next, refnode, nodes, labels, rot, accels);',
'structure.negotiate();',
'structure.n_f[0] = 0.;',
'structure.n_f[1] = 0.;',
'structure.n_f[2] = 0.;',
'structure.n_m[0] = 0.;',
'structure.n_m[1] = 0.;',
'structure.n_m[2] = 0.;'
};

pymex(script);

m = 1.;
omega = 10*2*pi;
xi = 0.1;
k = m*omega*omega;
d = 2*xi*m*omega;
maxiter = 5;
i = 0;

while 1
    pymex('matlab.push("aa",(structure.recv()));');
    if aa break; end
    
    pymex('matlab.push("x0",structure.n_x[0]);');
    pymex('matlab.push("x1",structure.n_x[1]);');
    pymex('matlab.push("x2",structure.n_x[2]);');
    pymex('matlab.push("xp0",structure.n_xp[0]);');
    pymex('matlab.push("xp1",structure.n_xp[1]);');
    pymex('matlab.push("xp2",structure.n_xp[2]);');
    
    f0 = -k*x0 - d*xp0;
    f1 = -k*x1 - d*xp1;
    f2 = -k*x2 - d*xp2;
    
    pymex('structure.n_f[0] = matlab.pull("f0");');
    pymex('structure.n_f[1] = matlab.pull("f1");');
    pymex('structure.n_f[2] = matlab.pull("f2");');

    i = i + 1;
    if (i==maxiter); i = 0; end
    pymex('ii = matlab.pull("i");');
    pymex('matlab.push("cc",(structure.send(ii == 0)));');
    if cc break; end

end
pymex('structure.destroy(); os.remove(tmpdir);');

mov = load(['output.mov']);

tt = (0:length(mov)-1)*0.01;

aa = length(mov);

figure
subplot(1,2,1)
hold on
plot(tt(1:aa),mov(1:aa,4),'-b','linewidth',2)
xlim([0 2])
xlabel('t (sec)','fontsize',14)
ylabel('q','fontsize',14)
grid on
box on

subplot(1,2,2)
hold on
plot(tt(1:aa),mov(1:aa,10),'-b','linewidth',2)
xlim([0 2])
ylim([-1 1])
xlabel('t (sec)','fontsize',14)
ylabel('dq / dt','fontsize',14)
grid on
box on
