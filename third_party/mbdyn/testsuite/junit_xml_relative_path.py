#!/usr/bin/env python3

# MBDyn (C) is a multibody analysis code.
# http://www.mbdyn.org
#
# Copyright (C) 1996-2025
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
# Copyright (C) 2025(-2025) all rights reserved.

# The copyright of this code is transferred
# to Pierangelo Masarati and Paolo Mantegazza
# for use in the software MBDyn as described
# in the GNU Public License version 2.1

# Convert all absolute paths from Google-Test to relative paths which can be handled by GitLab.

import xml.etree.ElementTree as ET
import sys

file_names = [];
prefix = ''

i = 1;

while i < len(sys.argv):
    match sys.argv[i]:
        case '-p':
            i = i + 1;
            prefix = sys.argv[i];
        case _:
            file_names.append(sys.argv[i]);
    i = i + 1;

if len(prefix) == 0:
    print('missing argument -p');
    exit(1);

for file_name in file_names:
    tree = ET.parse(file_name)
    root = tree.getroot()

    for testsuites in root.iter('testsuites'):
        for testsuite in testsuites.iter('testsuite'):
            for testcase in testsuite.iter('testcase'):
                fn = testcase.get('file');
                if (not fn is None):
                    fn = fn.replace(prefix, '');
                    testcase.set('file', fn);
                name = testcase.get('name');
                if (not name is None):
                    name = name.replace(prefix, '');
                    testcase.set('name', name);
                valpar = testcase.get('value_param');
                if (not valpar is None):
                    valpar = valpar.replace(prefix, '');
                    testcase.set('value_param', valpar);
                for failure in testcase.iter('failure'):
                    msg = failure.get('message');
                    if (not msg is None):
                        msg = msg.replace(prefix, '');
                        failure.set('message', msg);
                    tp = failure.get('type');
                    if (not tp is None):
                        tp = tp.replace(prefix, '');
                        failure.set('type', tp);
                    data = failure.text;
                    if (not data is None):
                        data = data.replace(prefix, '');
                        failure.text = data;

    tree.write(file_name);
    print('File %s converted' % file_name);
