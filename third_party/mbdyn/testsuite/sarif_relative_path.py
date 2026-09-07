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

## Convert absolute paths in *.sarif files into relative paths which can be handled by GitLab.

import json
import os
import sys

root = "";

file_names = [];

i = 1;

while i < len(sys.argv):
    match sys.argv[i]:
        case "--root":
            i = i + 1;
            root = sys.argv[i];
        case _:
            file_names.append(sys.argv[i]);
    i = i + 1;

for file_name in file_names:
    print("processing %s" % file_name);

    with open(file_name) as f:
        sarif = json.load(f);

    for run in sarif.get("runs", []):
        for art in run.get("artifacts", []):
            uri = art["location"]["uri"];
            if uri.startswith(root):
                art["location"]["uri"] = uri[len(root):];
        for res in run.get("results", []):
            for cf in res.get("codeFlows", []):
                for tf in cf.get("threadFlows", []):
                    for loc in tf["locations"]:
                        uri = loc["location"]["physicalLocation"]["artifactLocation"]["uri"];
                        if uri.startswith(root):
                            loc["location"]["physicalLocation"]["artifactLocation"]["uri"] = uri[len(root):];
            for loc in res.get("locations", []):
                uri = loc["physicalLocation"]["artifactLocation"]["uri"];
                if uri.startswith(root):
                    loc["physicalLocation"]["artifactLocation"]["uri"] = uri[len(root):];

    with open(file_name, "w") as f:
        json.dump(sarif, f, indent=2);
