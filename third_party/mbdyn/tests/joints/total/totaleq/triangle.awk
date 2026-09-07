BEGIN {
	nodeprop[nodeprop_num, "name"] = "actuator";
	nodeprop[nodeprop_num, "color"] = 8;
	nodeprop[nodeprop_num, "radius"] = 4.;
	nodeprop_num++;

	nodeprop[nodeprop_num, "name"] = "hinge";
	nodeprop[nodeprop_num, "color"] = 8;
	nodeprop[nodeprop_num, "radius"] = 3.;
	nodeprop_num++;

	for (i = 1; i < nodeprop_num; i++) {
		if (nodeprop[i, "name"] == "default") {
			nodeprop[i, "radius"] = 0;
		}
	}

	for (i = 1; i < edgeprop_num; i++) {
		if (edgeprop[i, "name"] == "rod_edge") {
			edgeprop[i, "color"] = 12;
		} else if (edgeprop[i, "name"] == "rod_offset") {
			edgeprop[i, "radius"] = 0;
		}
	}

	edgeprop[edgeprop_num, "name"] = "actuator";
	edgeprop[edgeprop_num, "color"] = 5;
	edgeprop[edgeprop_num, "radius"] = 2.;
	edgeprop_num++;

	sideprop[sideprop_num, "name"] = "chassis";
	sideprop[sideprop_num, "color"] = 7;
	sideprop_num++;

	sideprop[sideprop_num, "name"] = "end_effector";
	sideprop[sideprop_num, "color"] = 14;
	sideprop_num++;

	AD = 1.;
	HD = .5;

	# chassis
	for (i = 0; i < 3; i++) {
		node[node_num] = "chassis_" i + 1;
		node[node_num, "relative"] = 1000;
		node[node_num, 1] = AD*cos(i*2./3.*2.*3.14159);
		node[node_num, 2] = AD*sin(i*2./3.*2.*3.14159);
		node[node_num, 3] = 0.;
		node[node_num, "prop"] = "actuator";
		node_num++;
	}

	side[side_num] = "chassis";
	side[side_num, "N"] = 3;
	side[side_num, 1] = "chassis_1";
	side[side_num, 2] = "chassis_2";
	side[side_num, 3] = "chassis_3";
	side[side_num, "prop"] = "chassis";
	side_num++;

	# end effector
	for (i = 0; i < 3; i++) {
		node[node_num] = "end_effector_" i + 1;
		node[node_num, "relative"] = 2000;
		node[node_num, 1] = HD*cos(i*2./3.*2.*3.14159);
		node[node_num, 2] = HD*sin(i*2./3.*2.*3.14159);
		node[node_num, 3] = .01*AD;
		node[node_num, "prop"] = "hide";
		node_num++;
	}

	side[side_num] = "end_effector";
	side[side_num, "N"] = 3;
	side[side_num, 1] = "end_effector_1";
	side[side_num, 2] = "end_effector_2";
	side[side_num, 3] = "end_effector_3";
	side[side_num, "prop"] = "end_effector";
	side_num++;

	# actuators
	for (i = 0; i < 3; i++) {
		node[node_num] = "act_" i + 1 "_in";
		node[node_num, "relative"] = 4100 + 100*i;
		node[node_num, 1] = -HD;
		node[node_num, 2] = 0.;
		node[node_num, 3] = 0.;
		node[node_num, "prop"] = "hinge";
		node_num++;

		node[node_num] = "act_" i + 1 "_out";
		node[node_num, "relative"] = 4100 + 100*i;
		node[node_num, 1] = HD;
		node[node_num, 2] = 0.;
		node[node_num, 3] = 0.;
		node[node_num, "prop"] = "hide";
		node_num++;

		edge[edge_num] = "act_" i + 1;
		edge[edge_num, 1] = "act_" i + 1 "_in";
		edge[edge_num, 2] = "act_" i + 1 "_out";
		edge[edge_num, "prop"] = "actuator";
		edge_num++;
	}
}
