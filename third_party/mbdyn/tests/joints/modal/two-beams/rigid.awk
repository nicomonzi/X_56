BEGIN {
	edgeprop[edgeprop_num, "name"] = "hinge";
	edgeprop[edgeprop_num, "color"] = 5;
	edgeprop[edgeprop_num, "radius"] = 4.;
	edgeprop_num++;

	edgeprop[edgeprop_num, "name"] = "body";
	edgeprop[edgeprop_num, "color"] = 7;
	edgeprop[edgeprop_num, "radius"] = 2.;
	edgeprop_num++;

	L = .5;
	H = .05;

	# body 1
	node[node_num] = "body1_1";
	node[node_num, "relative"] = 2;
	node[node_num, 1] = -L/2.;
	node[node_num, 2] = 0.;
	node[node_num, 3] = 0.;
	node[node_num, "prop"] = "hide";
	node_num++;

	node[node_num] = "body1_2";
	node[node_num, "relative"] = 2;
	node[node_num, 1] = L/2.;
	node[node_num, 2] = 0.;
	node[node_num, 3] = 0.;
	node[node_num, "prop"] = "hide";
	node_num++;

	# body 2
	node[node_num] = "body2_1";
	node[node_num, "relative"] = 5;
	node[node_num, 1] = -L/2.;
	node[node_num, 2] = 0.;
	node[node_num, 3] = 0.;
	node[node_num, "prop"] = "hide";
	node_num++;

	node[node_num] = "body2_2";
	node[node_num, "relative"] = 5;
	node[node_num, 1] = L/2.;
	node[node_num, 2] = 0.;
	node[node_num, 3] = 0.;
	node[node_num, "prop"] = "hide";
	node_num++;

	# hinge 1
	node[node_num] = "hinge1_1";
	node[node_num, "relative"] = 2;
	node[node_num, 1] = -L/2.;
	node[node_num, 2] = 0.;
	node[node_num, 3] = -H/2.;
	node[node_num, "prop"] = "hide";
	node_num++;

	node[node_num] = "hinge1_2";
	node[node_num, "relative"] = 2;
	node[node_num, 1] = -L/2.;
	node[node_num, 2] = 0.;
	node[node_num, 3] = H/2.;
	node[node_num, "prop"] = "hide";
	node_num++;

	# hinge 2
	node[node_num] = "hinge2_1";
	node[node_num, "relative"] = 5;
	node[node_num, 1] = -L/2.;
	node[node_num, 2] = 0.;
	node[node_num, 3] = -H/2.;
	node[node_num, "prop"] = "hide";
	node_num++;

	node[node_num] = "hinge2_2";
	node[node_num, "relative"] = 5;
	node[node_num, 1] = -L/2.;
	node[node_num, 2] = 0.;
	node[node_num, 3] = H/2.;
	node[node_num, "prop"] = "hide";
	node_num++;

	# body 1
	edge[edge_num] = "body1";
	edge[edge_num, 1] = "body1_1";
	edge[edge_num, 2] = "body1_2";
	edge[edge_num, "prop"] = "body";
	edge_num++;

	# body 2
	edge[edge_num] = "body2";
	edge[edge_num, 1] = "body2_1";
	edge[edge_num, 2] = "body2_2";
	edge[edge_num, "prop"] = "body";
	edge_num++;

	# hinge 1
	edge[edge_num] = "hinge1";
	edge[edge_num, 1] = "hinge1_1";
	edge[edge_num, 2] = "hinge1_2";
	edge[edge_num, "prop"] = "hinge";
	edge_num++;

	# hinge 2
	edge[edge_num] = "hinge2";
	edge[edge_num, 1] = "hinge2_1";
	edge[edge_num, 2] = "hinge2_2";
	edge[edge_num, "prop"] = "hinge";
	edge_num++;
}
