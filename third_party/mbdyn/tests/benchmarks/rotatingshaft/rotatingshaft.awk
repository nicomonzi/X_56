BEGIN {
	if (DISK_NODE == 0) {
		DISK_NODE = 16;
	}

	RD = 0.24;
	TD = 0.05;
	L = 6;

	node_absolute_add("end_1", L, -RD, -RD, "hide");
	node_absolute_add("end_2", L, RD, -RD, "hide");
	node_absolute_add("end_3", L, RD, RD, "hide");
	node_absolute_add("end_4", L, -RD, RD, "hide");

	nodes[1] = "end_1";
	nodes[2] = "end_2";
	nodes[3] = "end_3";
	nodes[4] = "end_4";
	side_add("end", 4, nodes, "disk");

	NP = 16;
	for (i = 0; i < NP; i++) {
		THETA = 6.28*i/NP;
		node_add("disk_top_" i, DISK_NODE, TD/2, RD*cos(THETA), RD*sin(THETA), "hide");
		node_add("disk_bottom_" i, DISK_NODE, -TD/2, RD*cos(THETA), RD*sin(THETA), "hide");
	}
	node_add("disk_top", DISK_NODE, TD/2, 0., 0., "hide");
	node_add("disk_bottom", DISK_NODE, -TD/2, 0., 0., "hide");
	edgeprop_add("disk", 12, -RD);
	edge_add("disk_side", "disk_bottom", "disk_top", "disk");

	edgeprop_add("disk_radius", 0, 1);
	edge_add("disk_radius", "disk_top", "disk_top_" 0, "disk_radius");

	sideprop_add("disk", 12);

	for (i = 0; i < NP; i++) {
		nodes[i + 1] = "disk_top_" i;
	}
	side_add("disk_top", NP, nodes, "disk");
	for (i = 0; i < NP; i++) {
		nodes[i + 1] = "disk_bottom_" i;
	}
	side_add("disk_bottom", NP, nodes, "disk");
}
