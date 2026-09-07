BEGIN {
	edgeprop_add("arm_1", 7, -10.);
	edgeprop_add("arm_2", 8, -2);
	l_1 = .2;
	l_2 = .2;
	node_add("arm_1_top", 1, 0., 0., l_1/2, "hide");
	node_add("arm_1_bottom", 1, 0., 0., -l_1/2, "hide");
	node_add("arm_2_end", 1, -l_2, 0., 0., "hide");
	edge_add("arm_1", "arm_1_top", "arm_1_bottom", "arm_1");
	edge_add("arm_2", "arm_2_end", 1, "arm_2");
}
