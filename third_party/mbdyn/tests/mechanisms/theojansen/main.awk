BEGIN {
	# DEF
	x = -7.8880318341652327e+01 + 3.8792267270681663e+01;
	y = 0;
	z = -9.7885715569134935e-01 + 0.00000000000000000+00;
	node_add("front_def", 50, x, y, z);
	node_add("rear_def", 150, -x, y, z);

	# BEJ
	x = -4.2503343798562575e+01 + 3.8792267270681663e+01;
	y = 0;
	z = 4.1333738168767304e+01 + 0.00000000000000000+00;
	node_add("front_bej", 50, x, y, z);
	node_add("rear_bej", 150, -x, y, z);

	# HI
	x = -2.3003592298551194e+01 + 3.1074706088454597e+01;
	y = 0;
	z = -8.6865485527572289e+01 + 3.8534779737252705e+01;
	node_add("front_hi", 60, x, y, z);
	node_add("rear_hi", 160, -x, y, z);

	# FGH
	x = -6.7773879688160434e+01 + 3.1074706088454597e+01;
	y = 0;
	z = -3.8781066356162967e+01 + 3.8534779737252705e+01;
	node_add("front_fgh", 60, x, y, z);
	node_add("rear_fgh", 160, -x, y, z);

	nodes[1] = "front_def";
	nodes[2] = "front_bej";
	nodes[3] = 50;
	side_add("front_top", 3, nodes, 14);

	nodes[2] = "rear_def";
	nodes[1] = "rear_bej";
	nodes[3] = 150;
	side_add("rear_top", 3, nodes, 14);

	nodes[1] = "front_hi";
	nodes[2] = "front_fgh";
	nodes[3] = 60;
	side_add("front_bottom", 3, nodes, 14);

	nodes[2] = "rear_hi";
	nodes[1] = "rear_fgh";
	nodes[3] = 160;
	side_add("rear_bottom", 3, nodes, 14);
}
