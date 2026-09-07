BEGIN {
	for (bay = 0; bay < 30000; bay += 10000) {
		# DEF
		if (bay % 20000) {
			x = -6.9110362548148387e+01 + 3.8792267270681663e+01;
			y = 0;
			z = 2.6245439579981319e+01 + 0.00000000000000000+00;
		} else {
			x = -7.8880318341652327e+01 + 3.8792267270681663e+01;
			y = 0;
			z = -9.7885715569134935e-01 + 0.00000000000000000+00;
		}
		node_add("front_def_" bay, bay + 50, x, y, z);
		node_add("rear_def_" bay, bay + 150, -x, y, z);

		# BEJ
		if (bay % 20000) {
			x = -1.3730239356693771e+01 + 3.8792267270681663e+01;
			y = 0;
			z = 3.3077859012313525e+01 + 0.00000000000000000+00;
		} else {
			x = -4.2503343798562575e+01 + 3.8792267270681663e+01;
			y = 0;
			z = 4.1333738168767304e+01 + 0.00000000000000000+00;
		}
		node_add("front_bej_" bay, bay + 50, x, y, z);
		node_add("rear_bej_" bay, bay + 150, -x, y, z);

		# HI
		if (bay % 20000) {
			x = -8.5550520149076192e+01 + 5.9007027445073035e+01;
			y = 0;
			z = -7.4890323948268247e+01 + 3.3702425299848080e+01;
		} else {
			x = -2.3003592298551194e+01 + 3.1074706088454597e+01;
			y = 0;
			z = -8.6865485527572289e+01 + 3.8534779737252705e+01;
		}
		node_add("front_hi_" bay, bay + 60, x, y, z);
		node_add("rear_hi_" bay, bay + 160, -x, y, z);

		# FGH
		if (bay % 20000) {
			x = -8.6324982186682846e+01 + 5.9007027445073035e+01;
			y = 0;
			z = -9.1948887297364656e+00 + 3.3702425299848080e+01;
		} else {
			x = -6.7773879688160434e+01 + 3.1074706088454597e+01;
			y = 0;
			z = -3.8781066356162967e+01 + 3.8534779737252705e+01;
		}
		node_add("front_fgh_" bay, bay + 60, x, y, z);
		node_add("rear_fgh_" bay, bay + 160, -x, y, z);

		nodes[1] = "front_def_" bay;
		nodes[2] = "front_bej_" bay;
		nodes[3] = bay + 50;
		side_add("front_top_" bay, 3, nodes, 14);

		nodes[2] = "rear_def_" bay;
		nodes[1] = "rear_bej_" bay;
		nodes[3] = bay + 150;
		side_add("rear_top_" bay, 3, nodes, 14);

		nodes[1] = "front_hi_" bay;
		nodes[2] = "front_fgh_" bay;
		nodes[3] = bay + 60;
		side_add("front_bottom_" bay, 3, nodes, 14);

		nodes[2] = "rear_hi_" bay;
		nodes[1] = "rear_fgh_" bay;
		nodes[3] = bay + 160;
		side_add("rear_bottom_" bay, 3, nodes, 14);
	}
}
