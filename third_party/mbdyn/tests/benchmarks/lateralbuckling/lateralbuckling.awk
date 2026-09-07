BEGIN {
	h = 100e-3;
	b = 10e-3;

	for (n = 1000; n <= 1010; n++) {
		node_add(n "_C1", n, 0., b/2, h/2, "hide");
		node_add(n "_C2", n, 0., -b/2, h/2, "hide");
		node_add(n "_C3", n, 0., -b/2, -h/2, "hide");
		node_add(n "_C4", n, 0., b/2, -h/2, "hide");
	}

	sideprop_add("beam", 12);
	for (n = 1000; n < 1010; n++) {
		label1 = n;
		label2 = n + 1

		nodelabels[1] = label1 "_C1";
		nodelabels[2] = label1 "_C2";
		nodelabels[3] = label2 "_C2";
		nodelabels[4] = label2 "_C1";
		side_add(n "_C1_C2", 4, nodelabels, "beam");

		nodelabels[1] = label1 "_C2";
		nodelabels[2] = label1 "_C3";
		nodelabels[3] = label2 "_C3";
		nodelabels[4] = label2 "_C2";
		side_add(n "_C2_C3", 4, nodelabels, "beam");

		nodelabels[1] = label1 "_C3";
		nodelabels[2] = label1 "_C4";
		nodelabels[3] = label2 "_C4";
		nodelabels[4] = label2 "_C3";
		side_add(n "_C3_C4", 4, nodelabels, "beam");

		nodelabels[1] = label1 "_C4";
		nodelabels[2] = label1 "_C1";
		nodelabels[3] = label2 "_C1";
		nodelabels[4] = label2 "_C4";
		side_add(n "_C4_C1", 4, nodelabels, "beam");
	}

	nodelabels[1] = 1010 "_C1";
	nodelabels[2] = 1010 "_C2";
	nodelabels[3] = 1010 "_C3";
	nodelabels[4] = 1010 "_C4";
	side_add(1010, 4, nodelabels, "beam");

	edgeprop_add("rh", 7, 2.);

	node_add("C_p", 2000, 0., b/2, 0., "hide");
	node_add("C_m", 2000, 0., -b/2, 0., "hide");
	edge_add("C", "C_p", "C_m", "rh");

	node_add("B_p", 3000, 0., b/2, 0., "hide");
	node_add("B_m", 3000, 0., -b/2, 0., "hide");
	edge_add("B", "B_p", "B_m", "rh");

	node_add("G_p", 3002, 0., b/2, 0., "hide");
	node_add("G_m", 3002, 0., -b/2, 0., "hide");
	edge_add("G", "G_p", "G_m", "rh");
}
