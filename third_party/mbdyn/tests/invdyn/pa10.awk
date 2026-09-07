BEGIN {
	edgeprop_add("arm_2", 7, -.5);
	sideprop_add("motor_caps", 0);

	np = 12;

	edgeprop_add("base", 7, -1.);
	edge_add("base", 0, 100, "base");

	l = .1;
	d = .06;
	link_label = 100;
	node_add("motor_" link_label "_l", link_label, 0., 0., l/2, "hide");
	node_add("motor_" link_label "_r", link_label, 0., 0., -l/2, "hide");
	edgeprop_add("motor_" link_label, 7, 100*d);
	edge_add("motor_" link_label, "motor_" link_label "_l", "motor_" link_label "_r", "motor_" link_label);
	for (i = 0; i < np; i++) {
		label_l = "motor_" link_label "_l_" i;
		node_add(label_l, link_label, d/2*cos(6.28*i/np), d/2*sin(6.28*i/np), l/2, "hide");
		nodes_l[i+1] = label_l;
		label_r = "motor_" link_label "_r_" i;
		node_add(label_r, link_label, d/2*cos(6.28*i/np), d/2*sin(6.28*i/np), -l/2, "hide");
		nodes_r[i+1] = label_r;
	}
	side_add("motor_" link_label "_l", np, nodes_l, "motor_caps");
	side_add("motor_" link_label "_r", np, nodes_r, "motor_caps");

	l = .3;
	d = .08;
	link_label = 200;
	node_add("motor_" link_label, link_label, 0., 0., l, "hide");
	edgeprop_add("motor_" link_label, 7, 100*d);
	edge_add("motor_" link_label, link_label, "motor_" link_label, "motor_" link_label);
	for (i = 0; i < np; i++) {
		label = "motor_" link_label "_" i;
		node_add(label, link_label, d/2*cos(6.28*i/np), d/2*sin(6.28*i/np), l, "hide");
		nodes[i+1] = label;
	}
	side_add("motor_" link_label, np, nodes, "motor_caps");

	link_label = 300;
	edgeprop_add("motor_" link_label, 7, 60*d);
	edge_add("motor_" link_label, link_label, "motor_" link_label - 100, "motor_" link_label);

	l = .1;
	d = .06;
	link_label = 300;
	node_add("motor_" link_label "_l", link_label, 0., 0., l/2, "hide");
	node_add("motor_" link_label "_r", link_label, 0., 0., -l/2, "hide");
	edgeprop_add("motor_" link_label, 7, 100*d);
	edge_add("motor_" link_label, "motor_" link_label "_l", "motor_" link_label "_r", "motor_" link_label);
	for (i = 0; i < np; i++) {
		label_l = "motor_" link_label "_l_" i;
		node_add(label_l, link_label, d/2*cos(6.28*i/np), d/2*sin(6.28*i/np), l/2, "hide");
		nodes_l[i+1] = label_l;
		label_r = "motor_" link_label "_r_" i;
		node_add(label_r, link_label, d/2*cos(6.28*i/np), d/2*sin(6.28*i/np), -l/2, "hide");
		nodes_r[i+1] = label_r;
	}
	side_add("motor_" link_label "_l", np, nodes_l, "motor_caps");
	side_add("motor_" link_label "_r", np, nodes_r, "motor_caps");

	l = .2;
	d = .08;
	link_label = 400;
	node_add("motor_" link_label, link_label, 0., 0., l, "hide");
	edgeprop_add("motor_" link_label, 7, 100*d);
	edge_add("motor_" link_label, link_label, "motor_" link_label, "motor_" link_label);
	for (i = 0; i < np; i++) {
		label = "motor_" link_label "_" i;
		node_add(label, link_label, d/2*cos(6.28*i/np), d/2*sin(6.28*i/np), l, "hide");
		nodes[i+1] = label;
	}
	side_add("motor_" link_label, np, nodes, "motor_caps");

	link_label = 500;
	edgeprop_add("motor_" link_label, 7, 60*d);
	edge_add("motor_" link_label, link_label, "motor_" link_label - 100, "motor_" link_label);

	l = .12;
	d = .05;
	link_label = 500;
	node_add("motor_" link_label "_l", link_label, 0., 0., l/2, "hide");
	node_add("motor_" link_label "_r", link_label, 0., 0., -l/2, "hide");
	edgeprop_add("motor_" link_label, 7, 100*d);
	edge_add("motor_" link_label, "motor_" link_label "_l", "motor_" link_label "_r", "motor_" link_label);
	for (i = 0; i < np; i++) {
		label_l = "motor_" link_label "_l_" i;
		node_add(label_l, link_label, d/2*cos(6.28*i/np), d/2*sin(6.28*i/np), l/2, "hide");
		nodes_l[i+1] = label_l;
		label_r = "motor_" link_label "_r_" i;
		node_add(label_r, link_label, d/2*cos(6.28*i/np), d/2*sin(6.28*i/np), -l/2, "hide");
		nodes_r[i+1] = label_r;
	}
	side_add("motor_" link_label "_l", np, nodes_l, "motor_caps");
	side_add("motor_" link_label "_r", np, nodes_r, "motor_caps");

	l = .04
	d = .04;
	link_label = 600;
	node_add("motor_" link_label, link_label, 0., 0., l, "hide");
	edgeprop_add("motor_" link_label, 7, 100*d);
	edge_add("motor_" link_label, link_label, "motor_" link_label, "motor_" link_label);
	for (i = 0; i < np; i++) {
		label = "motor_" link_label "_" i;
		node_add(label, link_label, d/2*cos(6.28*i/np), d/2*sin(6.28*i/np), l, "hide");
		nodes[i+1] = label;
	}
	side_add("motor_" link_label, np, nodes, "motor_caps");

	link_label = 700;
	edgeprop_add("motor_" link_label, 7, 60*d);
	edge_add("motor_" link_label, link_label, "motor_" link_label - 100, "motor_" link_label);

	for (i = 0; i < np; i++) {
		label = "motor_" link_label "_" i;
		node_add(label, link_label, .6*d/2*cos(6.28*i/np), .6*d/2*sin(6.28*i/np), 0., "hide");
		nodes[i+1] = label;
	}
	side_add("motor_" link_label, np, nodes, "motor_caps");
}
