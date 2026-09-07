# $Header: /var/cvs/mbdyn/mbdyn/mbdyn-1.0/tests/joints/drvhinge/invkin.awk,v 1.2 2007/02/26 00:34:46 masarati Exp $

BEGIN {
	nodeprop[nodeprop_num, "name"] = "member_joint";
	nodeprop[nodeprop_num, "color"] = 5;
	nodeprop[nodeprop_num, "radius"] = 2.;
	nodeprop_num++;

	edgeprop[edgeprop_num, "name"] = "member";
	edgeprop[edgeprop_num, "color"] = 13;
	edgeprop[edgeprop_num, "radius"] = 1.;
	edgeprop_num++;
	
	node[node_num] = "DX_FOOT_TOES";
	node[node_num, "relative"] = 2000;
	node[node_num, 1] = .1;
	node[node_num, 2] = .0;
	node[node_num, 3] = .0;
	node[node_num, "prop"] = "member_joint";
	node_num++;
	
	node[node_num] = "DX_ANKLE";
	node[node_num, "relative"] = 2000;
	node[node_num, 1] = -.05;
	node[node_num, 2] = .0;
	node[node_num, 3] = .0;
	node[node_num, "prop"] = "member_joint";
	node_num++;
	
	node[node_num] = "DX_KNEE";
	node[node_num, "relative"] = 4000;
	node[node_num, 1] = -.25;
	node[node_num, 2] = .0;
	node[node_num, 3] = .0;
	node[node_num, "prop"] = "member_joint";
	node_num++;
	
	edge[edge_num] = "DX_FOOT_ANKLE";
	edge[edge_num, 1] = "DX_FOOT_TOES";
	edge[edge_num, 2] = "DX_ANKLE";
	edge[edge_num, "prop"] = "member";
	edge_num++;

	edge[edge_num] = "DX_LOWER_LEG_ANKLE";
	edge[edge_num, 1] = 4000;
	edge[edge_num, 2] = "DX_ANKLE";
	edge[edge_num, "prop"] = "member";
	edge_num++;
	
	edge[edge_num] = "DX_LOWER_LEG_KNEE";
	edge[edge_num, 1] = 4000;
	edge[edge_num, 2] = "DX_KNEE";
	edge[edge_num, "prop"] = "member";
	edge_num++;

	edge[edge_num] = "DX_UPPER_LEG_KNEE";
	edge[edge_num, 1] = 6000;
	edge[edge_num, 2] = "DX_KNEE";
	edge[edge_num, "prop"] = "member";
	edge_num++;
}
