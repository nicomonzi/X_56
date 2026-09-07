	body: LINK_O_F, LINK_O_F, O_F_M, reference, LINK_O_F, O_F_CM_X, O_F_CM_Y, 0., diag, 0., 0., O_F_IZ;
	body: LINK_E_F, LINK_E_F, E_F_M, reference, LINK_E_F, E_F_CM_X, E_F_CM_Y, 0., diag, 0., 0., E_F_IZ;
	body: LINK_H_E, LINK_H_E, H_E_M, reference, LINK_H_E, H_E_CM_X, H_E_CM_Y, 0., diag, 0., 0., H_E_IZ;
	body: LINK_G_E, LINK_G_E, G_E_M, reference, LINK_G_E, G_E_CM_X, G_E_CM_Y, 0., diag, 0., 0., G_E_IZ;
	body: LINK_A_G, LINK_A_G, A_G_M, reference, LINK_A_G, A_G_CM_X, A_G_CM_Y, 0., diag, 0., 0., A_G_IZ;
	body: LINK_A_H, LINK_A_H, A_H_M, reference, LINK_A_H, A_H_CM_X, A_H_CM_Y, 0., diag, 0., 0., A_H_IZ;
	body: LINK_E_B_D, LINK_E_B_D, E_B_D_M, reference, LINK_E_B_D, E_B_D_CM_X, E_B_D_CM_Y, 0., diag, 0., 0., E_B_D_IZ;

	joint: POINT_C, clamp, POINT_C, node, node;

	joint: POINT_O, total pin joint,
		LINK_O_F,
			position, reference, POINT_O, null,
			position orientation, reference, POINT_O, eye,
			rotation orientation, reference, POINT_O,
				# eye,
				3, 0., 0., 1.,
				1, POINT_F_X - POINT_O_X, POINT_F_Y - POINT_O_Y, 0.,
		# POINT_O
			position, reference, POINT_O, null,
			position orientation, reference, POINT_O, eye,
			rotation orientation, reference, POINT_O, eye,
		position constraint, 1, 1, 1, null,
		orientation constraint, 1, 1, 0, null;

	joint: POINT_F, total joint,
		LINK_O_F,
			position, reference, POINT_F, null,
			position orientation, reference, POINT_F, eye,
			rotation orientation, reference, POINT_F, eye,
		LINK_E_F,
			position, reference, POINT_F, null,
			position orientation, reference, POINT_F, eye,
			rotation orientation, reference, POINT_F, eye,
		position constraint, 1, 1, 0, null,
		orientation constraint, 0, 0, 0, null;

	joint: POINT_E + 1, total joint,
		LINK_E_B_D,
			position, reference, POINT_E, null,
			position orientation, reference, POINT_E, eye,
			rotation orientation, reference, POINT_E, eye,
		LINK_E_F,
			position, reference, POINT_E, null,
			position orientation, reference, POINT_E, eye,
			rotation orientation, reference, POINT_E, eye,
		position constraint, 1, 1, 1, null,
		orientation constraint, 1, 1, 0, null;

	joint: POINT_E + 2, total joint,
		LINK_E_B_D,
			position, reference, POINT_E, null,
			position orientation, reference, POINT_E, eye,
			rotation orientation, reference, POINT_E, eye,
		LINK_G_E,
			position, reference, POINT_E, null,
			position orientation, reference, POINT_E, eye,
			rotation orientation, reference, POINT_E, eye,
		position constraint, 1, 1, 1, null,
		orientation constraint, 1, 1, 0, null;

	joint: POINT_E + 3, total joint,
		LINK_E_B_D,
			position, reference, POINT_E, null,
			position orientation, reference, POINT_E, eye,
			rotation orientation, reference, POINT_E, eye,
		LINK_H_E,
			position, reference, POINT_E, null,
			position orientation, reference, POINT_E, eye,
			rotation orientation, reference, POINT_E, eye,
		position constraint, 1, 1, 1, null,
		orientation constraint, 1, 1, 0, null;

	joint: POINT_G, total joint,
		LINK_A_G,
			position, reference, POINT_G, null,
			position orientation, reference, POINT_G, eye,
			rotation orientation, reference, POINT_G, eye,
		LINK_G_E,
			position, reference, POINT_G, null,
			position orientation, reference, POINT_G, eye,
			rotation orientation, reference, POINT_G, eye,
		position constraint, 1, 1, 1, null,
		orientation constraint, 0, 0, 0, null;

	joint: POINT_H, total joint,
		LINK_A_H,
			position, reference, POINT_H, null,
			position orientation, reference, POINT_H, eye,
			rotation orientation, reference, POINT_H, eye,
		LINK_H_E,
			position, reference, POINT_H, null,
			position orientation, reference, POINT_H, eye,
			rotation orientation, reference, POINT_H, eye,
		position constraint, 1, 1, 1, null,
		orientation constraint, 0, 0, 0, null;

	joint: POINT_A + 1, total pin joint,
		LINK_A_G,
			position, reference, POINT_A, null,
			position orientation, reference, POINT_A, eye,
			rotation orientation, reference, POINT_A, eye,
		# POINT_A
			position, reference, POINT_A, null,
			position orientation, reference, POINT_A, eye,
			rotation orientation, reference, POINT_A, eye,
		position constraint, 1, 1, 0, null,
		orientation constraint, 1, 1, 0, null;

	joint: POINT_A + 2, total pin joint,
		LINK_A_H,
			position, reference, POINT_A, null,
			position orientation, reference, POINT_A, eye,
			rotation orientation, reference, POINT_A, eye,
		# POINT_A
			position, reference, POINT_A, null,
			position orientation, reference, POINT_A, eye,
			rotation orientation, reference, POINT_A, eye,
		position constraint, 1, 1, 0, null,
		orientation constraint, 1, 1, 0, null;

	joint: POINT_B, total pin joint,
		LINK_E_B_D,
			position, reference, POINT_B, null,
			position orientation, reference, POINT_B, eye,
			rotation orientation, reference, POINT_B, eye,
		# POINT_B
			position, reference, POINT_B, null,
			position orientation, reference, POINT_B, eye,
			rotation orientation, reference, POINT_B, eye,
		position constraint, 1, 1, 1, null,
		orientation constraint, 1, 1, 0, null;

	joint: POINT_D, rod,
		POINT_C,
			position, reference, POINT_C, null,
		LINK_E_B_D,
			position, reference, POINT_D, null,
		L0,
			# linear elastic, 0;
			linear elastic, K*L0; # remember, the "rod" wants EA = K * L0

	couple: LINK_O_F, absolute,
		LINK_O_F,
			position, reference, POINT_O, null,
		0., 0., 1.,
			# null;
			const, TAU;

/*
drive caller: 0,
	string, "model::element::body(LINK_O_F, \"E\") \
		+ model::element::body(LINK_E_F, \"E\") \
		+ model::element::body(LINK_H_E, \"E\") \
		+ model::element::body(LINK_G_E, \"E\") \
		+ model::element::body(LINK_A_G, \"E\") \
		+ model::element::body(LINK_A_H, \"E\") \
		+ model::element::body(LINK_E_B_D, \"E\")";
drive caller: 1,
	string, "K/2*(model::element::joint(POINT_D, \"L\") - L0)^2";

drive caller: 2,
	string, "TAU*(model::element::joint(POINT_O, \"rz\") - BETA0)";

	force: POINT_C, absolute,
		POINT_C,
			position, null,
		component,
			reference, 0,
			reference, 1,
			reference, 2;
*/