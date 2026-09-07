	body: CRANK, CRANK,
		0.12,
		null,
		diag, 0.0001, 0.00001, 0.0001;

	body: ROD, ROD,
		0.5,
		null,
		diag, 0.004, 0.0004, 0.004;

	body: BLOCK, BLOCK,
		2.,
		null,
		diag, 0.0001, 0.0001, 0.0001;

	joint: POINT_A, revolute pin,
		CRANK,
			position, reference, POINT_A, null,
			orientation, reference, POINT_A, 3, 1., 0., 0., 2, 0., 1., 0.,
		position, reference, POINT_A, null,
		orientation, reference, POINT_A, 3, 1., 0., 0., 2, 0., 1., 0.;

	joint: POINT_B, spherical hinge,
		CRANK,
			position, reference, POINT_B, null,
			orientation, reference, POINT_B, eye,
		ROD,
			position, reference, POINT_B, null,
			orientation, reference, POINT_B, eye;

	joint: POINT_C, cardano hinge,
		BLOCK,
			position, reference, POINT_C, null,
			orientation, reference, node, 3, 1., 0., 0., 2, guess,
		ROD,
			position, reference, POINT_C, null,
			orientation, reference, node, 2, 1., 0., 0., 3, guess;

	joint: POINT_C + 1, total pin joint,
		BLOCK,
			position, reference, POINT_C, null,
			position orientation, reference, POINT_C, eye,
			rotation orientation, reference, POINT_C, eye,
		position, reference, POINT_C, null,
		position orientation, reference, POINT_C, eye,
		rotation orientation, reference, POINT_C, eye,
		position constraint, 0, 1, 1, null,
		orientation constraint, 1, 1, 1, null;

	gravity: 0., 0., -1., const, 9.81;