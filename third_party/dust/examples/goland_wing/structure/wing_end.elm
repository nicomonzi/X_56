
body: 2*curr_elem-1, 2*curr_elem-1,
    (dL/2.)*m,
    reference, node, Off_CG, 0.0, 0.0,
    diag, 1./12.*(dL/2.)^3*m, (dL/2.)*j, 1./12.*(dL/2.)^3*m;

body: 2*curr_elem, 2*curr_elem,
    (dL/4.)*m,
    reference, node, Off_CG, -dL/8, 0.0,
    diag,  1./12.*(dL/4.)^3*m, (dL/4.)*j, 1./12.*(dL/4.)^3*m;

beam3: curr_elem,
    2*curr_elem-2, null,
    2*curr_elem-1, null,
    2*curr_elem  , null,
    eye,
    linear time variant viscoelastic generic,
	    diag, 
            EA, GAy, GAz, 
            EJx, GJ, EJz,
            const, 1.,
        diag, 0,0,0,
            0, 0*3e-3*GJ, 0,
	    const, 1,
            #ramp, -2.5 , 0.0, 0.4, 1.,
    same,
    same;


#aerodynamic beam3: curr_elem, curr_elem,
#    reference, node, Off_CA, 0.0, 0.0, 
#    1, -1, 0, 0, 
#    2,  0, 0, 1,
#    reference, node, Off_CA, 0.0, 0.0,
#    1, -1, 0, 0, 
#    2,  0, 0, 1,
#    reference, node, Off_CA, 0.0, 0.0,
#    1, -1, 0, 0, 
#    2,  0, 0, 1,
#    const,  Chord,           #chord
#    const,  0.0,           #AC offset
#    const, -Chord/2,          #collocation point
#    const,  0.0,           #twist
#    5, theodorsen, c81, naca0012, 
#    jacobian, yes;      #gauss points and c81 table