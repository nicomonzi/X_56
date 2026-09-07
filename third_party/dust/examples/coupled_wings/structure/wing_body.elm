# central body
body: GROUND + RIGHT, GROUND + RIGHT,
    (dL/2.)*m,
    reference, CENTRAL,  0, -Off_CG, 0.0,
    diag, (dL/2.)*j, 1./12.*(dL/2.)^3*m, 1./12.*(dL/2.)^3*m;
# right
body: GROUND + RIGHT + 1, GROUND + RIGHT + 1,
    (dL/2.)*m,
    reference, node,  0, -Off_CG, 0.0,
    diag, (dL/2.)*j, 1./12.*(dL/2.)^3*m, 1./12.*(dL/2.)^3*m;

body: GROUND + RIGHT + 2, GROUND + RIGHT + 2,
    (dL/2.)*m,
    reference, node,  0, -Off_CG, 0.0,
    diag, (dL/2.)*j, 1./12.*(dL/2.)^3*m, 1./12.*(dL/2.)^3*m;

body: GROUND + RIGHT + 3, GROUND + RIGHT + 3,
    (dL/2.)*m,
    reference, node,  0, -Off_CG, 0.0,
    diag, (dL/2.)*j, 1./12.*(dL/2.)^3*m, 1./12.*(dL/2.)^3*m;

body: GROUND + RIGHT + 4, GROUND + RIGHT + 4,
    (dL/4.)*m,
    reference, node,  -dL/8, -Off_CG, 0.0,
    diag, (dL/4.)*j, 1./12.*(dL/4.)^3*m, 1./12.*(dL/4.)^3*m;
# left
body: GROUND + LEFT + 1, GROUND + LEFT + 1,
    (dL/2.)*m,
    reference, node,  0, -Off_CG, 0.0,
    diag, (dL/2.)*j, 1./12.*(dL/2.)^3*m, 1./12.*(dL/2.)^3*m;

body: GROUND + LEFT + 2, GROUND + LEFT + 2,
    (dL/2.)*m,
    reference, node,  0, -Off_CG, 0.0,
    diag, (dL/2.)*j, 1./12.*(dL/2.)^3*m, 1./12.*(dL/2.)^3*m;

body: GROUND + LEFT + 3, GROUND + LEFT + 3,
    (dL/2.)*m,
    reference, node,  0, -Off_CG, 0.0,
    diag, (dL/2.)*j, 1./12.*(dL/2.)^3*m, 1./12.*(dL/2.)^3*m;

body: GROUND + LEFT + 4, GROUND + LEFT + 4,
    (dL/4.)*m,
    reference, node,  +dL/8, -Off_CG, 0.0,
    diag, (dL/4.)*j, 1./12.*(dL/4.)^3*m, 1./12.*(dL/4.)^3*m;