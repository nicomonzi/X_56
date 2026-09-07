
#include "mbconfig.h"
#include "matvec3.h"

#include <cmath>

int main(void) {
	int NumberOfTheta = 16;
	int NumberOfPhi   = 16;
	const doublereal HigherPhi = 0.8 * std::asin(1.);
	const doublereal Radius = 10.;
    	const doublereal For = 400.;
	

	std::cout << "begin: data;\n";
	std::cout << "	problem: initial value;\n";
	std::cout << "end: data;\n";

	std::cout << "begin: initial value;\n";
	std::cout << "	initial time: 0.;\n";
	std::cout << "	final time: 1.;\n";
	std::cout << "	time step: 1.e-2;\n";

	std::cout << "	tolerance: 1.e-6;\n";
	std::cout << "	max iterations: 10;\n";
	std::cout << "	derivatives tolerance: 1.E100;\n";

	std::cout << "	output: iterations;\n";
	std::cout << "	#output: residual;\n";
	std::cout << "	#output: jacobian;\n";
	std::cout << "	threads: disable;\n";
	// std::cout << "	output: residual, jacobian;\n";
	std::cout << "end: initial value;\n";

	std::cout << "begin: control data;\n";
	std::cout << "	structural nodes:\n";
	std::cout << (NumberOfTheta + 1) * (NumberOfPhi + 1) + 1 << "\n";
	std::cout << "	;\n";
	std::cout << "	plates:\n";
	std::cout << NumberOfTheta * NumberOfPhi << "\n";
	std::cout << "	;\n";
	std::cout << "	joints:\n";
	std::cout << 2 * 2 * (NumberOfPhi + 1) + 1 + 1<< "\n";
	std::cout << "	;\n";
	std::cout << "	forces:\n";
	std::cout << "		+2\n";
	std::cout << "	;\n";
	std::cout << "end: control data;\n";

	std::cout << "set: real E = 6.825E7;\n";
	std::cout << "set: real nu = 0.3;\n";
	std::cout << "set: real h = 0.04;\n";
	std::cout << "set: real as = 1.;\n";
	std::cout << "set: real at = 0.01;\n";
	std::cout << "set: real C = E*h/(1. - nu^2);\n";
	std::cout << "set: real D = E*h^3/12./(1. - nu^2);\n";

	std::cout << "begin: nodes;\n";

	// Nodes and constraints:
	// sphere with polar holes, axis along Z, centered
	// the eighth with X, Y, Z positive is modeled
	const doublereal dtheta = ( std::asin(1.) ) / NumberOfTheta;
	const doublereal dphi   = HigherPhi / NumberOfPhi;
	const bool Dofs111[3] = {1,1,1};
	std::cout << "structural: " << 1 << ", static,\n\t"
		<< 0. << ", " << 0. << ", " << 0. << ", \n\t"
		<< "eye, \n\t"
		<< "null, \n\t"
		<< "null; \n\n";
	for (int n1 = 0; n1 <= NumberOfTheta; n1++) {
	   	doublereal theta = dtheta * n1;
		for (int n2 = 0; n2 <= NumberOfPhi; n2++) {
			doublereal phi = dphi * n2;
			int LabelNode = (n1+1)*10000+(n2+1);
			Vec3 x(Radius*std::cos(theta)*std::cos(phi),Radius*std::sin(theta)*std::cos(phi),Radius*std::sin(phi));
			Vec3 alpha1, alpha2, alpha3;
			Mat3x3 alpha;
			alpha1 = Vec3(-std::sin(theta),std::cos(theta),0.);
			alpha3 = x/x.Norm();
			alpha2 = alpha3.Cross(alpha1);
			alpha = Mat3x3(alpha1, alpha2, alpha3);
			
			std::cout << "structural: " << LabelNode << ", static,\n\t"
				<< x(1) << ", " << x(2) << ", " << x(3) << ", \n\t"
				<< "eye,\n\t"
// 				<< "matr, \n\t\t"
// 				<< alpha(1,1) << ", " << alpha(1,2) << ", " << alpha(1,3) << ", \n\t\t"
// 				<< alpha(2,1) << ", " << alpha(2,2) << ", " << alpha(2,3) << ", \n\t\t"
// 				<< alpha(3,1) << ", " << alpha(3,2) << ", " << alpha(3,3) << ", \n\t"
				<< "null, \n\t"
				<< "null; \n\n";
		}
	}
				
std::cout << "end: nodes;\n";
std::cout << "begin: elements;\n";

	std::cout << "joint: 1, clamp, 1, node, node;\n\n";
	for (int n1 = 0; n1 <= NumberOfTheta; n1++) {
		for (int n2 = 0; n2 <= NumberOfPhi; n2++) {
			int LabelNode = (n1+1)*10000+(n2+1);
			int LabelConstr = (n1+1)*10000+(n2+1);
			int LabelConstr1 = 1000000000 + (n1+1)*10000+(n2+1);
			int LabelConstr2 = 2000000000 + (n1+1)*10000+(n2+1);
			if ( n1 == 0 ) {
				if ( n2 == 0 ) {
					//vincola uy uz rx rz
					std::cout << "joint: "<< LabelConstr << ", in plane,\n\t"
						<< "1,\n\t"
						<< "null,\n\t"
						<< "0., 1., 0.,\n\t"
						<< LabelNode << ";\n\n";
					std::cout << "joint: "<< LabelConstr2 << ", in plane,\n\t"
						<< "1,\n\t"
						<< "null,\n\t"
						<< "0., 0., 1.,\n\t"
						<< LabelNode << ";\n\n";
					std::cout << "joint: "<< LabelConstr1 << ", revolute rotation,\n\t"
						<< "1,\n\t"
						<< "hinge, reference, node, " << "3, 0., 1., 0., 1, guess,\n\t"
						<< LabelNode << ", "
						<< "hinge, reference, other node, " << "3, 0., 1., 0., 1, guess;\n\n";
// 					std::cout << "joint: " << LabelConstr << ", total joint,\n\t" 
// 						<< "1,\n\t"
// 						<< LabelNode << " ,\n\t"
// 						<< "position, reference, other position, null,\n\t"
// 						<< "rotation orientation, reference, other orientation, eye,\n\t"
// 						<< "position constraint, inactive, active, active,\n\t"
// 						<< "	null,\n\t"
// 						<< "orientation constraint, active, inactive, active,\n\t"
// 						<< "	null;\n\n";
				} else {
					//vincola uy rx rz
					std::cout << "joint: "<< LabelConstr << ", in plane,\n\t"
						<< "1,\n\t"
						<< "null,\n\t"
						<< "0., 1., 0.,\n\t"
						<< LabelNode << ";\n\n";
					std::cout << "joint: "<< LabelConstr1 << ", revolute rotation,\n\t"
						<< "1,\n\t"
						<< "hinge, reference, node, " << "3, 0., 1., 0., 1, guess,\n\t"
						<< LabelNode << ", "
						<< "hinge, reference, other node, " << "3, 0., 1., 0., 1, guess;\n\n";

// 					std::cout << "joint: " << LabelConstr << ", total joint,\n\t" 
// 						<< "1,\n\t"
// 						<< LabelNode << " ,\n\t"
// 						<< "position, reference, other position, null,\n\t"
// 						<< "rotation orientation, reference, other orientation, eye,\n\t"
// 						<< "position constraint, inactive, active, inactive,\n\t"
// 						<< "	null,\n\t"
// 						<< "orientation constraint, active, inactive, active,\n\t"
// 						<< "	null;\n\n";
				}
			} else if ( n1 == NumberOfTheta ) {
					//vincola ux ry rz
					std::cout << "joint: "<< LabelConstr << ", in plane,\n\t"
						<< "1,\n\t"
						<< "null,\n\t"
						<< "1., 0., 0.,\n\t"
						<< LabelNode << ";\n\n";
					std::cout << "joint: "<< LabelConstr1 << ", revolute rotation,\n\t"
						<< "1,\n\t"
						<< "hinge, reference, node, " << "3, 1., 0., 0., 1, guess,\n\t"
						<< LabelNode << ", "
						<< "hinge, reference, other node, " << "3, 1., 0., 0., 1, guess;\n\n";

// 					std::cout << "joint: " << LabelConstr << ", total joint,\n\t" 
// 						<< "1,\n\t"
// 						<< LabelNode << " ,\n\t"
// 						<< "position, reference, other position, null,\n\t"
// 						<< "rotation orientation, reference, other orientation, eye,\n\t"
// 						<< "position constraint, active, inactive, inactive,\n\t"
// 						<< "	null,\n\t"
// 						<< "orientation constraint, inactive, active, active,\n\t"
// 						<< "	null;\n\n";
			}
		}
	}

	// Shells:
	int ShellNodeConn[4];
	for (int n1 = 0; n1 < NumberOfTheta; n1++) {
		for (int n2 = 0; n2 < NumberOfPhi; n2++) {
			int LabelShell = 1000000 + (n1+1) * 10000 + (n2+1);
			ShellNodeConn[0] = (n1 + 1    ) * 10000 + (n2 + 1    );
			ShellNodeConn[1] = (n1 + 1 + 1) * 10000 + (n2 + 1    );
			ShellNodeConn[2] = (n1 + 1 + 1) * 10000 + (n2 + 1 + 1);
			ShellNodeConn[3] = (n1 + 1    ) * 10000 + (n2 + 1 + 1);
			std::cout << "shell4: " << LabelShell << ",\n\t"
				<< ShellNodeConn[0] << ", " 
				<< ShellNodeConn[1] << ", " 
				<< ShellNodeConn[2] << ", " 
				<< ShellNodeConn[3] << ", \n";
				std::cout << "	       C,		0,		 0,		  0,		C*nu,		    0,  	     0, 	      0,	       0,		0,		 0,		  0,\n";
				std::cout << "	       0,	 C*(1.-nu),		 0,		  0,		   0,		    0,  	     0, 	      0,	       0,		0,		 0,		  0,\n";
				std::cout << "	       0,		0, 1./2.*as*C*(1.-nu),		  0,		   0,		    0,  	     0, 	      0,	       0,		0,		 0,		  0,\n";
				std::cout << "	       0,		0,		 0,	   C*(1.-nu),		   0,		    0,  	     0, 	      0,	       0,		0,		 0,		  0,\n";
				std::cout << "	    C*nu,		0,		 0,		  0,		   C,		    0,  	     0, 	      0,	       0,		0,		 0,		  0,\n";
				std::cout << "	       0,		0,		 0,		  0,		   0, 1./2.*as*C*(1.-nu),  	     0, 	      0,	       0,		0,		 0,		  0,\n";
				std::cout << "	       0,		0,		 0,		  0,		   0,		    0,  	     D, 	      0,	       0,		0,	      D*nu,		  0,\n";
				std::cout << "	       0,		0,		 0,		  0,		   0,		    0,  	     0,        D*(1.-nu),	       0,		0,		 0,		  0,\n";
				std::cout << "	       0,		0,		 0,		  0,		   0,		    0,  	     0, 	      0,     at*D*(1.-nu),		0,		 0,		  0,\n";
				std::cout << "	       0,		0,		 0,		  0,		   0,		    0,  	     0, 	      0,	       0,	 D*(1.-nu),		 0,		  0,\n";
				std::cout << "	       0,		0,		 0,		  0,		   0,		    0,  	  D*nu, 	      0,	       0,		0,		 D,		  0,\n";
				std::cout << "	       0,		0,		 0,		  0,		   0,		    0,  	     0, 	      0,	       0,		0,		 0,	at*D*(1.-nu);\n";
		}
	}

	// Loads:
        int LoadNodeConn[2];
        LoadNodeConn[0] = (1) * 10000 + (1);
        LoadNodeConn[1] = (NumberOfTheta + 1) * 10000 + (1);

	std::cout << "force: 1, absolute,\n\t"
		<< LoadNodeConn[0] << ",\n\t"
		<< "position, null,\n\t"
		<< "1., 0., 0.,\n\t"
		<< "ramp, " << For / 2. << ",0., forever, 0.;\n\n";
	std::cout << "force: 2, absolute,\n\t"
		<< LoadNodeConn[1] << ",\n\t"
		<< "position, null,\n\t"
		<< "0., -1., 0.,\n\t"
		<< "ramp, " << For / 2. << ",0., forever, 0.;\n\n";


	std::cout << "end: elements;\n";
	return 0;
}
