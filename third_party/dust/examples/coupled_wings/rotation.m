clear all, close all


%%
syms lambda delta theta real

R_rot = [0  1   0;
        -1  0   0;
         0  0   1];
     
R_lambda = [ cos(lambda) -sin(lambda) 0; 
            sin(lambda) cos(lambda) 0;
            0       0       1];
R_delta = [ 1          0             0;
          0 cos(delta) sin(delta);
          0 -sin(delta) cos(delta)];
  
R_theta = [cos(theta) 0 -sin(theta) ; 
            0       1         0;
          sin(theta)      0 cos(theta)];
      
R = R_rot * R_theta *  R_lambda * R_delta * R_rot' ; % for basic mesh 
R = R_rot * R_theta *  R_lambda * R_delta; % for wing mesh  

R_mir = R_rot  * R_theta * R_lambda' * R_delta';
%%
R = matlabFunction(R);
R_mir = matlabFunction(R_mir);

%%
lambda = deg2rad(30);   %sweep
theta = deg2rad(10);     %twist
delta = deg2rad(15);    %dihedral

R_right = R(delta,lambda,theta)
R_left = R_mir(delta,lambda,theta)