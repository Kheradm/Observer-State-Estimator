function  [y,x,P] = lin_sys_kalman(x0,sys,u,y_plant,P,Q,R)
% Adaptive Kalman filter for linear and linearized models 
%Get the model size psrameters
[n,m] = size(sys.b);
%Initilize all the vectors and Matrices
l = size(sys.c,1);
C = sys.c;
A = sys.a;B = sys.b;
%Predict state
y = C*x0 + sys.d*u;
x = A*x0 + B*u;
P = A*P*(A.')+ Q;
e_k = y_plant-y;
% Update state estimation
S_k = C*P*(C.')+R;
K_k = P*(C.')/S_k;
x = x + (K_k*e_k);
P = (eye(n)-K_k*C)*P;
end
