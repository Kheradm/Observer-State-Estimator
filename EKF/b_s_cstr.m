function xdot = b_s_cstr(~,x,u)
%%%Enhanced Stability Regions for Model Predictive Control of Nonlinear Process Systems
%%%Maaz Mahmood and Prashant Mhaskar
%%%%%%%%%%%%%%%%%%%%%%Table 2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
V = 0.1;                                                     %m^3
R = 8.314;                                                   %KJKmol^-1.K
C_A0_s = 1.0;                                                %Kmol/m^3
T_A0_s = 310.0;                                              %K
Q_s = 0.0;                                                   %KJ/min
NDH = 4.78e4;                                                %KJ/Kmol
k_0 = 72e9;                                                  %1/min
E = 8.314e4;                                                 %KJ/Kmol
c_p = 0.239;                                                 %KJKg/K
rho = 1e3;                                                   %kg/m^3
F = 100e-3;                                                  %m^3/min
T_R_s = 395.33;                                              %K
C_A_s = 0.57;                                                %Kmol/m^3
%%%%%%%%%%%%%%%%%%%%%%%Unstable equilibrium point%%%%%%%%%%%%%%%%%%%%%%%%%%
% T_R_s = 395.326752696822;                                    %K
% C_A_s = 0.573366236515890;                                   %Kmol/m^3
%%%%%%%%%%%%%%%%%%%%%%%Assign Inputs and States%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Q = u(1);C_A0 = u(2);
C_A = x(1);T_R = x(2);
%%%%%%%%%%%%%%%%%%%%  Equations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C_A_dot = F/V*(C_A0-C_A)-k_0*(exp(-E/(R*T_R)))*C_A;
T_R_dot = F/V*(T_A0_s-T_R)+NDH/(rho*c_p)*k_0*exp(-E/(R*T_R))*C_A+Q/(rho*c_p*V);
xdot =[C_A_dot;T_R_dot];
end