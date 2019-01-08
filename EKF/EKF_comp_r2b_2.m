function EKF_comp_r2b
clc;clear
No_Ts = 80;                         %No of Iterations=simulation time
Ts = .5;                           %Sample Time
%%%%%%%%%%%%%%%%%%%%%%Model constant parameters%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u = [-10,1];
n_y = 1;
n_x = 2;
n_u = 2;
dh_dx = [0,1];
%%%%%%%%%%%%%%%%%%%%%%Initial condition%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T_R_s = 395.326752696822;                                    %K
C_A_s = 0.573366236515890;                                   %Kmol/m^3
x0 = [C_A_s;T_R_s];
P0 = 1e3*eye(n_x);
xP0 = [1.5*x0;reshape(P0,4,1)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Noise%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mu_C_A = 0;sigma_C_A = 0.01;
mu_T = 0;sigma_T = 5;
pn1 = normrnd(mu_C_A,sigma_C_A,[No_Ts,1]);
pn2 = normrnd(mu_T,sigma_T,[No_Ts,1]);
n2 = normrnd(mu_T,sigma_T,[No_Ts,1]);
noise = 1*[n2];
pnoise = 1*[pn1,pn2];
%%%%%%%%%%%%%%%%%%%%%%Filter Parameters%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Q_Kalman = 0*eye(n_x);                %Process Noise
R_Kalman = 1e-19*eye(n_y);            %Measurment Noise
% R_Kalman = cov(noise);
for k_time=1:No_Ts
    if k_time == 40
        u = [-400,.5];
    end
    tspani = [(k_time-1)*Ts,k_time*Ts];
    [tspani, x_planti] = ode45(@(t,x) b_s_cstr(t,x,u),tspani,x0);
%     x0 = x_planti(end,:).';
    x_k(k_time,:) = x_planti(end,:);
    Tsim(k_time) =tspani(end);
    x_plant(k_time,:) = x_planti(end,:)+pnoise(k_time,:);
    x0 = x_plant(k_time,:).';
    y_plant(k_time,:) = dh_dx*x_plant(k_time,:).'+noise(k_time);
    [~, xP_pred] = ode45(@(t,xP) EKF_predict_r2b(t,xP,u,Q_Kalman),tspani,xP0);
    xP_pred = xP_pred(end,:).';
    x_pred = xP_pred(1:2);
    P_pred = reshape(xP_pred(3:6),2,2);
    K = P_pred*(dh_dx.')/(dh_dx*P_pred*dh_dx.'+R_Kalman);
    x_upd(k_time,:) = (x_pred+K*(y_plant(k_time,:).'-dh_dx*x_pred)).';
    P_upd = (eye(2)-K*dh_dx)*P_pred;
    xP0 = [x_upd(k_time,:).';reshape(P_upd,4,1)];
    k_time
end
x_upd = [x0.';x_upd];
x_plant = [x0.';x_plant];
x_k = [x0.';x_k];
Tsim = [0,Tsim];
subplot(2,1,1)
plot(Tsim,x_upd(:,1),Tsim,x_plant(:,1),'--g',Tsim,x_k(:,1),'--r')
subplot(2,1,2)
plot(Tsim,x_upd(:,2),Tsim,x_plant(:,2),'--g',Tsim,x_k(:,2),'--r')
clc
end