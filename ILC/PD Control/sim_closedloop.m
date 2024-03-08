%% init
clear
close all
clc
% get reference signal
load('traj_7pp.mat')
% time = traj_7pp.Time;
traj_KGT = traj_7pp.Variables;
traj_KGT = [traj_KGT;zeros(500,4)];

Ts = 1/4e3; % sample time
time = Ts*(0:length(traj_KGT)-1)';



xdes = traj_KGT(:,1);
vdes = traj_KGT(:,2);
ades = traj_KGT(:,3);
jdes = traj_KGT(:,4);

% model parameters
% two-mass no fric
m1 = 142.4748;
m2 = 470.4535;
k = 5.9679e7;
d = 1.1731e4;

% friction
fc1 = 15.9184;
fv1 = 303.4944;
fc2 = 64.4997;
fv2 = 383.2807;

% transmission ratio
gamma = 0.04/(2*pi);
beta = 1e3;
theta = [m1 m2 k d fc1 fv1 fc2 fv2 gamma beta];
clear m1 m2 k d fc1 fv1 fc2 fv2 gamma beta





%% Simulation with PPI and FF
x_sim = zeros(4,length(time));

% P position controller
Kv = 50;

% PI velocity controller
KpV = 180;
KiV = KpV*23*Ts;

% C_pi = pid(KpV,KiV,0,0,Ts);

% torque feedforward
lin_acc_2_torque = 3.724433384330698;  % equivalent "mass"  % a_table to tau_m

% init u_PI
e_vm_previous = 0;
u_PI_previous = 0;

for k = 1:length(time)-1
  
  % read init condition
  xk = x_sim(:,k);
  
  % outer loop
  vm_des = vdes(k) + Kv*(xdes(k)-xk(2));
  
  % motor velocity error
  e_vm = vm_des - xk(3);
  % PI controller
  u_PI = u_PI_previous + KpV  *(e_vm - e_vm_previous) + KiV*e_vm;

  % control action
  uk = ades(k)*lin_acc_2_torque + u_PI;
  
  % one step ahead
  x_sim(:,k+1) = lpred_RK4(xk,uk,theta,Ts);
%   x_sim(:,k+1) = lpred_ode45(xk,uk,theta,Ts);
  
  % update for PI 
  e_vm_previous =  e_vm;
  u_PI_previous =  u_PI;

end

figure
subplot(211)
plot(time,xdes,time,x_sim(1,:))
ylabel('Position')
title('Closed loop control FF+FB')
legend('x_{des}','x_{act}')
subplot(212)
plot(time,xdes'-x_sim(1,:))
title('Position error')

%% PD controller

m_all = 500;

x_sim = zeros(4,length(time));

% P position controller
% Kv = 50;

Kp = 10;
Kd = 0.5;

gamma = 0.04/2/pi;

for k = 1:length(time)-1
  
  % read init condition
  xk = x_sim(:,k);
  
  % table position
  xT_act = xk(2);
  vT_act = xk(4);
  
  uFF = m_all*ades(k)*gamma^2;
  uFB_k = m_all* ( Kp* (xdes(k) - xT_act) + Kd*(vdes(k) - vT_act) );

  uk = uFF +uFB_k;
  
  % one step ahead
  x_sim(:,k+1) = lpred_RK4(xk,uk,theta,Ts);
%   x_sim(:,k+1) = lpred_ode45(xk,uk,theta,Ts);
  
  % update for PI 
  e_vm_previous =  e_vm;
  u_PI_previous =  u_PI;

end

figure
subplot(211)
plot(time,xdes,time,x_sim(1,:))
ylabel('Position')
title('Closed loop control PD control')
legend('x_{des}','x_{act}')
subplot(212)
plot(time,xdes'-x_sim(1,:))
title('Position error')

%% ILC


m_all = 500;





% P position controller
% Kv = 50;

Kp = 2;
Kd = 0.5;

gamma = 0.04/2/pi;

nIter = 5;

% init FF
uFF = m_all*ades*gamma^2;

uFB = zeros(nIter,length(uFF));

figure(10) 

for iter = 1:nIter

x_sim = zeros(4,length(time));

    for k = 1:length(time)-1
      
      % read init condition
      xk = x_sim(:,k);
      
      % table position
      xT_act = xk(2);
      vT_act = xk(4);
      
%       uFF = m_all*ades(k)*gamma^2;
      uFB_k = m_all* ( Kp* (xdes(k) - xT_act) + Kd*(vdes(k) - vT_act) );
      
      uFB(iter,k) = uFB_k;

      uk = uFF(k) + uFB_k;
      
      % one step ahead
      x_sim(:,k+1) = lpred_RK4(xk,uk,theta,Ts);
    %   x_sim(:,k+1) = lpred_ode45(xk,uk,theta,Ts);
      
      % update for PI 
      e_vm_previous =  e_vm;
      u_PI_previous =  u_PI;
    
    end

    % update feedforward
    eta = 0.5;
%     uFF(1:end-1) = uFF(1:end-1) + eta * uFB(iter,2:end)' ; 

    uFF(1:end-3) = uFF(1:end-3) + eta *( uFB(iter,2:end-2)' ...
                                       + 0.1*  uFB(iter,3:end-1)'...
                                       + 0.05 * uFB(iter,4:end)') ;    

    figure(10)
    subplot(211)
    plot(time,xdes,time,x_sim(1,:))
    hold on
    ylabel('Position')
    title('Closed loop control PD control')
    legend('x_{des}','x_{act}')
    subplot(212)
    plot(time,xdes'-x_sim(1,:))
    hold on
    title('Position error')



end

% figure
% subplot(211)
% plot(time,xdes,time,x_sim(1,:))
% ylabel('Position')
% title('Closed loop control PD control')
% legend('x_{des}','x_{act}')
% subplot(212)
% plot(time,xdes'-x_sim(1,:))
% title('Position error')



%% Prediction function

function x_next = lpred_RK4(xk,uk,theta,Ts)

% xk: x(k)
% uk: u(k)
% theta: model parameters
% Ts: sampling time

% selection dynamics function
f_fun = @f_twomass_fric;

s1 = f_fun(xk,uk,theta);
s2 = f_fun(xk+0.5*Ts*s1,uk,theta);
s3 = f_fun(xk+0.5*Ts*s2,uk,theta);
s4 = f_fun(xk+s3*Ts,uk,theta);

% predicted x(k+1)
x_next = xk + (1/6)*(s1+2*s2+2*s3+s4)*Ts;

end

function x_next = lpred_ode45(xk,uk,theta,Ts)

% xk: x(k)
% uk: u(k)
% theta: model parameters
% Ts: sampling time

% selection dynamics function
f_fun = @f_twomass_fric;

% setup solver
odeopts = odeset('RelTol',1e-2);  
integrator = @ode45;

% Gridded data interpolation, to get the function
% u(t)
griddedU = griddedInterpolant([0 Ts],[uk uk],'nearest');
% griddedU = uk; % zero-order hold

% define system ode function 
SystemODEFun = @(t,x)f_fun(x,griddedU(t),theta);
% simulation time
tspan = [0,Ts/2,Ts];

% simulate
[~,x_ode] = integrator(SystemODEFun, tspan, xk, odeopts);

% predicted x(k+1)
x_next = x_ode(end,:)';

end
