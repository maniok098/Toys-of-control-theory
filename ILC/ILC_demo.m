clear 
close all
clc

%% get reference signal

load('traj_7pp.mat')
time = traj_7pp.Time;
traj_KGT = traj_7pp.Variables;

Ts = 1/4e3; % sample time
%% construct nominal model

omega0 = 215.1584;
D0 = 0.44;

A = [0 1 0;
    0 0 1;
    0 -omega0^2 -2*D0*omega0];
B = [0 0 omega0^2]';
C = [1 0 0];

Sys_nominal = ss(A,B,C,0);
Sys_nominal_discrete = c2d(Sys_nominal,Ts);

Ad = Sys_nominal_discrete.A;
Bd = Sys_nominal_discrete.B;

%% read xd, vd
xd = traj_KGT(:,1);
vd = traj_KGT(:,2);
ad = traj_KGT(:,3);
jd = traj_KGT(:,4);

N = length(xd);

nIter = 5;

u_FF = vd;% + 2*D0/omega0 *ad + 1/omega0^2*jd;

u_FB = zeros(nIter,N);
err_k = zeros(nIter,N);

%% low pass filter for Q-filter
% cut-off freq
fc = 20*2*pi;
% fc.Maximum = 300;
Tc = 1/fc; 
% differentiator
% s = tf('s');
% PT3 for synchronization
PT1 = tf(1,[Tc,1]);%tf(1,[Tc^3,3*Tc^2,3*Tc,1]);%tf(1,[Tc,1]);
PT1_discrete = c2d(PT1,0.001)

b = PT1_discrete.Numerator{1,1};
a = PT1_discrete.Denominator{1,1};

% xdes =  [zeros(1,100), ones(1,100)];
% xdes_filter = filter(b, a, xdes);
% figure
% plot(xdes)
% hold on
% plot(xdes_filter)

%% closed loop simulation
figure

for iter = 1:nIter

    x_sim = zeros(3,length(time));
    % x_sim(:,1) = [0 0 0]';
    
    % P position controller
    Kv = 50;
    
    for k = 1:length(time)-1
      
      % read init condition
      x0 = x_sim(:,k);
      % read control action

      u = u_FF(k) + Kv*(xd(k)-x0(1));
%       u = vd(k) + Kv*(xd(k)-x0(1));

      % save feedback
      u_FB(iter,k) = Kv*(xd(k)-x0(1));
    
      % one step ahead
      x_sim(:,k+1) = Ad*x0 + Bd*u; 
    
    end


    % update feedforward
    gamma = 0.2;
    % causal learning
%     u_FF(1:end-1) = u_FF(1:end-1) + gamma * u_FB(iter,2:end)' ; 

    % non-causal learning
%     u_FF(1:end-2) = u_FF(1:end-2) + gamma *( u_FB(iter,2:end-1)'+ 0.1*  u_FB(iter,3:end)') ; 
    u_FF(1:end-3) = u_FF(1:end-3) + gamma *( u_FB(iter,2:end-2)' ...
                                           + 0.1*  u_FB(iter,3:end-1)'...
                                           + 0.05 * u_FB(iter,4:end)') ;    
% add Q-filter (bad)
%     u_FF = filter(b, a, u_FF);  

    plot(time,xd'-x_sim(1,:))
    hold on
end

legend('1','2','3','4','5')