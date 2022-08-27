clc
clear
close all

%% init
% initial condition
%  e.g. postion = -10,  velocity = 1
x0 = [-10;
        1];

% sampling time of simulation
Ts = 0.001;
% maximum simulation time
t_end = 100;
% time
time = 0:Ts:t_end;
% number of samples
nSamples = length(time);
% init x_sim
x_sim = [x0, zeros(2,nSamples-1)];
%% design desired trajectories
[x_des, v_des] = getTrajectories(time);
x_state_desired = [x_des; v_des];  % store in one matrix

%% simulation
for idx = 1:length(time)-1 
    % read t, state vector, and desired state
    t_previous = time(idx);
    x_previous = x_sim(:,idx);
    x_des_previous = x_state_desired(:,idx);
    % design input 
    u_previous =  myController(x_previous,x_des_previous);
    % one-step-ahead integration
    x_next = integrate_RK4(t_previous,x_previous,u_previous,Ts);
    x_sim(:,idx+1) = x_next;
end

%% plot result
figure
ax(1) = subplot(211);
plot(time,x_sim(1,:))
hold on
plot(time,x_des)
ylabel('position [m]')
legend('actual','desired')
title('Result of your controller')

ax(2) = subplot(212);
plot(time,x_sim(2,:))
hold on
plot(time,v_des)
legend('actual','desired')
xlabel('time [s]')
ylabel('velocity [m/s]')

linkaxes(ax,'x');



