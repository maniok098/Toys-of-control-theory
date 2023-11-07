clc
clear
close all

% initial condition
%  e.g. postion = -10,  velocity = 1
x0 = [-10;
        1];

% sampling time of simulation
Ts = 0.001;
% maximum simulation time
t_end = 50;
% time
time = 0:Ts:t_end;
% number of samples
nSamples = length(time);
% init x_sim
x_sim = [x0, zeros(2,nSamples-1)];
% simulation
for idx = 1:length(time)-1 
    t = time(idx);
    x_previous = x_sim(:,idx);
    u_previous =  myController(x_previous);
    x_next = integrate_RK4(t,x_previous,u_previous,Ts);
    x_sim(:,idx+1) = x_next;
end

% plot result
figure
plot(time,x_sim(1,:))
hold on
plot(time,x_sim(2,:))
legend('position [m]','velocity [m/s]')
xlabel('time [s]')
title('Result of your controller')


