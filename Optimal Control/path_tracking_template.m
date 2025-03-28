% function path_tracking_template()
% PATH_TRACKING_CVX - Time-energy optimal path tracking
%                     numerical solution by SOCP using CVX
%
%   Verscheure, D. and Demeulenaere, B. and Swevers, J.
%   and De Schutter, J. and Diehl, M.:
%   Time-Optimal Path Tracking for Robots:
%   a Convex Optimization Approach
%   IEEE Trans. on Automatic Control 54(2009)10,2318-2327
%
% E. Arnold   2016-06-28
%             2019-11-01 CVX
% M. Mrochen  2019-12-16 Template Trajectory Generation
clear all
close all
clc
%% Install CVX
% if 0
%     cwd = pwd;
%     addpath(genpath('cvx-w64'));
%     cvx_setup
% end % if
% 
% %% Add CVX to Matlab path
% if ~exist('cvx_clear')
%     cwd = pwd;
%     addpath(genpath('cvx-w64'));
%     cvx_startup
%     cd(cwd)
% end
% 
% % Delete result files and figures
% try
%     delete('path_tracking_res.dat', 'path_tracking_resx.dat')
%     delete([1 2 3 4 10 21])
%     drawnow
% catch
% end

%% Parameters
N = 100;            % number of discretization intervals
gamma1 = 0;       % weighting coefficient tau^2
gamma2 = 1e-6;      % weighting coefficient tau'
m_x = 1;
m_y = 1;
F_x_minmax = 0.1;
F_y_minmax = 0.25;
A_x = 1.0;
A_y = 1.0;
om_x = 2.0*pi;
om_y = 4.0*pi;
phi_x = pi/2.0;
phi_y = pi/4.0;

%% Discretization grid
ds = 1.0/N;
s = (0:N-1)'*ds+ds/2.0;

%% Build m(s), c(s), g(s) from a)
q_x = A_x*sin(om_x*s + phi_x);
q_y = A_y*sin(om_y*s + phi_y);
q_x_s = A_x*om_x*cos(om_x*s + phi_x);
q_y_s = A_y*om_y*cos(om_y*s + phi_y);
q_x_ss = -A_x*om_x^2*sin(om_x*s + phi_x);
q_y_ss = -A_y*om_y^2*sin(om_y*s + phi_y);
ms = [m_x*q_x_s   m_y*q_y_s];
cs = [m_x*q_x_ss  m_y*q_y_ss];
gs = zeros(length(s), 2);

%% Solve optimization problem
cvx_quiet false
tic
[b, F_x, F_y, s, dt, t] = solve_problem(ms, cs, gs, F_x_minmax, F_y_minmax, gamma1, gamma2, N);
toc
%% Plot solution
% Path velocity s_dot
% ..

% Forces F_x, F_y
% ..
figure
plot(t,[0;F_x], t, [0;F_y])

% Trajectory q_y vs. q_x
% ..
kk = 1
%% Write solution file and start Simulink simulation
% dlmwrite('path_tracking_res.dat', [t(1) F_x(1) F_y(1); t(1:N)+dt/2.0 F_x F_y; t(N+1) F_x(N) F_y(N)], ' ')
% dlmwrite('path_tracking_resx.dat', [t A_x*sin(om_x*s+phi_x) A_y*sin(om_y*s+phi_y) A_x*om_x*cos(om_x*s+phi_x).*sqrt(b) A_y*om_y*cos(om_y*s+phi_y).*sqrt(b)], ' ')
% 
% % Load parameters
% tF = textread('path_tracking_res.dat');
% tx = textread('path_tracking_resx.dat');
% A_x = 1.0;
% om_x = 2.0*pi;
% phi_x = pi/2.0;
% A_y = 1.0;
% om_y = 4.0*pi;
% phi_y = pi/4.0;

% Simulation and plot without controller
% ..

% Solve problem of deviation of actual and reference
% ..


% Function solve_problem
% Setup and solve SOCP
function [b, F_x, F_y, s, dt, t] = solve_problem(ms, cs, gs, F_x_minmax, F_y_minmax, gamma1, gamma2, N)
    ds = 1.0/N;
%     sgamma1 = sqrt(gamma1);
    cvx_begin % cvx specification
        % solver settings
        cvx_solver sedumi
        cvx_precision low
        
        % optimization variables
        variable a(N) % [0:N-1]
        variable b(N+1) nonnegative % [0:N]
        variable c(N+1)
        variable d(N)
        variable e(N-1,2) nonnegative
        variable F_x(N) % [0:N-1]
        variable F_y(N) % [0:N-1]

        % objective
        minimize (2*ds*norm(d,1) + gamma2*norm(e,1)) % (2.12 a)
        
        % constraints
        subject to
            % (2.12c)-(2.12e)
            b(1) == 0
            % (2.12g)
            b(N+1) == 0
            % (2.12b)
            b(2:N+1) - b(1:N) == 2*ds*a
            % (2.12h)
            -F_x_minmax <= F_x <= F_x_minmax
            -F_y_minmax <= F_y <= F_y_minmax
            % (2.12i)-% (2.12j)
            F_x == ms(:,1).*a + cs(:,1).*(b(1:N) + b(2:N+1))/2 + gs(:,1)
            F_y == ms(:,2).*a + cs(:,2).*(b(1:N) + b(2:N+1))/2 + gs(:,2)
            % (2.12f) declared with 'nonnegative'
            % (2.12h)
            -e <= [F_x(2:N)-F_x(1:N-1)  F_y(2:N)-F_y(1:N-1)] <= e
            % (2.12i)
%             for k= 1:N
%                 norm([2; 2*sqrt(gamma1)*F_x(k); ...
%                          2*sqrt(gamma1)*F_y(k); ...
%                          c(k+1)+c(k)-d(k)], 2)  ... 
%                      <=  c(k+1)+c(k)+d(k)
%             end
%             %  (2.12j)
%             for k= 1:N+1
%                norm([2*c(k); b(k)-1] , 2) <= b(k)+1
%             end

%      % vectorized

            norms([2*ones(N,1),...
                   2*sqrt(gamma1)*F_x, ...
                   2*sqrt(gamma1)*F_y, ...
                   c(2:N+1)+c(1:N)-d(1:N)], 2, 2)  ... 
                 <=c(2:N+1)+c(1:N)+d(1:N)
            
            norms([2*c, b-1],2,2) <= b+1


            
    cvx_end  % call solver
    if strcmp(cvx_status, 'Solved')
        s = (0:N)'*ds;
        dt = 2.0*ds./(sqrt(b(1:N))+sqrt(b(2:N+1)));
        t = [0; cumsum(dt)];
    else
        fprintf('%s\n', cvx_status)
        b = []; F_x = []; F_y = []; s = []; dt = []; t = [];
    end
end
