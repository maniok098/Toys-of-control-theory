clear
close all
clc
%% import Casadi
import casadi.*

% init NLP variables for Casadi
w = {};       % decision variable lbw <= w <= ubw
w0 = [];      % initial guess for w
lbw = [];     % lower bound of w
ubw = [];     % upper bound of w
J = 0;        % cost function
g = {};       % constraints lbg <= g(w) <= ubg
lbg = [];     % lower bound of g(w)
ubg = [];     % upper bound of g(w)

%% construct optimal control problem

nStates = 3;
nInput = 1;

x0 = [-1.2;0;0]; % init condition of state 
xf = [2;0;0];  % end condition
u0 = 0;               

% constraints of input and states
umin = -0.3; 
umax = 0.3 ;
xmin =[-inf;-0.5;-0.3]; 
xmax = [inf;0.5;0.3];

% constraints of end time
tfmin = 1 ; 
tfmax = 100;

% number of discretization points
N = 500;

% time transformation
T = 1; % normalized time
Ts = T/N; % sampling time

%% define optimization variables

Xk = SX.sym('Xk', [nStates,N]);
Uk = SX.sym('Uk', [nInput,N-1]);
tf = SX.sym('tf');

% write decision variables w    
w = [w(:)' {Xk(:)} ];
w = [w(:)' {Uk(:)} ];
w = [w(:)' {tf} ];

% collocation
% describe system dynamics via equality constraints 
Xprevious = Xk(:,1:end-1);
Xnext = Xk(:,2:end);
Uprevious = Uk;
Xnext_mittelpunkt = mittelpunkt_vec(Xprevious,Uprevious,Xnext,Ts,tf);

% write equality constraints
% g = [g(:)' {Xnext(:) - Xnext_trapezoidal(:)}];
g = [g(:)' {Xnext(:) - Xnext_mittelpunkt(:)}];
lbg = [lbg; zeros(length(Xnext)*nStates, 1)];
ubg = [ubg; zeros(length(Xnext)*nStates, 1)];  


% write bounds of w
% bounds of Xk (3 * N) with fixed x0 and xf
lb_xk = [x0;repmat(xmin,N-2,1);xf];
ub_xk = [x0;repmat(xmax,N-2,1);xf];
% bounds of Uk 
lb_uk = repmat(umin,N-1,1);
ub_uk = repmat(umax,N-1,1);
% bounds of tf 
lb_tf = 0;
ub_tf = 1000;

% % lb and ub of w (put these together)
lbw = [lb_xk;lb_uk;lb_tf];
ubw = [ub_xk;ub_uk;ub_tf];


% write init guess of w
w0 = [x0;repmat(x0,N-2,1);xf; ...
      zeros(N-1,1);...
      10];


% cost function
gamma = 1e-1; % penalty for input changing

J = tf + gamma*sum(Uk.^2)*tf*T/N;

% create NLP problem
tic
disp('Create NLP problem in ipopt...')
prob = struct('f', J, 'x', vertcat(w{:}), 'g', vertcat(g{:}));
%  solver = nlpsol('solver', 'ipopt', prob);

% % restrict max iteration
options = struct;
options.ipopt.max_iter = 300;
% % options.ipopt.linear_solver = 'wsmp';
% % options.ipopt.mumps.ICNTL(23) = 102400; 
options.ipopt.mumps_mem_percent = 1e5;
solver = nlpsol('solver', 'ipopt', prob, options);

fprintf('time %f s\n', toc);


% solve NLP problem
tic
disp('Solve NLP problem in ipopt...')
sol = solver('x0', w0, 'lbx', lbw, 'ubx', ubw, 'lbg', lbg, 'ubg', ubg);
w_opt = full(sol.x);
fprintf('time %f s\n', toc);


%% plot results
% read results
xopt = w_opt(1:3*N);
x1 = xopt(1:3:end-3);
x2 = xopt(2:3:end-3);
x3 = xopt(3:3:end-3);
uopt = w_opt(3*N+1:3*N+N-1);
topt = w_opt(end);
fprintf('Optimal time %f s\n', topt);

% time vector for plot
dt = topt/N;
t = (1:N-1)'*dt;

figure
subplot(221)
plot(t,x1)
title('x1')
subplot(222)
plot(t,x2)
hold on
plot(t,xmin(2)*ones(size(t)),'r.-')
plot(t,xmax(2)*ones(size(t)),'r.-')
title('x2')

subplot(223)
plot(t,x3)
hold on
plot(t,xmin(3)*ones(size(t)),'r.-')
plot(t,xmax(3)*ones(size(t)),'r.-')
title('x3')
xlabel('t')

subplot(224)
plot(t,uopt)
hold on
plot(t,umin*ones(size(t)),'r.-')
plot(t,umax*ones(size(t)),'r.-')
title('u')
xlabel('t')
%% numerical integration

function Xnext_mittelpunkt = mittelpunkt_vec(Xprevious,Uprevious,Xnext,Ts,tf)
% vectorized integration rule

% 1. integration rule
% 2. transformation from free end time to fixed end time
fd = f_fun((Xprevious+Xnext)/2,Uprevious)*tf;
Xnext_mittelpunkt = Xprevious + fd*Ts;

end

function xdot = f_fun(x,u)
% define the original dynamics here

A = [0 1 0; 0 0 1; 0 0 0];
B = [0 0 1]';
xdot = A*x + B*u;

end
