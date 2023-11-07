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

nStates = 2;
nInput = 1;

x0 = [1;0]; % init condition of state 
u0 = 0;               

% % constraints of input and states
umin = -inf; 
umax =  inf ;
xmin =[-inf;-inf]; 
xmax = [inf;inf];

% % constraints of end time
% tfmin = 1 ; 
% tfmax = 100;

% number of discretization points
N = 100;

T = 5;
Ts = T/N;


%% define optimization variables

Xk = SX.sym('Xk', [nStates,N]);
Uk = SX.sym('Uk', [nInput,N-1]);

% write decision variables w    
w = [w(:)' {Xk(:)} ];
w = [w(:)' {Uk(:)} ];

% collocation
% describe system dynamics via equality constraints 
Xprevious = Xk(:,1:end-1);
Xnext = Xk(:,2:end);
Uprevious = Uk;
Xnext_mittelpunkt = mittelpunkt_vec(Xprevious,Uprevious,Xnext,Ts);

% write equality constraints
% g = [g(:)' {Xnext(:) - Xnext_trapezoidal(:)}];
g = [g(:)' {Xnext(:) - Xnext_mittelpunkt(:)}];
lbg = [lbg; zeros(length(Xnext)*nStates, 1)];
ubg = [ubg; zeros(length(Xnext)*nStates, 1)];  

% equality constraints of end state
%  x2(tf) - x1(tf) = 1
g = [g(:)' {Xk(2,end) - Xk(1,end)}];
lbg = [lbg; 1];
ubg = [ubg; 1];

% write bounds of w
% bounds of x0
lb_xk = [x0;repmat(xmin,N-1,1)];
ub_xk = [x0;repmat(xmax,N-1,1)];
% bounds of Uk 
lb_uk = repmat(umin,N-1,1);
ub_uk = repmat(umax,N-1,1);


% % lb and ub of w (put these together)
lbw = [lb_xk;lb_uk];
ubw = [ub_xk;ub_uk];


% write init guess of w
w0 = [x0;repmat(x0,N-1,1); ...
      zeros(N-1,1)];

% cost function
J = sum(Uk.^2)+ sum( Xk(:).^2) ;

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
xopt = w_opt(1:2*N);
x1 = xopt(1:2:end-2);
x2 = xopt(2:2:end-2);
uopt = w_opt(2*N+1:2*N+N-1);

t = (1:N-1)'*Ts;

figure
subplot(221)
plot(t,x1)
title('x1')
subplot(222)
plot(t,x2)
title('x2')
subplot(223)
plot(t,uopt)
title('u')

%% numerical integration

function Xnext_mittelpunkt = mittelpunkt_vec(Xprevious,Uprevious,Xnext,Ts)
% vectorized integration rule

%  integration rule
fd = f_fun((Xprevious+Xnext)/2,Uprevious);
Xnext_mittelpunkt = Xprevious + fd*Ts;

end

function xdot = f_fun(x,u)
% define the original dynamics here
%  vectorized
x1 = x(1,:);
x2 = x(2,:);

dx1 = x2;
dx2 = -x1 + (1-x1.^2).*x2 +u;

xdot = [dx1;dx2];

end

