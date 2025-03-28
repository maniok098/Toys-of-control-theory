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

%% def
nStates = 2;  %  pos and vel
nInput = 1;   %  acc

% boundary conditions
% x0

% constraints of input and states
%   kinematic limitation of Ax1-x
umin_1 = -10; 
umax_1 =  10 ;
xmin_1 = [-0.3;
          -0.2]; 
xmax_1 = [0.3;
          0.2];
%   kinematic limitation of Ax1-y
umin_2 = -10; 
umax_2 =  10 ;
xmin_2 = [-0.3;
          -0.2]; 
xmax_2 = [0.3;
          0.2];


%   kinematic limitation of Ax2-x
umin_3 = -10; 
umax_3 =  10 ;
xmin_3 = [-0.01;
          -2]; 
xmax_3 = [0.01;
            2];
%   kinematic limitation of Ax2-y
umin_4 = -10; 
umax_4 =  10 ;
xmin_4 = [-0.01;
          -2]; 
xmax_4 = [0.01;
           2];

% number of points to be considered
N = 1000;
% sampling time
Ts = 0.01;

%% define optimization variables
Xk_Ax1_x = SX.sym('Xk_Ax1_x', [nStates,N]);
Uk_Ax1_x = SX.sym('Uk_Ax1_x', [nInput,N-1]);

Xk_Ax1_y = SX.sym('Xk_Ax1_y', [nStates,N]);
Uk_Ax1_y = SX.sym('Uk_Ax1_y', [nInput,N-1]);

Xk_Ax2_x = SX.sym('Xk_Ax2_x', [nStates,N]);
Uk_Ax2_x = SX.sym('Uk_Ax2_x', [nInput,N-1]);

Xk_Ax2_y = SX.sym('Xk_Ax2_y', [nStates,N]);
Uk_Ax2_y = SX.sym('Uk_Ax2_y', [nInput,N-1]);

% write decision variables w    
w = [w(:)' {Xk_Ax1_x(:)} ];
w = [w(:)' {Uk_Ax1_x(:)} ];
w = [w(:)' {Xk_Ax1_y(:)} ];
w = [w(:)' {Uk_Ax1_y(:)} ];
w = [w(:)' {Xk_Ax2_x(:)} ];
w = [w(:)' {Uk_Ax2_x(:)} ];
w = [w(:)' {Xk_Ax2_y(:)} ];
w = [w(:)' {Uk_Ax2_y(:)} ];

% system dynamics
Ad = [1 0.01;
      0    1];
Bd = [5e-5;
      0.01];

% describe system dynamics as equlity constraints
% Ax1 - x
Xprevious = Xk_Ax1_x(:,1:end-1);
Xnext = Xk_Ax1_x(:,2:end);
Uprevious = Uk_Ax1_x;
Xnext_Analytic = Ad*Xprevious + Bd*Uprevious;
g = [g(:)' {Xnext(:) - Xnext_Analytic(:)}];
lbg = [lbg; zeros(length(Xnext)*nStates, 1)];
ubg = [ubg; zeros(length(Xnext)*nStates, 1)];  

% Ax1 - y
Xprevious = Xk_Ax1_y(:,1:end-1);
Xnext = Xk_Ax1_y(:,2:end);
Uprevious = Uk_Ax1_y;
Xnext_Analytic = Ad*Xprevious + Bd*Uprevious;
g = [g(:)' {Xnext(:) - Xnext_Analytic(:)}];
lbg = [lbg; zeros(length(Xnext)*nStates, 1)];
ubg = [ubg; zeros(length(Xnext)*nStates, 1)];  

% Ax2 - x
Xprevious = Xk_Ax2_x(:,1:end-1);
Xnext = Xk_Ax2_x(:,2:end);
Uprevious = Uk_Ax2_x;
Xnext_Analytic = Ad*Xprevious + Bd*Uprevious;
g = [g(:)' {Xnext(:) - Xnext_Analytic(:)}];
lbg = [lbg; zeros(length(Xnext)*nStates, 1)];
ubg = [ubg; zeros(length(Xnext)*nStates, 1)];  

% Ax2 - y
Xprevious = Xk_Ax2_y(:,1:end-1);
Xnext = Xk_Ax2_y(:,2:end);
Uprevious = Uk_Ax2_y;
Xnext_Analytic = Ad*Xprevious + Bd*Uprevious;
g = [g(:)' {Xnext(:) - Xnext_Analytic(:)}];
lbg = [lbg; zeros(length(Xnext)*nStates, 1)];
ubg = [ubg; zeros(length(Xnext)*nStates, 1)];  

% write bounds of w
% bounds of Xk (nSates * N) 
lb_xk = [x0_1;repmat(xmin_1,N-2,1);xf_1];
ub_xk = [x0_1;repmat(xmax_1,N-2,1);xf_1];
% bounds of Uk 
lb_uk = repmat(umin_1,N-1,1);
ub_uk = repmat(umax_1,N-1,1);
% % lb and ub of w (put these together)
lbw = [lb_xk;lb_uk];
ubw = [ub_xk;ub_uk];
