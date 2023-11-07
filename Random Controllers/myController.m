function u = myController(x)
% define the controller

x1 = x(1); % position
x2 = x(2); % velocity

%% P controller
u = -10*x1 -2*x2;

%% SMC controller
% % % % No cheating.

% c = 1;
% sigma = x2 + c*x1;
% rho = 2;
% 
% v = -rho*sign(sigma);
% u = -c*x2 + v;

end

