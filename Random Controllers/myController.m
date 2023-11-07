function u = myController(x,x_des)
% define the controller

e1 = x(1) - x_des(1); % position error
e2 = x(2) - x_des(2); % velocity error

%% P controller
% u = -10*e1 -2*e2;

%% SMC controller
% % % % No cheating.

c = 5;
sigma = e2 + c*e1;
rho = 4;

v = -rho*sign(sigma);
u = -c*e2 + v;

end

