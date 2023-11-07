function xdot = f_fun(t,x,u)
% define the system dynamics via differential equation
%   xdot = f_fun(x,u,theta)

% Inputs:
%   t         :   current simulation time
%   x         :   state vector
%   u         :   input (applied force)
%   theta     :   parameters of the model

%% underwater vehicle

x1 = x(1); % position 
x2 = x(2); % velocity

umax= 10;
% saturation of u
if u>umax
    u = umax;
elseif u<-umax
    u = -umax;
else 
    % do nothing
end

k = 1;     % damping ratio
m = 1;     % mass

disturbance = sin(t);  % disturbance force
% disturbance = 0;  % disturbance force

dx1 = x2;  % dot_x_1, velocity
dx2 = ( -k*abs(x2)*x2 + u + disturbance ) / m; 

xdot = [dx1;
        dx2];

end

