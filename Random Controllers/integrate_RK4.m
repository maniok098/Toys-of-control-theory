function x_next = integrate_RK4(t,x_previous,u_previous,Ts)
% one-step-ahead integration of differential equation 
%   method: 4th-order Runge-Kutta 

% Inputs:
%   t         :   current simulation time
%   x_previous:   x(k), state vector at time k
%   u_previous:   u(k), input at time k
%   Ts        :   sampling time

% Output:
%   x_next:   x(k+1), state vector at time k+1

% Auxiliary function:
%   f_fun: continuous differential equation with xdot = f_fun(x,u,theta)

s1 = f_fun(t,x_previous,             u_previous);
s2 = f_fun(t,x_previous + 0.5*Ts*s1, u_previous);
s3 = f_fun(t,x_previous + 0.5*Ts*s2, u_previous);
s4 = f_fun(t,x_previous + s3*Ts,     u_previous);
x_next = x_previous + (1/6)*(s1+2*s2+2*s3+s4)*Ts;
end



