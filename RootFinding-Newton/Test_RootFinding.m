clear 
close all
clc


func = @(x) x^3 + x-11  ;

x0 = -10;
maxiter = 15000;
x_tol  = 1e-20;
root_tol = 1e-20;

tic
x_root_1 = newton_secant(func, x0,maxiter, x_tol, root_tol)
toc

func(x_root_1)

tic
x_root_2 = Steffensen(func, x0,maxiter, x_tol, root_tol)
toc

func(x_root_2)



tic
x_root_00 = fzero(func,x0)
toc
func(x_root_00)
% 
% figure
% ezplot(func)

%% root finder

function x_root = newton_secant(func, x0,maxiter, x_tol, root_tol)

% p: root point,  q: function value of root point

% init
p0 = x0;
% iter = 0;

% Secant method
eps = 1e-3;

p1 = x0*(1+eps);

if  p1>=0
    p1 = p1 + eps;
else
    p1 = p1 - eps;
end

% get values
q0 = func(p0);
q1 = func(p1);

% if decrease, then update
if abs(q1) < abs(q0)
    [p0, p1, q0, q1] = deal(p1, p0, q1, q0);
end

% begin iteration
for iter = 1:maxiter

    if q1 == q0
       p = (p1 + p0) / 2;
       
       % break
       break;

    else
       
        if abs(q1) > abs(q0)
             p = (-q0 / q1 * p1 + p0) / (1 - q0 / q1);
        else
             p = (-q1 / q0 * p0 + p1) / (1 - q1 / q0);
        end

    end

    % update variable     
     p0 = p1;
     q0 = q1;

     p1 = p;
     q1 = func(p1);

     % check convergence
     if abs(p0-p1)<x_tol | abs(q1)<root_tol
         break;
     end

end

x_root = p;

end



%% Steffensen

function x_root = Steffensen(func, x0,maxiter, x_tol, root_tol)
% This function takes as inputs: a fixed point iteration function, f, 
% and initial guess to the fixed point, p0, and a tolerance, tol.
% The fixed point iteration function is assumed to be input as an
% inline function. 
% This function will calculate and return the fixed point, p, 
% that makes the expression f(x) = p true to within the desired 
% tolerance, tol.

p0 = x0;

format compact   % This shortens the output.
format long      % This prints more decimal places.

for i = 1:maxiter   % get ready to do a large, but finite, number of iterations.
                 % This is so that if the method fails to converge, we won't
                 % be stuck in an infinite loop.
    p1 = func(p0);  % calculate the next two guesses for the fixed point.
    p2 = func(p1);
    p = p0-(p1-p0)^2/(p2-2*p1+p0); % use Aitken's delta squared method to
                                  % find a better approximation to p0.

    if abs(p - p0) < x_tol | abs(p0)<root_tol  % test to see if we are within tolerance.
        break             % if we are, stop the iterations, we have our answer.
    end
    p0 = p;               % update p0 for the next iteration.
end

x_root = p;

end