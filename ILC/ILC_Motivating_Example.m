close all 
clear all
clc

A = [-0.8 -0.22; 1 0];
B = [0.5; 1];
C = [1 0.5];
D = 0;

T = 100;
k = 1:1:T;
r = sin(0.8*k/10);

x0 = [0; 0];
u = zeros(1,T);
iterations = 10;

gamma_vec = [0.5 0.8 1.1 1.5];


for i = 1:4
    
    gamma = gamma_vec(i)    
    
    x = zeros(2,T);
    u = zeros(1,T);
    y = zeros(1,T);

    for t = 1:iterations

        % Simulate data
        for k=1:T
            x(:,k+1) = A*x(:,k) + B*u(1,k);
            y(1,k) = C*x(:,k)+ D*u(1,k);
        end

        % ILC Update
        e = r-y;
        u = u + gamma*[e(2:end) 0];
        
        % Plot result
        figure(1)
        subplot(2,2,i)
        plot(y)
        xlabel('k')
        ylabel('y')
        ylim([-1.1 1.1])
        title(['\gamma =',num2str(gamma)])
        hold on
        
    end

end