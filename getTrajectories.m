function [x_des, v_des] = getTrajectories(time)
% get desired trajectories


time = time(:)'; % convert to row vector

% x_des : desired position 
% v_des : desired velocity 


%  e.g. 
a0 = 1;
omega = 3;

x_des = a0*sin(omega*time);
v_des = a0*omega*cos(omega*time);

end

