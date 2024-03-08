

u_rigid = vd;

% nominal 
omega0 = 457.7;
D0 = 0.2578;
u_FF_nominal = vd+ 2*D0/omega0 *ad + 1/omega0^2*jd;


% correct 
omega0 = 382.8217 ;
D0 = 0.0516;
u_FF_correct = vd+ 2*D0/omega0 *ad + 1/omega0^2*jd;


% wrong 
omega0 = 502.0145;
D0 = 0.4266;
u_FF_wrong = vd+ 2*D0/omega0 *ad + 1/omega0^2*jd;


figure
plot(ad)

figure
% hold on
plot(u_FF_nominal - vd)
hold on
plot(u_FF_correct - vd)
plot(u_FF_wrong - vd)

legend('nominal', 'correct', 'wrong')

% figure
% plot(u_FF_nominal - vd)
