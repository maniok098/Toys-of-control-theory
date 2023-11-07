clear 
clc

A = [0 1;
     0 0];

B = [0;
     1];

C = eye(2);
D = [0;0];

sysc = ss(A,B,C,D)
Ts = 0.01;
sysd = c2d(sysc,0.01)