%% Example 1: directSingleShooting

clear all;
clc;

% functions

lagrange = @(x,u,t)(u^2);
mayer = @(x,tf)(0);
k = 1;
b = 1;
mass = 1;
f = @(x,u,t)([x(2) ; -(k/mass)*x(1)-(b/mass)*x(2)+(1/mass)*u]);

% constraints

eqPathCon = @(x,u)(0);
inPathCon = @(x,u)([x(1)-7 ; -x(1)-7]);
eqTerCon = @(xN,tf)(xN-[0;0]);
inTerCon = @(xN,tf)(0);

x0 = [5;0];
tf = 2;
N = 20;
m = 1; % number of control input signal

method = 'rk4';
[X,U,t,J] = directSingleShooting(lagrange,mayer,f,eqPathCon,inPathCon,eqTerCon,inTerCon,x0,tf,N,m,method);
