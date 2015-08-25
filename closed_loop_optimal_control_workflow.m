%{
 This is an example of a complete workflow for an optimal control feedback 
 policy design. 
 
 I used a fully-actuated pendulum as an example.
%}

%% clear workspace

clear all;
clc;
elaspedTime = 0;
%% define optimal control problem

% functions
lagrange = @(x,u,t)(u'*u);
mayer = @(x,tf)(0);
f = @(x,u,t)([x(2);-sin(x(1))+u]);

% constraints
uMax = 10;
eqPathCon = @(x,u)(0);
inPathCon = @(x,u)(u*[1;-1]-uMax*[1;1]);
eqTerCon = @(xN,tf)(xN-[0;0]);
inTerCon = @(xN,tf)(0);

x0 = [-pi;0];
tf = 10;
N = 40;
m = 1; % number of control input signal
method = 'rk4';

%% Obtain optimal state and control trajectories
tic;
disp('solving optimal control problem using fmincon');

[~,~,~,trajectory,~] = directCollocation(lagrange,mayer,f,eqPathCon,inPathCon,eqTerCon,inTerCon,x0,tf,N,m,method);

endTime = toc;
str = sprintf('Optimal trajectories are obtained. Took %.3f seconds',endTime);
disp(str);
elaspedTime = elaspedTime + endTime;

%% Linearize nonlinear system to get A(t) and B(t)
tic;
disp('linearizing the system along the obtained optimal trajectories');

[At,Bt] = linearizeAB(f,trajectory.X,trajectory.U);

endTime = toc;
str = sprintf('Linearization is completed. Took %.3f seconds',endTime);
disp(str);
elaspedTime = elaspedTime + endTime;

%% Define Time-Varying LQR problem and obtain close-loop feedback gain, K(t)

Q = eye(2);
Qf = eye(2);
R = 1;
Q = @(t)Q;
R = @(t)R;

tic;
disp('calculating for S(t) and K(t) for time-varying LQR');

[St,Kt] = timeVaryingLQR(At,Bt,Q,R,Qf,tf,10);

endTime = toc;
str = sprintf('S(t) and K(t) are calculated. Took %.3f seconds',endTime);
disp(str);
elaspedTime = elaspedTime + endTime;

%% Simulate for verification

%%

xt = trajectory.X;
ut = trajectory.U;
u = @(x,t)(ut(t)-Kt(t)*(x-xt(t)));

tic;
disp('Simulating the closed-loop control system')

[T,X] = ode45(@(t,x)f(x,u(x,t),t),[0 tf],x0+[0.01 -0.01]');

endTime = toc;
str = sprintf('Time T and state trajectory X are obtained. Took %.3f seconds',endTime);
disp(str);
elaspedTime = elaspedTime + endTime;

str = sprintf('It took %f minutes in total to run this script.',elaspedTime/60);
disp(str);
%% visualization

step = 1;
for i = 1:step:numel(T)-step
    plot([0 sin(X(i,1))],[0 cos(X(i,1))]);
    axis([-1 1 -1 1]);
    axis equal;
    pause(T(i+step)-T(i))
end
