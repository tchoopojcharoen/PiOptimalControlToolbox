function [X,U,t,trajectory,J] = directCollocation(lagrange,mayer,f,eqPathCon,inPathCon,eqTerCon,inTerCon,x0,tf,N,m,intMethod,varargin)
%DIRECTCOLLOCATION - solves a given optimal control problem using direct
%collocation method
%   [X,U,T,J] = DIRECTCOLLOCATION(LAGRANGE,MAYER,F,EQPC,INPC,EQTC,INTC,X0,TF,N,M,INT)
%   solves a continuous time optimal control problem based on given Bolza's
%   formulation,
%
%      J* = min { int(L(x(t),u(t),t),t = 0:T) + M(x(T),T) }
%           u(.)
%
%   u*(t) = argmin { int(L(x(t),u(t),t),t = 0:T) + M(x(T),T) }
%           u(.)
%
%   s.t.           x'(t)  =  f(x(t),u(t),t)     : dynamic constraint
%           g(x(t),u(t))  =  0                  : equality path constraint
%           h(x(t),u(t)) <=  0                  : inequality path constraint
%              q(x(T),T)  =  0                  : equality terminal constraint
%              r(x(T),T) <=  0                  : inequality terminal constraint
%                   x(0)  =  x0                 : initial condition
%
%   ,which includes the following arguments.
%
%   Langrange's term (L)          : LAGRANGE - @langrange(x,u,t) -> scalar
%   Mayer's term (M)              : MAYER    - @mayer(xN,tf) -> scalar
%   dynamic system (f)            : F        - @f(x,u,t) -> dim(x)x1 column vector
%   eq. path constraint (g)       : EQPC     - @eqPathCon(x,u) -> num_eq_p column vector
%   ineq. path constraint (h)     : INPC     - @inPathCon(x,u) -> num_in_p column vector
%   eq. terminal constraint (q)   : EQTC     - @eqTerCon(xN,tf) -> num_eq_t column vector
%   ineq. terminal constraint (r) : INTC     - @inTerCon(xN,tf) -> num_in_t column vector
%   initial state (x0)            : X0       - n x 1 column vector
%   final time(T)                 : TF       - positive finite scaler or empty
%   number of sample              : N        - positive integer
%   number of control input       : m        - positive integer
%   integration method (int)      : INT      - 'euler' or 'rk4'
%
%   The function returns and plot state (x(t)) and control input (u(t))
%   trajectories, X and U respectively, with respect to a time vector T.
%   The function also returns the minimum cost (J*) J.
%
%   [X,U,T] = DIRECTCOLLOCATION(LAGRANGE,MAYER,F,EQPC,INPC,EQTC,INTC,X0,TF,N,M,'given',MET)
%   solves a continuous time optimal control problem based on given Bolza's
%   formulation with the given integration menthod MET.

%% Validate attributes :

% TO DO : write the validation for each argument

%% Guess for w
n = size(x0,1);
w0 = rand(n*N+m*N,1);
freeTF = isempty(tf);

if freeTF
    w0 = [w0;10];
end

%% Formulate problem for fmincon

FUN = @(w) bolza(w,lagrange,mayer,x0,tf,N,intMethod,varargin);
A = [];
B = [];
Aeq = [];
Beq = [];
LB = [];
UB = [];
NONLCON = @(w) constraint(w,f,eqPathCon,inPathCon,eqTerCon,inTerCon,x0,tf,N,intMethod,varargin);

OPTIONS = optimoptions('fmincon','Algorithm','sqp','display','off');

%% Optimize using fmincon (or other method)

[w,J,flag] = fmincon(FUN,w0,A,B,Aeq,Beq,LB,UB,NONLCON,OPTIONS);
exitFlagfmincon(flag);

%% Obtain and visulaize state and control input trajectories

if isempty(tf)
    tf = w(end);
    u_vector = w(n*N+1:end-1);
else
    u_vector = w(n*N+1:end);
end
x_vector = w(1:n*N);

t = linspace(0,tf,N+1);
X = [x0 reshape(x_vector,length(x_vector)/N,N)];
U = reshape(u_vector,length(u_vector)/N,N);
plotXU(t,X,U);

% [coefficient,XCol,UCol] = collocationPoints(f,X,U,tf,N);

trajectory.X = piecewiseCubic(f,X,U,tf,N);
trajectory.U = piecewiseLinear(U,tf,N);

end

%% help functions

function b = bolza(w,lagrange,mayer,x0,tf,N,intMethod,varargin)
%BOLZA - calculates the total cost based on given Bolza's formulation
%   B = BOLZA(W,L,M,X0,TF,N,INT) calculates the total cost based on given
%   guess W, Lagrange's term L, Mayer's term M, an
%   initial states X0, a final time TF, and number of sample N with the
%   integration method INT.
%
%   B = BOLZA(W,L,M,F,X0,TF,N,'given',MET) calculates the total cost based
%   on given Bolza's formulation with given integration method MET.

% Obtain guess X and U

n = size(x0,1);
if isempty(tf)
    tf = w(end);
    u_vector = w(n*N+1:end-1);
else
    u_vector = w(n*N+1:end);
end

x_vector = [x0 ; w(1:n*N)];
XGuess = reshape(x_vector,n,N+1);
UGuess = reshape(u_vector,length(u_vector)/N,N);

% Numerically evaluate Bolza's formulation using the guess

[b,~,~] = bolzaCost(lagrange,mayer,XGuess,UGuess,tf,intMethod,varargin{:});

end

function [inCon,eqCon] = constraint(w,f,eqPathCon,inPathCon,eqTerCon,inTerCon,x0,tf,N,~,varargin)
%CONSTRAINT - returns evaluated inequality and equality constaints
%   [INC,EQC] = CONSTRAINT(W,F,X0,EQPC,INPC,EQTC,INTC,TF,N,INT) evaluates
%   and forms inequality and equality constraints , INC and EQC respectively,
%   based on given guess W, dyanmic system F, an initial states X0,
%   equality path constraint EQPC, inequality path constraint INPC,
%   equality terminal constraint EQTC, inequality terminal constraint INTC,
%   final time TF, number of sample N, and integration method INT.
%
%   [INC,EQC] = CONSTRAINT(W,F,X0,EQPC,INPC,EQTC,INTC,TF,N,'given',MET) evaluates
%   and forms inequality and equality constraints using the given
%   integration method MET for forward simulation of state trajectory.

%% Obtain guess X and U
n = size(x0,1);
freeTF = isempty(tf);
if freeTF
    tf = w(end);
    u_vector = w(n*N+1:end-1);
    dimTF = 1;
else
    u_vector = w(n*N+1:end);
    dimTF = 0;
end
x_vector = [x0 ; w(1:n*N)];
n = length(x_vector)/(N+1);
m = length(u_vector)/N;
XGuess = reshape(x_vector,n,N+1);
UGuess = reshape(u_vector,m,N);

%% Evaluate constraints

% Allocate space for equality constraint

dimEqPathCon = size(eqPathCon(XGuess(:,1),UGuess(:,1)),1);
dimEqTerCon = size(eqTerCon(XGuess(:,end),tf),1);
dimDynCon = n*N;
eqCon = zeros(dimEqPathCon*(N+1)+dimEqTerCon+dimDynCon,1);

% Allocate space for inequality constraint

dimInPathCon = size(inPathCon(XGuess(:,1),UGuess(:,1)),1);
dimInTerCon = size(inTerCon(XGuess(:,end),tf),1);
inCon = zeros(dimInPathCon*(N+1)+dimInTerCon+dimTF,1);

% Path Constraint
for i = 1:N,
    eqCon(dimEqPathCon*(i-1)+1:dimEqPathCon*i,1) = eqPathCon(XGuess(:,i),UGuess(:,i));
    inCon(dimInPathCon*(i-1)+1:dimInPathCon*i,1) = inPathCon(XGuess(:,i),UGuess(:,i));
end
eqCon(dimEqPathCon*N+1:dimEqPathCon*(N+1),1) = eqPathCon(XGuess(:,N+1),UGuess(:,N));
colCon = collocationConstraint(f,XGuess,UGuess,tf,N);
eqCon(dimEqPathCon*(N+1)+dimEqTerCon + 1:end) = colCon;

% Terminal Constraint

eqCon(dimEqPathCon*(N+1)+1:dimEqPathCon*(N+1)+dimEqTerCon,1) = eqTerCon(XGuess(:,end),tf);
inCon(dimInPathCon*N+1:dimInPathCon*N+dimInTerCon,1) = inTerCon(XGuess(:,end),tf);

% Positive Time Constraint

if freeTF
    inCon(dimInPathCon*(N+1)+dimInTerCon+1,1) = -tf;
end

end

function colCon = collocationConstraint(f,X,U,tf,N)
%COLLOCATIONCONSTRAINT - internal function : evaluate collocation
%constraints

n = size(X,1);
dt = tf/N;
t = linspace(0,tf,N+1);
A = [   eye(n)      zeros(n)    zeros(n)    zeros(n)    ;
    eye(n)      eye(n)      eye(n)      eye(n)      ;
    zeros(n)    eye(n)/dt   zeros(n)    zeros(n)    ;
    zeros(n)    eye(n)/dt   eye(n)*2/dt eye(n)*3/dt ];


colCon = zeros(n*N,1);

for i = 1:N,
    b = zeros(4*n,1);
    b(1:n) = X(:,i);
    b(n+1:2*n) = X(:,i+1);
    b(2*n+1:3*n) = f(X(:,i),U(:,i),t(i));
    if i ~= N
        b(3*n+1:end) = f(X(:,i+1),U(:,i+1),t(i+1));
    else
        b(3*n+1:end) = f(X(:,i+1),U(:,i),t(i+1));
    end
    
    ci = A\b;
    coefficient = reshape(ci,n,4);
    tCol = (t(i)+t(i+1))/2;
    tDiff = tCol-t(i);
    d = [1;tDiff/dt;(tDiff/dt)^2;(tDiff/dt)^3];
    dDot = [0;1/dt;2*tDiff/dt^2;3*tDiff^2/dt^3];
    
    xCol = coefficient*d;
    xdotCol = coefficient*dDot;
    if i ~=N
        uCol = U(:,i)+(tCol-t(i))*(U(:,i+1)-U(:,i))/dt;
    else
        uCol = U(:,i)+(tCol-t(i))*(U(:,i)-U(:,i))/dt;
    end
    
    colCon((i-1)*n+1:i*n) = f(xCol,uCol,tCol)-xdotCol;
    
end

end

function h = plotXU(t,X,U)
%PLOTXU - internal function : plots state and control input trajectories

n = size(X,1);
m = size(U,1);
scrsz = get(groot,'ScreenSize');

% State Trajectory
figure('Position',[50 scrsz(4)/5 scrsz(3)/2-50 3*scrsz(4)/5],'Name','State Trajectory');

labelY = cell(1,n);
for i = 1:n,
    labelY{i} = sprintf('x_{%d}',i);
end
h1 = plotVector(t,X,'State Trajectory : x(t)','t : time',labelY);

% Control Input Trajectory
figure('Position',[scrsz(3)/2 scrsz(4)/5 scrsz(3)/2-50 3*scrsz(4)/5],'Name','Control Input Trajectory');

labelY = cell(1,m);
for i = 1:m,
    labelY{i} = sprintf('u_{%d}',i);
end
h2 = plotVector(t(1:end-1),U,'Control Input Trajectory','t : time',labelY);

h(1) = h1;
h(2) = h2;
end

