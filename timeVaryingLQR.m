function [St,Kt] = timeVaryingLQR(A,B,Q,R,Qf,tf,N)
%TIMEVARYINGLQR - solve time-varying linear quadratic regulator problem by
%solving the corresponding Differential Algebraic Ricatti Equation.
%   [ST,KT] = TIMEVARYINGLQR(A,B,Q,R,Qf,TF,N) solve time-varying linear
%   quadratic regulator problem by a Differential Algebraic
%   Ricatti Equation (DARE) based on give linearized state matrix A ,
%   input matrix B, matrices Q, R, Qf, a final TF and number of samples N.
%   The given DARE formulation is used for solving an optimal control
%   problem that minimize the following cost J.
%
%   J(u) = x(tf)'*Qf*x(tf) + int(x(t)'*Q(t)*x(t)+u(t)'*R(t)*u(t), t=0,tf)
%
%   with given linearized dynamic system
%
%   x'(t) = A(t)*x(t)+B(t)*u(t)
%
%   The DARE can written in the following form.
%
%   S'(t) = S(t)*B(t)*K(t) - S(t)*A(t) - A(t)'*S(t) - Q(t)
%
%   where K(t) = inv(R(t))*B'*S(t)
%         S(tf) = Qf
%
%   A,B,Q, and R can be matrices or function handles that take time (t) as
%   an argument, and return the matrices.
%
%   The function returns function handles of S(t) and K(t).

S = DAREsolver(A,B,Q,R,Qf,tf,N);
[St,Kt] = piecewiseCubicS(S,A,B,Q,R,tf,N);
end

%% helper functions

function [S] = DAREsolver(A,B,Q,R,Qf,tf,N)
%DARESOLVER - solves a given Differential Algebraic Ricatti Equation using
%direct collocation method.
%   [S,K] = DARESOLVER(A,B,Q,R,Qf,TF,N) solves a Differential Algebraic
%   Ricatti Equation (DARE) based on give linearized state matrix A ,
%   input matrix B, matrices Q, R, Qf, a final TF and number of samples N.
%   The given DARE formulation is used for solving an optimal control
%   problem that minimize the following cost J.
%
%   J(u) = x(tf)'*Qf*x(tf) + int(x(t)'*Q(t)*x(t)+u(t)'*R(t)*u(t), t=0,tf)
%
%   with given linearized dynamic system
%
%   x'(t) = A(t)*x(t)+B(t)*u(t)
%
%   The DARE can written in the following form.
%
%   S'(t) = S(t)*B*K(t) - S(t)*A - A'*S(t) - Q
%
%   where K(t) = inv(R)*B'*S(t)
%         S(tf) = Qf
%
%   The function returns discretize matrix S, and feedback matrix K

%% Validate attributes :

% TO DO : write the validation for each argument

%% Determine time dependency

if ~isa(A,'function_handle')
    A = @(t)A;
end

if ~isa(B,'function_handle')
    B = @(t)B;
end

if ~isa(Q,'function_handle')
    Q = @(t)Q;
end

if ~isa(R,'function_handle')
    R = @(t)R;
end

%% Guess for w
n = size(Qf,1);
w0 = zeros(n*n*N,1); % each point along te curve (not including collocation)

%% Formulate problem for fmincon

FUN = @(w) 0;
As = [];
Bs = [];
Aeq = [];
Beq = [];
LB = [];
UB = [];
NONLCON = @(w) constraint(w,A,B,Q,R,Qf,tf,N);

OPTIONS = optimoptions('fmincon','Algorithm','sqp','display','off');

%% Optimize using fmincon (or other method)

[w,~,flag] = fmincon(FUN,w0,As,Bs,Aeq,Beq,LB,UB,NONLCON,OPTIONS);
exitFlagfmincon(flag);

% obtain S
wMatrix = permute(reshape(w,n,n,N),[2 1 3]);
S = cell(1,N+1);
for i = 1:N,
    S{i} = wMatrix(:,:,i);
end
S{N+1} = Qf;
%
% % fix this part for matrix
% trajectory.X = piecewiseCubic(f,X,U,tf,N);
% trajectory.U = piecewiseLinear(U,tf,N);
t = linspace(0,tf,N+1);
for i = 1:N
    
    Si = S{i};
    Sn = S{i+1};
    subplot(4,1,1);
    hold on;
    plot([t(i) t(i+1)],[Si(1,1) Sn(1,1)])
    
    subplot(4,1,2);
    hold on;
    plot([t(i) t(i+1)],[Si(1,2) Sn(1,2)])
    
    subplot(4,1,3);
    hold on;
    plot([t(i) t(i+1)],[Si(2,1) Sn(2,1)])
    
    subplot(4,1,4);
    hold on;
    plot([t(i) t(i+1)],[Si(2,2) Sn(2,2)])
    
    
end
end

%% help functions

function [inCon,eqCon] = constraint(w,A,B,Q,R,Qf,tf,N)
%CONSTRAINT - returns evaluated collocation constraint
%   [INC,EQC] = CONSTRAINT(W,A,B,Q,R,Qf,TF,N) evaluates
%   and forms collocation constraints
%   one from the final condition + the rest is from collocation method

%% Obtain guess S
n = size(Qf,1);
wMatrix = permute(reshape(w,n,n,N),[2 1 3]);
S = cell(1,N+1);
for i = 1:N,
    S{i} = wMatrix(:,:,i);
end
S{N+1} = Qf;

%% Evaluate constraints

eqCon = collocationConstraint(S,A,B,Q,R,tf,N);
inCon = [];
end

function colCon = collocationConstraint(S,A,B,Q,R,tf,N)
%COLLOCATIONCONSTRAINT - internal function : evaluate collocation
%constraints

% T*c = b;
% c = T\b;

n = size(S{1},1);
dt = tf/N;
t = linspace(0,tf,N+1);
T = [   eye(n*n)    zeros(n*n)  zeros(n*n)      zeros(n*n)      ;
    eye(n*n)    eye(n*n)    eye(n*n)        eye(n*n)        ;
    zeros(n*n)  eye(n*n)/dt zeros(n*n)      zeros(n*n)      ;
    zeros(n*n)  eye(n*n)/dt eye(n*n)*2/dt   eye(n*n)*3/dt   ];

colCon = zeros(n*n*N,1); % collocation cosntraint vector at each point i

for i = 1:N,
    b = zeros(4*n*n,1);
    b(1:n*n) = reshape(permute(S{i},[2,1]),n*n,1);
    b(n*n+1:2*n*n) = reshape(permute(S{i+1},[2,1]),n*n,1);
    b(2*n*n+1:3*n*n) = diffS(S{i},A(t(i)),B(t(i)),Q(t(i)),R(t(i))); % nested function diffS(Si,Ai,Bi,Qi,Ri)
    b(3*n*n+1:end) = diffS(S{i+1},A(t(i+1)),B(t(i+1)),Q(t(i+1)),R(t(i+1)));
    
    ci = T\b;
    %   coefficient = permute(reshape(reshape(ci,4,n*n),n,n,4),[2 1 3 4];
    coefficient = reshape(ci,n*n,4);
    tCol = (t(i)+t(i+1))/2;
    tDiff = tCol-t(i);
    d = [1;tDiff/dt;(tDiff/dt)^2;(tDiff/dt)^3];
    dDot = [0;1/dt;2*tDiff/dt^2;3*tDiff^2/dt^3];
    
    SCol = permute(reshape(coefficient*d,n,n),[2,1]);
    SdotCol = permute(reshape(coefficient*dDot,n,n),[2,1]);
    
    colCon((i-1)*n*n+1:i*n*n) = diffS(SCol,A(tCol),B(tCol),Q(tCol),R(tCol))-SdotCol;
    
end

end


function [St,Kt] = piecewiseCubicS(S,A,B,Q,R,tf,N)
%PIECEWISECUBICS - create a piecewise cubic S mstrix for algebraic riccati
%equation


S = @(t)piecewise(t,S,A,B,Q,R,tf,N);
St = @(t)S(t);
Kt = @(t)(R(t)\(B(t)'*S(t)));
    function st = piecewise(time,S,A,B,Q,R,tf,N)
        %COLLOCATIONCONSTRAINT - internal function : evaluate collocation
        %constraints
        
        % T*c = b;
        % c = T\b;
        
        n = size(S{1},1);
        dt = tf/N;
        t = linspace(0,tf,N+1);
        T = [   eye(n*n)    zeros(n*n)  zeros(n*n)      zeros(n*n)      ;
            eye(n*n)    eye(n*n)    eye(n*n)        eye(n*n)        ;
            zeros(n*n)  eye(n*n)/dt zeros(n*n)      zeros(n*n)      ;
            zeros(n*n)  eye(n*n)/dt eye(n*n)*2/dt   eye(n*n)*3/dt   ];
        
        
        coefficient = cell(1,N);
        for i = 1:N,
            b = zeros(4*n*n,1);
            b(1:n*n) = reshape(permute(S{i},[2,1]),n*n,1);
            b(n*n+1:2*n*n) = reshape(permute(S{i+1},[2,1]),n*n,1);
            b(2*n*n+1:3*n*n) = diffS(S{i},A(t(i)),B(t(i)),Q(t(i)),R(t(i))); % nested function diffS(Si,Ai,Bi,Qi,Ri)
            b(3*n*n+1:end) = diffS(S{i+1},A(t(i+1)),B(t(i+1)),Q(t(i+1)),R(t(i+1)));
            
            ci = T\b;
            %   coefficient = permute(reshape(reshape(ci,4,n*n),n,n,4),[2 1 3 4];
            coefficient{i} = reshape(ci,n*n,4);
            
        end
        
        if numel(time)>1
            if size(time,1)<= size(time,2) % row vector of time
                st = cell(1,length(time));
                for i = 1:length(time)
                    if time(i) <= 0
                        st{i} = S{1};
                    elseif time(i)>=tf
                        st{i} = S{N};
                    else
                        idx = floor(time(i)/tf*N)+1;
                        dt = tf/N;
                        tDiff = time(i)-dt*(idx-1);
                        d = [1;tDiff/dt;(tDiff/dt)^2;(tDiff/dt)^3];
                        st{i} = permute(reshape(coefficient{idx}*d,2,2),[2 1]);
                    end
                end
            end
        else
            if time <= 0
                st = S{1};
            elseif time>=tf
                st = S{N};
            else
                idx = floor(time/tf*N)+1;
                dt = tf/N;
                tDiff = time-dt*(idx-1);
                d = [1;tDiff/dt;(tDiff/dt)^2;(tDiff/dt)^3];
                st = permute(reshape(coefficient{idx}*d,2,2),[2 1]);
            end
        end
    end
end

function dS = diffS(S,A,B,Q,R)
dS = S*B*(R\(B'*S))-S*A-A'*S-Q;
end
