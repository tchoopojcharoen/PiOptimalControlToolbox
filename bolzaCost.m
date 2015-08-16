function [B, intL,M] = bolzaCost(lagrange,mayer,X,U,tf,intMethod,varargin)
%BOLZACOST - calculates the cost based on given Bolza's formulation, state
%trajectory, and control input trajectory over a specified amount of time.
%   B = BOLZACOST(LAGRANGE,MAYER,X,U,TF,INT) calculates the cost based 
%   on the given Bolza's formulation 
%
%   B = int(L(x(t),u(t),t),t = 0:T) + M(x(T),T)
%
%   which includes the following
%   
%   Langrange's term (L)                : LAGRANGE - @(x,u,t)langrange -> scalar
%   Mayer's term (M)                    : MAYER    - @(xN,tf)mayer -> scalar
%   state trajectory (x(t=0:T))         : X        - n x(N+1) matrix
%   control input trajectory (u(t=0:T)) : U        - m x N matrix
%   final time (T)                      : TF       - positive finite scaler
%   integration method (int)            : INT      - 'euler' or 'rk4'
%
%   B = BOLZACOST(LAGRANGE,MAYER,X,U,TF,'given',MET) calculates the cost based 
%   on the given Bolza's formulation and integration method MET in the 
%   form met(f,x,u,t).
%
%   [B,INTL] = BOLZACOST(...) returns the total cost B and integrated
%   Lagrage's term INTL
%
%   [B,INTL,M] = BOLZACOST(...) also returns evaluated Mayer's term M

%% Augment Langrange's term to perform forward simulation

n = size(X,1);
N = size(U,2);
    function lAug = langrangeAug(lagrange,~,uAug,t)
        x = uAug(1:n);
        u = uAug(n+1:end);
        lAug = lagrange(x,u,t);
    end
UAug = [X(:,1:end-1);U];
augmentedLagrange = @(l,uAug,t)langrangeAug(lagrange,l,uAug,t);

%% Forward Simualtion

[S,~] = forwardSimulation(augmentedLagrange,0,UAug,tf,N,intMethod,varargin{:});

%% evaluating total cost

intL = S(end);          % integrated Langrange's term
M = mayer(X(:,end),tf); % Mayer's term

B = intL+M;             % Bolza's form  

end