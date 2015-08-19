function trajectory = piecewiseCubic(f,X,U,tf,N)
%PIECEWISECUBIC - create a piecewise cubic trajectory based on given state
%and control input trajectories
%   T = PIECEWISECUBIC(F,X,U,TF,N) creates a piecewise cubic trajectory T
%   based on given dynamic system F, discretized state trajectory X, 
%   discretized control input trajectory U , final time, TF, and number of
%   samples N. The trajectory T is a function handle that accept time as 
%   its parameter.
%
%   Example :
%   traj = piecewiseCubic(f,X,U,tf,N);  % generates trajectory
%   value = traj(0.5);               % evaluate trajectory at time t = 0.5

trajectory = @(t)piecewise(t,f,X,U,tf,N);
trajectory = @(t)trajectory(t);
    function xt = piecewise(t,f,X,U,tf,N)
        
        % Compute coefficients
        n = size(X,1);
        dt = tf/N;
        timeVector = linspace(0,tf,N+1);
        A = [   eye(n)      zeros(n)    zeros(n)    zeros(n)    ;
            eye(n)      eye(n)      eye(n)      eye(n)      ;
            zeros(n)    eye(n)/dt   zeros(n)    zeros(n)    ;
            zeros(n)    eye(n)/dt   eye(n)*2/dt eye(n)*3/dt ];
        
        coefficient = cell(1,N);
        for i = 1:N,
            b = zeros(4*n,1);
            b(1:n) = X(:,i);
            b(n+1:2*n) = X(:,i+1);
            b(2*n+1:3*n) = f(X(:,i),U(:,i),timeVector(i));
            if i ~= N
                b(3*n+1:end) = f(X(:,i+1),U(:,i+1),timeVector(i+1));
            else
                b(3*n+1:end) = f(X(:,i+1),U(:,i),timeVector(i+1));
            end
            
            ci = A\b;
            coefficient{i} = reshape(ci,n,4);
        end
        
        
        if size(t,1)<= size(t,2) % row vector of time
            xt = zeros(size(X,1),length(t));
            for i = 1:length(t)
                if t(i) <= 0
                    xt(:,i) = X(:,1);
                elseif t(i)>=tf
                    xt(:,i) = X(:,end);
                else
                    idx = floor(t(i)/tf*N)+1;
                    dt = tf/N;
                    tDiff = t(i)-dt*(idx-1);
                    d = [1;tDiff/dt;(tDiff/dt)^2;(tDiff/dt)^3];
                    
                    xt(:,i) = coefficient{idx}*d;
                    
                end
            end
        else
            xt = zeros(length(t),size(X,1));
            for i = 1:length(t)
                if t(i) <= 0
                    xt(i,:) = X(:,1)';
                elseif t(i)>=tf
                    xt(i,:) = X(:,end)';
                else
                    idx = floor(t(i)/tf*N)+1;
                    dt = tf/N;
                    tDiff = t(i)-dt*(idx-1);
                    d = [1;tDiff/dt;(tDiff/dt)^2;(tDiff/dt)^3];
                    
                    xt(i,:) = (coefficient{idx}*d)';
                    
                end
            end
        end
    end
end