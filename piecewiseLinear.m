function trajectory = piecewiseLinear(X,tf,N)
%PIECEWISELINEAR - create a piecewise linear trajectory based on given points
%   T = PIECEWISELINEAR(X,TF,N) creates a piecewise linear trajectory T
%   based on given discretized trajectory X, final time, TF, and number of
%   samples N. The trajectory T is a function handle that accept time as 
%   its parameter.
%
%   Example :
%   traj = piecewiseLinear(X,tf,N);  % generates trajectory
%   value = traj(0.5);               % evaluate trajectory at time t = 0.5

trajectory = @(t)piecewise(t,X,tf,N); 
trajectory = @(t)trajectory(t);

    function xt = piecewise(t,X,tf,N)
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
                    
                    if (idx==N) && (size(X,2)<=N)
                        xt(:,i) = X(:,idx);
                    else
                        xt(:,i) = X(:,idx)+(t(i)-dt*(idx-1))*(X(:,idx+1)-X(:,idx))/dt;
                    end
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
                    
                    if (idx==N) && (size(X,2)<=N)
                        xt(i,:) = X(:,idx)';
                    else
                        xt(i,:) = (X(:,idx)+(t(i)-dt*(idx-1))*(X(:,idx+1)-X(:,idx))/dt)';
                    end
                end
            end
        end
    end
end