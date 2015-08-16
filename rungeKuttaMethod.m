function xNew = rungeKuttaMethod(f,x,u,t,dt)
%RUNGEKUTTAMETHOD - updates the next iteration of given dynamic system 
%using 4th-Order Runge-Kutta's method
%   XNEW = RUNGEKUTTAMETHOD(F,X,U,DT) updates the next iteration of a given
%   function handle of dynamic system F, a previous state X, and a control 
%   input U with the time step of DT.

% validate attributes
validateattributes(f,{'function_handle'},{});
validateattributes(x,{'numeric'},{'column','finite'});
validateattributes(u,{'numeric'},{'column','finite'});
validateattributes(t,{'numeric'},{'scalar','finite'});
validateattributes(dt,{'numeric'},{'scalar','finite'});
if (size(f(x,u,t),1)~= size(x,1))||(size(f(x,u,t),2)~= size(x,2))
    error('dimension of f(.) must be the same as x');
end

% update

k1 = f(x,u,t);
k2 = f(x + k1*dt/2,u,t + dt/2);
k3 = f(x + k2*dt/2,u,t + dt/2);
k4 = f(x + k3*dt,u,t + dt);

xNew = x + (k1+2*k2+2*k3+k4)*dt/6;

end