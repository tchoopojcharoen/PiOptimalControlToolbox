function xNew = eulerMethod(f,x,u,t,dt)
%EULERMETHOD - updates the next iteration of given dynamic system using 
%Euler's method
%   XNEW = EULERMETHOD(F,X,U,DT) updates the next iteration of a given
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

xNew = x + dt*f(x,u,t);

end