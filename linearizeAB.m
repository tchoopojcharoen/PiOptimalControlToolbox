function [At,Bt] = linearizeAB(f,xt,ut)
%LINEARIZEAB - linearize a given dynamic system along given state and control
%trajectories
%   [AT,BT] = LINEARIZEAB(F,XT,UT) linearize a given dynamic system F along a
%   state trajectory XT and control input trajectroy UT.
%   AT and BT are function handles that take t as a paremeter
%
%   Example : 
%
%   function f = myfun(x,u,t)
%       f = [x(2)^2;-u*cos(x(1))+u^3];
%   end
%
%   xt = @(t)[t^3;t];
%   ut = @(t)sin(t);
%
%   [At,Bt] = linearizeAB(myfun,xt,ut);
%
%   % call state matrix at time t = 10
%
%   >> At(10)
%
%   ans =
%
%           0   20.000
%      -0.4498       0            

n = size(xt(0),1);
m = size(ut(0),1);

x = sym('x',[n,1]);
u = sym('u',[m,1]);
t = sym('t');

dfdx = jacobian(f(x,u,t),x);
dfdu = jacobian(f(x,u,t),u);
A = @(time)stateMatrix(dfdx,x,u,t,xt,ut,time);
B = @(time)inputMatrix(dfdu,x,u,t,xt,ut,time);
At = @(t)A(t);
Bt = @(t)B(t);

    function A = stateMatrix(dfdx,x,u,t,xt,ut,time)
        Asub = subs(dfdx,x,xt(time));
        Asub = subs(Asub,u,ut(time));
        A = eval(subs(Asub,t,time));
    end

    function B = inputMatrix(dfdu,x,u,t,xt,ut,time)
        Bsub = subs(dfdu,x,xt(time));
        Bsub = subs(Bsub,u,ut(time));
        B = eval(subs(Bsub,t,time));
    end
end