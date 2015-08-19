function [X,t,h] = forwardSimulation(f,x0,U,tf,N,intMethod,varargin)
%FORWARDSIMULATION - forward simulate an IVP with given dynamic system and
%specified method at equal time step
%   [X,T] = FORWARDSIMULATION(F,X0,U,TF,N,INT) simulates a controlled
%   dynamic system F forward in time with given initial states X0,
%   a (m X N) control input trajectory, a final time TF, number of iterations N
%   , and an integration method INT. The function return state trajectory X
%   and time vector T.
%
%   [X,T] = FORWARDSIMULATION(F,X0,U,TF,N,'given',MET) simulates a controlled
%   dynamic system with the given integration method MET.
%
%   [X,T,H] = FORWARDSIMULATION(F,X0,U,TF,N,...,'plot') simulates the
%   controlled dynamic system, plot state and control input
%   trajectories, and return state trajectory X, time vector T, and plot
%   handle H.

%% validate attributes

validateattributes(f,{'function_handle'},{});
validateattributes(x0,{'numeric'},{'column','finite'});
validateattributes(U,{'numeric'},{'finite'});
validateattributes(tf,{'numeric'},{'scalar','finite'});
validateattributes(N,{'numeric'},{'integer','finite'});

% verify method's validity
listNotMATLABode = {'euler','runge-kutta','rk4','given'};

isMATLABode = true;
for i = 1:numel(listNotMATLABode)
    isMATLABode = isMATLABode && ~strcmpi(intMethod,listNotMATLABode{i});
end

if strcmpi(intMethod,'euler')
    update = @(f,x,u,t,dt)eulerMethod(f,x,u,t,dt);
elseif strcmpi(intMethod,'runge-kutta')||strcmpi(intMethod,'rk4')
    update = @(f,x,u,t,dt)rungeKuttaMethod(f,x,u,t,dt);
elseif strcmpi(intMethod,'given')
    if numel(varargin) == 2
        update = varargin{1};
    elseif numel(varargin) == 1
        if strcmpi(varargin{:},'plot')
            error('method is not given');
        else
            validateattributes(varargin{:},{'function_handle'},{});
            update = varargin{:};
        end
    elseif numel(varargin) == 0
        error('method is not given');
    else
        error('too many parameters');
    end
else
    error('invalid update method');
end

if size(U,2)~=N
    error('Given control input trajectory must have length of N');
end


t = linspace(0,tf,N+1);
dt = t(2)-t(1);

if (size(x0,1)~=size(update(f,x0,U(:,1),t(1),dt),1))||(size(x0,2)~=size(update(f,x0,U(:,1),t(1),dt),2))
    error('Update method must return vector f the same size of x0');
end

%% iterations
n = size(x0,1);
m = size(U,1);

X = zeros(n,N+1);
X(:,1) = x0;

for k = 1:N,
    X(:,k+1) = update(f,X(:,k),U(:,k),t(:,k),dt);
end

%% Visualization

% obatin plotting command
if numel(varargin) == 2
    if strcmpi(varargin{2},'plot')
        doPlot = true;
    else
        error('invalid plot command')
    end
elseif numel(varargin) == 1
    if strcmpi(intMethod,'given')
        doPlot = false;
    else
        if strcmpi(varargin{1},'plot')
            doPlot = true;
        elseif isempty(varargin{1})
            doPlot = false;
        else
            error('invalid plot command');
        end
    end
elseif numel(varargin) == 0
    doPlot = false;
else
    error('too many input parameters');
end


if doPlot
    h = plotXU(t,X,U);
else
    h = {};
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