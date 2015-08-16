function h = plotVector(x,y,varargin)
%PLOTVECTOR - plot a set of trajectory of given vector on different plots
%   PLOTVECTOR(X,Y) plots a set of trajectory from given 
%   pairs (X,Y) with default x-label and y-label. Then plots them on 
%   seperate plots in the given figure handle H
%
%   PLOTVECTOR(X,Y,LABELX,LABELY) plots a set of trajectory from given 
%   pairs (X,Y) with x-label LABELX and y-label LABELY. Then plots them on 
%   seperate plots in the given figure handle H
%
%   H = PLOTVECTOR(X,Y,LABELX,LABELY,H) also return the figure handle.

validateattributes(x,{'numeric'},{'vector'});
validateattributes(y,{'numeric'},{});

if size(x,2)~= size(y,2)
    error('x and y must have same length')
end

h = gcf;
n = size(y,1);

if numel(varargin)
    title = varargin{1};
    labelX = varargin{2};
    labelY = varargin{3};
else
    title = 'y_i[n] vs. n';
    labelX = 'n \{sample\}';
    labelY = cell(1,n);
    for i = 1:n,
        labelY{i} = sprintf('y_{%d}[n]',i);
    end
end

validateattributes(title,{'char'},{});
validateattributes(labelX,{'char'},{});
validateattributes(labelY,{'cell'},{});

for i = 1:numel(labelY)
    validateattributes(labelY{i},{'char'},{});
end

color = {};
suptitle(title);

for i = 1:n,
    ax = subplot(n,1,i);
    
    color_i = rand(1,3);
    if ~isempty(color)
        exist = false;
    else
        exist = false;
        for j = 1:numel(color)
            exist = exist || norm(color{j}-color_i)<=0.05;
        end
    end
    
    invisible = norm(color_i)>=1.0392;
    k = 1;
    while exist || invisible,
        if ~(mod(k,100))
            disp('finding colors for plots ...');
        end
        color_i = rand(1,3);
        if isempty(color)
            exist = false;
        else
            exist = false;
            for j = 1:numel(color)
                exist = exist || (norm(color{j}-color_i)<=0.5);
            end
        end
        invisible = norm(color_i)>=1.0392;
        k = k+1;
    end
    color{end+1} = color_i;
    plot(ax,x,y(i,:),'color',color_i);
    
    xlabel(labelX);
    ylabel(labelY{i});
    
end


end