function flag = exitFlagfmincon(flag)
%EXITFLAHFMINCON - (internal function) suppress text when the optimization is successful.
if flag~= 1
    switch flag
        case 0
            message = 'All algorithms : Too many function evaluations or iterations.';
        case -1
            message = 'All algorithms : Stopped by output/plot function.';
        case -2
            message = 'All algorithms : No feasible point found.';
        case 2
            message = 'Trust-region-reflective, interior-point, and sqp : Change in X too small.';
        case 3
            message = 'Trust-region-reflective : Change in objective function too small.';
        case 4
            message = 'Active-set only : Computed search direction too small.';
        case 5
            message = 'Active-set only : Predicted change in objective function too small.';
        case -3
            message = 'Interior-point and sqp : Problem seems unbounded.';
    end
    
    disp(message);
end
end