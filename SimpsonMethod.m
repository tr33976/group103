function distTravelled = SimpsonMethod(endTime, timeSpan, velocity, h)
% SimpsonMethod Uses the simpson's rule for numerical integration to
% numerically evaluate definite integral up to the end time range. 0 is
% assumed to be the beginning of the range.
%
% Inputs:
% endTime: End of integral range
% timeSpan: Discrete time range to integrate over
% velocity: Discrete range of velocity values
% h: Global step size
%
% Outputs:
% distTravelled: Estimated value for numerical integral up to end time. 

    %limit n to within required timespan
    n = length(timeSpan(timeSpan<=endTime));
    so = 0;
    for i = 1:n / 2
        so = so + abs(velocity(2*i - 1)); % sum all odd indices
    end
    se = 0;
    for i = 1:(n / 2) - 1
        se = se + abs(velocity(2*i)); % sum all even indices
    end
    distTravelled = h/3 * (abs(velocity(1)) + 4 * so + 2*se + abs(velocity(n)));
end
