function accel = CentralDifferentiation(vals, h)
% CENTRALDIFF Uses central difference method to differentiate values
% Uses forward and backward difference for first and last values
% respectively
%
% Inputs:
% vals: input values for numerical function
% h: global step size for calculations
% 
% Outputs:
% position: Vector of acceleration values indexed to timespan

    accel = zeros(length(vals),1);
    
    accel(1) = (vals(2) - vals(1)) / h; % forward diff
    accel(end) = (vals(end) - vals(end-1)) / h; %backward difff

    %central diff for all other values
    for i = 2:length(vals)-1
        accel(i) = (vals(i+1) - vals(i-1)) / (2*h);
    end
    
end