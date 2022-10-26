function [len, spring, time, acc] = WaterTouchParams(persHeight, jumpheight)
% WaterTouchParams Iterate values for rope length and spring elasticity to
% find optimal pair that satisfies a water touch experience and its
% restrcitions. 

% Input:
% persHeight: Height of jumper
% jumpheight: Total height of jump from river to platform

% Output:
% len: Optimal rope length
% spring: Optimal spring elasticity
% time: Jump time for optimal parameters
% acc: Max G's for optimal paramters

h = 0.01;
timeSpan = 0:h:100;

A=jumpheight-persHeight;
B=A+0.1;

ropeLengths = 30:0.5:70;
springCoefs = 50:0.5:100;

results = zeros((length(ropeLengths)*length(springCoefs)), 4);

%%%%%% DEFAULT PARAMETERS minus sweep vals %%%
c = 0.9; % drag coefficient (kg/m)
m = 80; % jumper mass (kg)
C = c/m; % drag / mass 
g = 9.8; % gravity (m/s^2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



idx = 1;
f = waitbar(0,'Searching water touch params');
for L = ropeLengths
    for k = springCoefs
        % velocity ode
        dvdt = @(y, v) g - C .* abs(v) .* v - max(0, k/m .*(y-L));
        [position, velocity] = RK4Coupled(dvdt, timeSpan, h, 0, 0);
        maxpos = max(position);
        idx = idx + 1;
        waitbar(idx/length(results),f, sprintf('Searching water touch param: %d %%', round(idx/length(results)*100)));
        if not(maxpos >= A && maxpos < B); continue; end
        maxacc = max(abs(CentralDifferentiation(velocity, h)/9.8506));
        if maxacc > 2; continue; end
        minimaIDX = islocalmin(position);
        pos10 = timeSpan(minimaIDX);
        jumpend = pos10(11);
        results(idx,1) = L;
        results(idx,2) = k;
        results(idx,3) = jumpend;
        results(idx,4) = maxacc;
    end
end
close(f);
opt =  results(results(:,3) > 0,:);
opt = opt(opt(:,3) == min(opt(:,3)), :);
len = opt(1);
spring = opt(2);
time = opt(3);
acc = opt(4);

end

