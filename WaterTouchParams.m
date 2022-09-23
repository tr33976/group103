function [len, spring, time] = WaterTouchParams(timeSpan, h, persHeight, jumpheight)

A=jumpheight-persHeight;
B=A+0.1;

ropeLengths = 40:0.5:60;
springCoefs = 60:0.5:85;

results = zeros((length(ropeLengths)*length(ropeLengths)), 3);

%%%%%% DEFAULT PARAMETERS minus sweep vals %%%
c = 0.9; % drag coefficient (kg/m)
m = 80; % jumper mass (kg)
C = c/m; % drag / mass 
g = 9.8; % gravity (m/s^2)
L=25;
k=90;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dvdt = @(y, v) g - C .* abs(v) .* v - max(0, k/m .*(y-L));
totalp = length(results);

idx = 1;
counter = 1;
for L = ropeLengths
    for k = springCoefs
        % velocity ode
        [position, velocity] = RK4Coupled(dvdt, timeSpan, h, 0, 0);
        maxpos = max(position);
        if mod(totalp, counter) == totalp*0.1;disp([num2str(counter/totalp*10), '% done with param sweep']);end
        counter = counter +1;
        if not(maxpos >= A && maxpos < B); continue; end
        maxacc = max(abs(CentralDifferentiation(velocity, h)/9.8506));
        if maxacc > 2; continue; end
        minimaIDX=islocalmin(position);
        pos10 = timeSpan(minimaIDX);
        jumpend = pos10(11);
        results(idx,1) = L;
        results(idx,2) = k;
        results(idx,3) = jumpend;
        idx = idx + 1;
    end
end

opt =  results(results(:,3) > 0,:);
opt = opt(opt(:,3) == min(opt(:,3)), :);
len = opt(1);
spring = opt(2);
time = opt(3);

end

