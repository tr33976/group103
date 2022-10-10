tic
clc
clear 

jumpHeight = 74;
persHeight = 1.75;
A=jumpHeight-persHeight;
B=A+0.1;

h=0.001;
timeSpan = 0:h:100;

%%%%%% DEFAULT PARAMETERS minus sweep vals %%%
c = 0.9; % drag coefficient (kg/m)
m = 80; % jumper mass (kg)
C = c/m; % drag / mass 
g = 9.8; % gravity (m/s^2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u = @(a,b) (b-a)* rand() + a;

rope = [10, 60];
spring =  [40, 200];

for ii = 1

q1_r_u = u(rope(1), median(rope));
q1_s_u = u(spring(1), median(spring));

q2_r_u = u(median(rope)+1E-8, rope(end));
q2_s_u = u(spring(1), median(spring));

q3_r_u = u(rope(1), median(rope));
q3_s_u = u(median(spring)+1E-8, spring(end));

q4_r_u = u(median(rope)+1E-8, rope(end));
q4_s_u = u(median(spring)+1E-8, spring(end));

maxpos = inf(4,1);
maxacc = inf(4,1);
jumpend = inf(4,1);
kval = inf(4,1);
lval = inf(4,1);

lengths = [q1_r_u,q2_r_u,q3_r_u,q4_r_u];
springC = [q1_s_u,q2_s_u,q3_s_u,q4_s_u];


counter = 1;
for i = 1:4
    cont = true;
    abancount = 0;
    while cont
    L = lengths(i);
    k = springC(i);
    dvdt = @(y, v) g - C .* abs(v) .* v - max(0, k/m .*(y-L));
    [position, velocity] = RK4Coupled(dvdt, timeSpan, h, 0, 0, false);
    minimaIDX = islocalmin(position);
    pos10 = timeSpan(minimaIDX);
    kval(i) = k;
    lval(i) = L;
    jumpend(i) = pos10(11);
    abancount = abancount + 1;
    if abancount == 50; break; end
    if maxacc(counter) > 2; continue; end
    if not(maxpos(counter) >= A-5 & maxpos(counter) < B+5); continue; end
    maxpos(counter) = max(position);
    maxacc(counter) = max(abs(CentralDifferentiation(velocity, h)/9.8506));
    cont = false;
    end
    counter = counter + 1;
end

idx = [1;2;3;4];


% maxidx = not(maxacc == max(maxacc));
% 
% 
% accidx =  find(maxacc < 2);
% maxacc = maxacc(accidx,:);
% maxpos = maxpos(accidx,:);
% jumpend = jumpend(accidx,:);
% idx = idx(accidx, :);
% kval = kval(accidx, :);
% lval = lval(accidx, :);
% 
% [closest, closesidx] = min(abs(maxpos-A));
% maxpos = maxpos(closesidx, :);
% maxacc = maxacc(closesidx,:);
% jumpend = jumpend(closesidx,:);
% idx = idx(closesidx, :);
% kval = kval(closesidx, :);
% lval = lval(closesidx, :);

if(lval > median(rope)); rope = [median(rope) rope(end)];
else; rope = [rope(1) median(rope)];end

if(kval > median(spring)); spring = [median(spring) spring(end)];
else; spring = [spring(1) median(spring)]; end

end

toc



