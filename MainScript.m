clear
format long g

% base parameters
h = 0.001;
timeSpan = 0:h:100;
maxBounces = 10;
defaultparams = true; 

switch defaultparams
    case true
        %%%%%% DEFAULT PARAMETERS DONT CHANGE %%%
        H = 74; % height of jump (m)
        DH = 31; % heights of deck from water (m)
        DECK =  H-DH; % distance from jump to deck
        c = 0.9; % drag coefficient (kg/m)
        m = 80; % jumper mass (kg)
        C = c/m; % drag / mass 
        L = 25; % bunge rope length (m)
        k = 90; % rope spring (N/m)
        K = k/m; % spring / mass
        g = 9.8; % gravity (m/s^2)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case false
        %%%%%% CUSTOMER PARAMETERS CHANGE AWAY %%%
        %%%%%% Can be used for playing/testing %%%
        H = 74; % height of jump (m)
        DH = 31; % heights of deck from water (m)
        DECK =  H-DH; % distance from jump to deck
        c = 0.9; % drag coefficient (kg/m)
        m = 80; % jumper mass (kg)
        C = c/m; % drag / mass 
        L = 25; % bunge rope length (m)
        k = 90; % rope spring (N/m)
        K = k/m; % spring / mass
        g = 9.8; % gravity (m/s^2)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    otherwise
        error('invalid parameter flag')
end

% velocity ode, captures parameters set above when run
% make sure to run the whole script together to avoid
% unrecognised changes
dvdt = @(y, v) g - C .* abs(v) .* v - max(0, K .*(y-L));


%%%%%%%%%%%%%%%%%%%%%%%%% TASK 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate position and velocity versus time
% ode evaluation using RK4 method
[position, velocity] = RK4Coupled(dvdt, timeSpan, h, 0, 0);

% plot position vs time with bounce annotations
figure('Position',[100 100 1500 500])
minimaIDX=islocalmin(position);
plot(timeSpan(minimaIDX),position(minimaIDX),'g*')
hold on
% this bit adds stars and numbers to bounces
iter = find(minimaIDX==1, maxBounces);
stopVal=0;
for ii = 1:length(find(minimaIDX==1, 10))
    item=iter(ii);
    text(timeSpan(item),position(item)-5,num2str(ii),'Color','k')
    if ii==length(find(minimaIDX==1, 10))
        stopVal=timeSpan(item);
    end
end
hold on
plot(timeSpan, position)
yline(H, 'b','River')
yline(0, 'k', 'Jump Point')
yline(DECK, 'k', 'Deck')
xline(stopVal, 'k', ['Stop (s): ',num2str(stopVal)], ...
    'LabelOrientation', 'aligned', 'LabelHorizontalAlignment', 'left')
title('Position (m from platform) vs Time')
subtitle('*Y Axis Reversed*')
xlabel('Time(s)')
ylabel('Relative Postion from Platform (m)')
ylim([0 80])
xlim([0 stopVal+5])
set(gca, 'YDir','reverse')
%%%%%%%%%%%%%%%%%%%%%%%%%%%% END 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%% TASK 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%to do: plot velocity and identify max speed


%%%%%%%%%%%%%%%%%%%%%%%%% END 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%% TASK 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%to do: numerically differentiate velocity to get acceleration
%suggest central difference method as the easiest with forward
%differece for the first point and backward diff for last

%acceleration = NumericalDiff(velocity, h)
%%%%%%%%%%%%%%%%%%%%%%%%% END 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%% TASK 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%to do: Numerically integrate to get distance travelled, 
%simpsons or trapezoid are the options we have it seems.

% distTravelled = DefIntegral(timeSpan, abs(velocity), h, 0, 60)
% sprintf('Total distance travelled: \n%.2f', DECK)
%%%%%%%%%%%%%%%%%%%%%%%%% END 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%% TASK 5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%to do: Implement function for divided diff interpolation
% div difference is easier to get the polynomial back as a 
% anon function for plotting

%subset timspan and position by the first bounce

%Interpolating function for 4 points either side of deck crossing

%iterate root finding method for exact crossing value to maybe 1E-5 error

%%%%%%%%%%%%%%%%%%%%%%%%% END 5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%% TASK 6 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure it out later.
%%%%%%%%%%%%%%%%%%%%%%%%% END 6 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
