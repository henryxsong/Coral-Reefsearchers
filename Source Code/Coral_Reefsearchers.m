% Team Coral Reefsearcher's Matlab Scripts


%--------------------------------------------------------------------------
%% Compartment - Single View - System of Differential Equations
clear; % Clears workspace
clc; % Clears Command Window

%---------------------------------------------
% Parameter Values
mu1 = 0.15; % mortality rate of coral reefs
mu2 = 0.22; % natural death rate of parrotfish
q = 0.47; % intrinsic growth rate for parrotfish
omega = 1; % %maximum grazing intensity
sigma = 0.01; % rate that parrotfish bite corals
r = 0.5; % rate that coral recruit to overgrow algal turfs
phi = 0.8; %rate that macroalgae spread vegetative over algal turfs
beta = 1; % carrying capacity

a0 = 0.99; % rate that coral is overgrown by macroalgae
h = 0.1; %<----CONTROL VARIABLE FOR GAME THEORY

%grazing intensity 'g'
g = @(P) (omega*P)/beta;

%sin function of 
a = @(t) abs((a0*(9*sin(pi*t)+1))/(10));
%---------------------------------------------

%---------------------------------------------
% Compartment Initial Conditions
C = 1/4;
P = 3/4;
T = 1/4;
M = 1/2;
Prop_Total = C + T + M;
IC = [C, P, T, M, Prop_Total];
%---------------------------------------------

% dC/dt = rTC + sigmaPC - C(aM + d)
% dP/dt = qP(1-P/betaC) - (h + mu)P 
% dT/dt = dC + (g(P)M)/(M + T) - (rC + phiM)T
% dM/dt = aMC + phiMT - (g(P)M)/(M + T)
%
% C = y(1), P = y(2), T = y(3), M = y(4), C+T+M = y(5)
f = @(t,y) [r*y(3)*y(1) + sigma*y(2)*y(1) - y(1)*(a(t)*y(4) + mu1),
        q*y(2)*(1-(y(2)/(beta*y(1)))) - (h+mu2)*y(2), 
        mu1*y(1) +  (g(y(2))*y(4))/(y(4)+y(3)) - (r*y(1) + phi*y(4))*y(3),
        a(t) * y(4)*y(1) + phi*y(4)*y(3) - (g(y(2))*y(4))/(y(4)+y(3)),
        y(1)+y(3)+y(4)];
            
[t,ya] = ode45(f, [0 5], IC);

figure
hold on
plot(t, ya(:,1), '+-.', 'Color', '#FFC996', 'Linewidth', 2.5)
plot(t, ya(:,2), 'x-.', 'Color', '#4974A5', 'Linewidth', 2.5)
plot(t, ya(:,3), 'o-.', 'Color', '#BDD2B6', 'Linewidth', 2.5)
plot(t, ya(:,4), '*-.', 'Color', '#CF0000', 'Linewidth', 2.5)
%legend('Coral (C)', 'Algal Turf (T)', 'Macroalgae (M)')
set(gca, 'FontSize',18);
ylim([0 1]);
legend('Coral (C)', 'Parrotfish (P)', 'Algal Turf (T)', 'Macroalgae (M)')
%text(0.25,0.05,txt, 'FontSize', 18);
xlabel('Time (Year)')
ylabel('Proportion of Population')
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
%% Compartment - Animation View - System of Differential Equations
clear; % Clears workspace
clc; % Clears Command Window

% Parameter Values
mu1 = 0.15; % mortality rate of coral reefs
mu2 = 0.22; % natural death rate of parrotfish
q = 0.47; % intrinsic growth rate for parrotfish
omega = 1; % %maximum grazing intensity
sigma = 0.01; % rate that parrotfish bite corals
r = 0.5; % rate that coral recruit to overgrow algal turfs
phi = 0.8; %rate that macroalgae spread vegetative over algal turfs
beta = 1; % carrying capacity

a0 = 0.99; % rate that coral is overgrown by macroalgae
h = 0.1; %<----CONTROL VARIABLE FOR GAME THEORY

%grazing intensity 'g'
g = @(P) (omega*P)/beta;

%sin function of 
a = @(t) abs((a0*(9*sin(pi*t)+1))/(10));
%---------------------------------------------

%---------------------------------------------
% Compartment Initial Conditions
C = 1/4;
P = 3/4;
T = 1/4;
M = 1/2;
Prop_Total = C + T + M;
IC = [C, P, T, M, Prop_Total]; %Initial Conditions array
%---------------------------------------------

% System of Differential Equations
    f = @(t,y) [r*y(3)*y(1) + sigma*y(2)*y(1) - y(1)*(a(t)*y(4) + mu1),
       q*y(2)*(1-(y(2)/(beta*y(1)))) - (h+mu2)*y(2), 
       mu1*y(1) +  (g(y(2))*y(4))/(y(4)+y(3)) - (r*y(1) + phi*y(4))*y(3),
       a(t) * y(4)*y(1) + phi*y(4)*y(3) - (g(y(2))*y(4))/(y(4)+y(3)),
       y(1)+y(3)+y(4)];
    
    % Solve using ODE45
    [t,ya] = ode45(f, [0 5], IC);
    const_t = t; % constant used since t changes each time ode45 is calculated
    
for i = 1:length(const_t)
    h = i/length(const_t); %variable to animate

   % System of Differential Equations
    f = @(t,y) [r*y(3)*y(1) + sigma*y(2)*y(1) - y(1)*(a(t)*y(4) + mu1),
       q*y(2)*(1-(y(2)/(beta*y(1)))) - (h+mu2)*y(2), 
       mu1*y(1) +  (g(y(2))*y(4))/(y(4)+y(3)) - (r*y(1) + phi*y(4))*y(3),
       a(t) * y(4)*y(1) + phi*y(4)*y(3) - (g(y(2))*y(4))/(y(4)+y(3)),
       y(1)+y(3)+y(4)];
    
    % Solve using ODE45
    [t,ya] = ode45(f, [0 5], IC);
    
    % Plot
    txt = ['h = ' num2str(h)]; % shows value of param value at iteration
    
    fig = figure;
    hold on
    plot(t, ya(:,1), '+-.', 'Color', '#FFC996', 'Linewidth', 2.5)
    plot(t, ya(:,2), 'x-.', 'Color', '#4974A5', 'Linewidth', 2.5)
    plot(t, ya(:,3), 'o-.', 'Color', '#BDD2B6', 'Linewidth', 2.5)
    plot(t, ya(:,4), '*-.', 'Color', '#CF0000', 'Linewidth', 2.5)
    
    set(gca, 'FontSize',18); % sets axis & legend font size to 18
    ylim([0 1]); % sets y-axis limit to always be 0-1
    legend('Coral (C)','Parrotfish (P)','Algal Turf (T)','Macroalgae (M)')
    text(0.25,0.05,txt, 'FontSize', 18); % displays text on plot
    
    xlabel('Time (Year)')
    ylabel('Proportion of Population')
    
    % automatically save figure into root directory (where this .m file is
    % stored)
    fname = append('Frame-', num2str(i)); %file name of current iteration
    %saveas(fig, fname, 'png'); %save figure as .png
end
%---------------------------------------------------------------
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
%% Disease Free Equilibrium
clear; % Clears workspace
clc; % Clears Command Window

%---------------------------------------------
% Symbolic Definitions
syms C P T M mu1 mu2 q omega sigma beta r h phi g a t
% g = g(P), a = a(t)

g = @(P) (omega*P)/beta; %grazing intensity 'g'

%a = @(t) abs((0.9*(9*sin(pi*t)+1))/(10)); %sin function of a(t)
%---------------------------------------------

M_0 = 0;

dTdt = mu1*C + (g(P)*M_0)/(M_0+T) - T*(r*C + phi*M_0);
T_0 = solve(dTdt == 0, T);

C_0 = 1 - T_0;

dPdt = q*P*(1-(P/(beta*C))) - P*(h+mu2);
P_0 = solve(dPdt == 0, P);
P_0 = subs(P_0, C, C_0);

%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
%% R0
%Note: RUN 'DISEASE FREE EQUILIBRIUM' SECTION FIRST

%---------------------------------------------
% R0
%script F
sF = [a*C*M + phi*T*M];
F = jacobian(sF, [M]); % jacobian matrix
F = subs(F, T, T_0); 

%script V
sV = [(g(P)*M)/(M+T)];
V = jacobian(sV, [M]);
V = subs(V, M, M_0);

%eigenvalues of F*V^-1
eigens = eig(F * inv(V));

% basic reproduction number
R0 = eigens(1);
R0 = subs(R0, T, T_0);
R0 = subs(R0, C, C_0);
R0 = subs(R0, P, P_0(2));
%---------------------------------------------

%---------------------------------------------
% Plot R0

%---------------------------------------------
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
%% Herd Immunity
% Note: MUST RUN 'DISEASE FREE EQUILIBRIUM' & 'R0' SECTION FIRST
param_array = [mu1, mu2, q, omega, sigma, r, beta, a, phi]; %when a0 =
% 0.99
param_values = [0.15, 0.22, 0.9, 1, 0.01, 0.5, 1, 0.5, 0.8]; %a(t) = 0.5

herd_immunity_equation = R0;
herd_immunity_value = solve(R0 == 1, h);
for i = 1:length(param_array)
    herd_immunity_equation = subs(herd_immunity_equation, param_array(i), param_values(i));
    herd_immunity_value = subs(herd_immunity_value, param_array(i), param_values(i));
end

herd_immunity_point = [double(herd_immunity_value), 1];

h_var = 0:0.001:1;

figure
hold on
fplot(herd_immunity_equation, [0 1])
plot(herd_immunity_point(1), herd_immunity_point(2), 'o')
plot([herd_immunity_point(1), herd_immunity_point(1)], [0, herd_immunity_point(2)], 'r-')
plot([0, herd_immunity_point(1)], [herd_immunity_point(2), herd_immunity_point(2)], 'r-')
xlim([0 1]);
ylim([0 2]);

set(gca, 'FontSize',18);
title('Herd Immunity');
xlabel('h')
ylabel('R_0')
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
%% Sensitivity Analysis
% Note: MUST RUN 'R0' SECTION FIRST

param_array = [mu1, mu2, q, omega, sigma, r, phi, beta, h, a];
param_values = [0.15, 0.22, 0.47, 1, 0.01, 0.5, 0.8, 1, 0.1, 0.5];
sens_analysis = [];

for i = 1:length(param_array)
    R0_diff = diff(R0, param_array(i));
    %sens_analysis = cat(1, sens_analysis, diff(R0, param_array(i)));
    for j = 1:length(param_array)
        R0_diff = subs(R0_diff, param_array(j), param_values(j));
    end
    sens_analysis = cat(1, sens_analysis, double(R0_diff));
    %subs_ans = cat(1, subs_ans, [param_array(i) double(R0_diff)]);
end 
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
%% Endemic Equilibrium
clear; % Clears workspace
clc; % Clears Command Window

%---------------------------------------------
% Symbolic Definitions
syms C P T M mu1 mu2 q omega sigma beta r h phi g a t
% g = g(P), a = a(t)

g = @(P) (omega*P)/beta; %grazing intensity 'g'

a = @(t) abs((0.9*(9*sin(pi*t)+1))/(10)); %sin function of a(t)
%---------------------------------------------

%---------------------------------------------
% Equation
T_E = (mu1 + a*M)/r;

C_E = 1 - (T+M);
C_E = subs(C_E, T, T_E)

P_E = (beta*C*(q-(h+mu2)))/q;
P_E = subs(P_E, C, C_E)

M_E = (omega*P)/(beta*(a*C+phi*T))-T - M;
M_E = subs(M_E, P, P_E);
M_E = subs(M_E, C, C_E);
M_E = subs(M_E, T, T_E);
M_E_equation = solve(M_E == 0, M)
%---------------------------------------------

%---------------------------------------------
% Value
param_array = [mu1, mu2, q, omega, sigma, r, phi, beta, h, t];
param_values = [0.15, 0.22, 0.47, 1, 0.01, 0.5, 0.8, 1, 0.1, 0.5];

M_E_Value = M_E_equation;

for i = 1:length(param_array)
    M_E_Value = subs(M_E_Value, param_array(i), param_values(i));
end

disp(M_E_Value)
%--------------------------------------------------------------------------


