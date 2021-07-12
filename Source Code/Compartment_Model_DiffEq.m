%%% Compartment Model Calculations
%%% Team Coral Reefsearchers
%
% Note: In order to use Animation View, function 'f' and variables [t, ya]
%       must be loaded into the workspace first by running Single View then
%       commenting out Single View.
%--------------------------------------------------------------------------
%% Parameter Values
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

%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
%% Compartment Initial Conditions
C = 1/4;
P = 3/4;
T = 1/4;
M = 1/2;
Prop_Total = C + T + M;
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
%% Single View - System of Differential Equations
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
            
[t,ya] = ode45(f, [0 5], [C, P, T, M, Prop_Total]);

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
%% Animation View - System of Differential Equations
% for i = 1:length(t)
%     a0 = i/length(t);
%     a = @(t) abs((a0*(9*sin(pi*t)+1))/(10));
%     
%     % System of Differential Equations
%     f = @(t,y) [r*y(3)*y(1) + sigma*y(2)*y(1) - y(1)*(a(t)*y(4) + mu1),
%             q*y(2)*(1-(y(2)/(beta*y(1)))) - (h+mu2)*y(2), 
%             mu1*y(1) + (g(y(2))*y(4))/(y(4)+y(3)) - (r*y(1) + phi*y(4))*y(3),
%             a(t) * y(4)*y(1) + phi*y(4)*y(3) - (g(y(2))*y(4))/(y(4)+y(3)),
%             y(1)+y(3)+y(4)];
%     
%     % Solve using ODE45
%     [t,ya] = ode45(f, [0 5], [C, P, T, M, 1]);
%     
%     % Plot
%     txt = ['a_0 = ' num2str(a0)]; % shows value of param value at iteration
%     
%     fig = figure;
%     hold on
%     plot(t, ya(:,1), '+-.', 'Color', '#FFC996', 'Linewidth', 2.5)
%     plot(t, ya(:,2), 'x-.', 'Color', '#4974A5', 'Linewidth', 2.5)
%     plot(t, ya(:,3), 'o-.', 'Color', '#BDD2B6', 'Linewidth', 2.5)
%     plot(t, ya(:,4), '*-.', 'Color', '#CF0000', 'Linewidth', 2.5)
%     
%     set(gca, 'FontSize',18); % sets axis & legend font size to 18
%     ylim([0 1]); % sets y-axis limit to always be 0-1
%     legend('Coral (C)','Parrotfish (P)','Algal Turf (T)','Macroalgae (M)')
%     text(0.25,0.05,txt, 'FontSize', 18); % displays text on plot
%     
%     xlabel('Time (Year)')
%     ylabel('Proportion of Population')
%     
%     % automatically save figure into root directory (where this .m file is
%     % stored)
%     fname = append('Frame-', num2str(i)); %file name of current iteration
%     saveas(fig, fname, 'png'); %save figure as .png
% end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
%% Sensitivity Analysis
R0 = -(beta*mu1*(a*(mu1/r - 1) - (mu1*phi)/r))/(P*omega*r);
param_array = [mu1, mu2, q, omega, sigma, r, phi, beta, a0];

for i = 1:length(param_array)
    p_

end
%--------------------------------------------------------------------------