% parameter values
mu1 = 0.15; % mortality rate of coral reefs
mu2 = 0.22; % natural death rate of parrotfish
q = 0.47; % intrinsic growth rate for parrotfish
alpha = 1; % %maximum grazing intensity
sigma = 0.01; % rate that parrotfish bite corals
r = 0.5; % rate that coral recruit to overgrow algal turfs
gamma = 0.8; %rate that macroalgae spread vegetative over algal turfs
beta = 1; % carrying capacity

a0 = 0.99; % rate that coral is overgrown by macroalgae
h = 0.1; %<----CONTROL VARIABLE FOR GAME THEORY

%grazing intensity 'g'
g = @(P) (alpha*P)/beta;

%sin function of 
a = @(t) abs((a0*(9*sin(pi*t)+1))/(10));

% set up SoDEx
% dC/dt = rTC + sigmaPC - C(aM + mu1)
% dP/dt = qP(1-P/betaC) + kapaP - (h + mu2)P #remove kapaP, so  = qP(1-P/betaC) - (h + mu)P
% dT/dt = mu1C + (g(P)M)/(M + T) - (rC + gammaM)T
% dM/dt = aMC + gammaMT - (g(P)M)/(M + T)
%
% C = y(1), P = y(2), T = y(3), M = y(4)
f = @(t,y) [r*y(3)*y(1) + sigma*y(2)*y(1) - y(1)*(a(t)*y(4) + mu1),
            q*y(2)*(1-(y(2)/(beta*y(1)))) - (h+mu2)*y(2), 
            mu1*y(1) + (g(y(2))*y(4))/(y(4)+y(3)) - (r*y(1) + gamma*y(4))*y(3),
            a(t) * y(4)*y(1) + gamma*y(4)*y(3) - (g(y(2))*y(4))/(y(4)+y(3)),
            y(1)+y(3)+y(4)];

C = 1/3;
P = 0.75;
T = 1/3;
M = 0;
        
  
% solve with ODE 45
[t,ya] = ode45(f, [0 5], [C, P, T, M, 1]);

% graph
figure
hold on
plot(t, ya(:,1), '+-.', 'Color', '#FFC996', 'Linewidth', 2.5)
plot(t, ya(:,2), 'x-.', 'Color', '#4974A5', 'Linewidth', 2.5)
plot(t, ya(:,3), 'o-.', 'Color', '#BDD2B6', 'Linewidth', 2.5)
plot(t, ya(:,4), '*-.', 'Color', '#CF0000', 'Linewidth', 2.5)
%legend('Coral (C)', 'Algal Turf (T)', 'Macroalgae (M)')
legend('Coral (C)', 'Parrotfish (P)', 'Algal Turf (T)', 'Macroalgae (M)')
xlabel('Time (Year)')
ylabel('Proportion of Population')