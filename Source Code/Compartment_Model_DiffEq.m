

r = 10;
d = 2; %unknown
a = 0.1; %unknown
gamma = 0.8; %unknown
g = 10;
beta = 21; %unknown
kapa = 0.6; %known but not used
mu = 0.22;
h = 0.3; %<----CONTROL VARIABLE FOR GAME THEORY
alpha = 1;
q = 0.47;
sigma = 0.01; %unknown

g = @(P) (alpha*P)/beta;

% set up DFE
% dC/dt = rTC + sigmaPC - C(aM + d) #ignore sigmaPC for now
% dP/dt = qP(1-P/betaC) + kapaP - (h + mu)P 
% #remove kapaP, so = qP(1-P/betaC) - (h + mu)P
% dT/dt = dC + (g(P)M)/(M + T) - (rC + gammaM)T
% dM/dt = aMC + gammaMT - (g(P)M)/(M + T)
%
% C = y(1), P = y(2), T = y(3), M = y(4)
f = @(t,y) [r*y(3)*y(1) + sigma*y(2)*y(1) - y(1)*(a*y(4) + d),
            q*y(2)*(1-(y(2)/(beta*y(1)))) - (h+mu)*y(2), 
            d*y(1) + (g(y(2))*y(4))/(y(4)+y(3)) - (r*y(1) + gamma*y(4))*y(3),
            a * y(4)*y(1) + gamma*y(4)*y(3) - (g(y(2))*y(4))/(y(4)+y(3))];

C = 1/4;
P = 1;
T = 1/2;
M = 1/4;
        
  
% solve with ODE 45
[t,ya] = ode45(f, [0 5], [C, P, T, M]);

% graph
figure
hold on
plot(t, ya(:,1), '+-.', 'Color', '#FFC996', 'Linewidth', 2.5)
plot(t, ya(:,2), 'x-.', 'Color', '#A7D0CD', 'Linewidth', 2.5)
plot(t, ya(:,3), 'o-.', 'Color', '#BDD2B6', 'Linewidth', 2.5)
plot(t, ya(:,4), '*-.', 'Color', '#CF0000', 'Linewidth', 2.5)
legend('Coral (C)', 'Parrotfish (P)', 'Algal Turf (T)', 'Macroalgae (M)')
xlabel('Time (Year)')
ylabel('Proportion of Population')