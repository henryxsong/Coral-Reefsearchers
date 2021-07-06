syms C P T M mu1 mu2 q alpha sigma r h gamma beta g a t
% g = g(P), a = a(t)

%grazing intensity 'g'
%g = @(P) (alpha*P)/beta;

%sin function of 
%a = @(t) abs((0.9*(9*sin(pi*t)+1))/(10));

%solving disease free equilibrium
dMdt = a*C*M + gamma*T*M - (g*M)/(M+T);
M0 = solve(dMdt == 0, M); % solves for M when dMdt = 0
%M0 = subs(M0, [C,P,T], [0,0,0]);

dPdt = q*P*(1-(P/(beta*C))) - P*(h*mu2);
P0 = solve(dPdt == 0, P);

dTdt = mu1*C + (g*M0(1))/(M+T) - T*(r*C + gamma*M0(1));
T0 = solve(dTdt == 0, T);

dCdt = r*T*C + sigma*P*C - (a*M0(1) + mu1)*C;
C0 = solve(dCdt == 0, C);

% scriptF
% F: Jacobian matrix of script F
% scriptVi
% V: Jacobial matrix of script VI
% the largest absolute eigenvalue of F*inverse(V) is the basic reproduction
% number

%script F
% sF = [beta * S * I, 0];
% F = jacobian(sF, [S I]); % jacobian matrix
% F = subs(F, S, S0); %S0 = Pi/mu
% F = subs(F, I, 0);

%script V
% sV = [(sigma+mu)*E, (gamma + mu)*I-sigma*E];
% V = jacobian(sV, [E I]);

% F*V^-1
% FVInv = F * inv(V);

%eigenvalues of F*V^-1
% eigen = eig(FVInv);

% basic reproduction number
% R0 = eigen(2);
