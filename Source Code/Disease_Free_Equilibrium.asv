syms C P T M mu1 mu2 q alpha sigma r h gamma beta g a t
% g = g(P), a = a(t)

%grazing intensity 'g'
%g = @(P) (alpha*P)/beta;

%sin function of 
%a = @(t) abs((0.9*(9*sin(pi*t)+1))/(10));

%-----------------------------------------------------
%solving DFE
dMdt = a*C*M + gamma*T*M - (g*M)/(M+T);
M_0 = solve(dMdt == 0, M); % solves for M when dMdt = 0
%M0 = subs(M0, [C,P,T], [0,0,0]);

dPdt = q*P*(1-(P/(beta*C))) - P*(h*mu2);
P_0 = solve(dPdt == 0, P);

dTdt = mu1*C + (g*M_0(1))/(M_0(1)+T) - T*(r*C + gamma*M_0(1));
T_0 = solve(dTdt == 0, T);

dCdt = r*T_0*C + sigma*P_0(2)*C - (a*M_0(1) + mu1)*C;
C_0 = solve(dCdt == 0, C);
%------------------------------------------------------


%------------------------------------------------------
%solving Endemic
dMdt = a*C*M + gamma*T*M - (g*M)/(M+T);
dPdt = q*P*(1-(P/(beta*C))) - P*(h*mu2);
dTdt = mu1*C + (g*M)/(M+T) - T*(r*C + gamma*M);
dCdt = r*T*C + sigma*P*C - (a*M + mu1)*C;

M_E = solve(dMdt == 0, M);
P_E = solve(dPdt == 0, P);
T_E = solve(dTdt == 0, T);
C_E = solve(dCdt == 0, C, ReturnConditions = true);
%------------------------------------------------------


% scriptF
% F: Jacobian matrix of script F
% scriptVi
% V: Jacobial matrix of script VI
% the largest absolute eigenvalue of F*inverse(V) is the basic reproduction
% number

%script F
sF = [beta * S * I, 0];
F = jacobian(sF, [S I]); % jacobian matrix
F = subs(F, S, S0); %S0 = Pi/mu
F = subs(F, I, 0);

%script V
% sV = [(sigma+mu)*E, (gamma + mu)*I-sigma*E];
% V = jacobian(sV, [E I]);

%eigenvalues of F*V^-1
% eigens = eig(F * inv(V));

% basic reproduction number
% R0 = eigen(2);
