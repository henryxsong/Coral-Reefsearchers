syms C P T M mu1 mu2 q alpha sigma r h gamma beta g a t
% g = g(P), a = a(t)

%grazing intensity 'g'
g = @(P) (alpha*P)/beta;

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
M_E = solve(dMdt == 0, M);

dPdt = q*P*(1-(P/(beta*C))) - P*(h*mu2);
P_E = solve(dPdt == 0, P);

%dTdt = mu1*C + (g*M)/(M+T) - T*(r*C + gamma*M);
dTdt = (mu1 + a*M_E(2))/(r);
T_E = solve(dTdt == 0, T);

%dCdt = r*T_E(2)*C + sigma*P_E(2)*C - (a*M_E(2) + mu1)*C - C;
dCdt = r*T*C + sigma*P*C - (a*M + mu1)*C - C;
C_E = solve(dCdt == 0, C, ReturnConditions = true);
%------------------------------------------------------


%------------------------------------------------------
% Solving R0
% scriptF
% F: Jacobian matrix of script F
% scriptVi
% V: Jacobial matrix of script VI

C_0 = q - h - mu2;
%script F
sF = [mu1*C + (g*M)/(M+T), 0];
F = jacobian(sF, [C T]); % jacobian matrix
F = subs(F, C, C_0); %S0 = Pi/mu
F = subs(F, T, 0);

%script V
sV = [r*C*T + gamma*M*T, -(a*C*M + gamma*T*M) + (g*M)/(M+T)];
V = jacobian(sV, [T M]);

%eigenvalues of F*V^-1
eigens = eig(F * inv(V));

% basic reproduction number
R0 = eigens(2);
%------------------------------------------------------
