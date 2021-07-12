syms C P T M mu1 mu2 q alpha sigma r h gamma beta g a t
% g = g(P), a = a(t)

%grazing intensity 'g'
g = @(P) (alpha*P)/beta;

%sin function of 
%a = @(t) abs((0.9*(9*sin(pi*t)+1))/(10));

%-----------------------------------------------------
%solving DFE
% dMdt = a*C*M + gamma*T*M - (g*M)/(M+T);
% %M_0 = solve(dMdt == 0, M); % solves for M when dMdt = 0
% %M0 = subs(M0, [C,P,T], [0,0,0]);
% M_0 = 0;
% 
% dTdt = mu1*C + (g(P)*M_0)/(M_0+T) - T*(r*C + gamma*M_0);
% T_0 = solve(dTdt == 0, T);
% 
% %dCdt = r*T_0*C + sigma*P_0*C - (a*M_0 + mu1)*C;
% %C_0 = solve(dCdt == 0, C);
% C_0 = 1 - T_0;
% 
% dPdt = q*P*(1-(P/(beta*C))) - P*(h+mu2);
% P_0 = solve(dPdt == 0, P);
% %P_0 = 0;
%------------------------------------------------------

% %------------------------------------------------------
% % % Solving R0
% 
% %script F
% sF = [a*C*M + gamma*T*M];
% F = jacobian(sF, [M]); % jacobian matrix
% F = subs(F, T, T_0); 
% %F = subs(F, M, 0);
% %F = subs(F, C, C_0);
% 
% %script V
% syms G;
% sV = [(G*M)/(M+T)];
% V = jacobian(sV, [M]);
% V = subs(V, M, M_0);
% 
% %eigenvalues of F*V^-1
% eigens = eig(F * inv(V));
% 
% % basic reproduction number
% R0 = eigens(1);
% R0 = subs(R0, C, C_0);
% R0 = subs(R0, T, T_0);
%------------------------------------------------------


%------------------------------------------------------
%solving Endemic
% dPdt = q*P*(1-(P/(beta*C))) - P*(h*mu2);
% P_E = solve(dPdt == 0, P);
% 
% dMdt = a*C*M + gamma*T*M - (g*M)/(M+T);
% M_E = solve(dMdt == 0, M);
% 
% dTdt = mu1*C + (g*M)/(M+T) - T*(r*C + gamma*M);
% T_E = solve(dTdt == 0, T);
% 
% P_E = P_E(2);
% M_E = M_E(2);
% 
% M_E = subs(M_E, P, P_E);
% 
% % T_E = subs(T_E, M, M_E);
% % T_E = subs(T_E(1), T, T_E);
% % T_E = subs(T_E(2), T, T_E);
% 
% %dCdt = r*T_E(2)*C + sigma*P_E(2)*C - (a*M_E(2) + mu1)*C - C;
% dCdt = r*T*C + sigma*P*C - (a*M + mu1)*C - C;
% C_E = solve(dCdt == 0, C, ReturnConditions = true);
% %------------------------------------------------------

T_E = (mu1 + a*M)/r;

C_E = 1 - (T+M);
C_E = subs(C_E, T, T_E);

P_E = (beta*C*(q-(h+mu2)))/q;
P_E = subs(P_E, C, C_E);
%testP = (beta*(1-((mu1+a*M)/r) + M)*(q-(h+mu2)))/(q);

M_E = (alpha*P)/(beta*(a*C+gamma*T))-T;
M_E = subs(M_E, P, P_E);
M_E = subs(M_E, C, C_E);
M_E = subs(M_E, T, T_E);
%testM = (alpha*(beta*(1-(mu1+a*M)/r + M)*(q-(h+mu2))))/(q*beta*(a*(1-((mu1+a*M)/r + M))) + gamma*((mu1+a*M)/r)) - (mu1+a*M)/r;
ansM = solve(M_E == 0, M);
l = latex(simplify(ansM(2)));