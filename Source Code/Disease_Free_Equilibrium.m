syms C P T M mu1 mu2 q omega sigma beta r h phi g a t
% g = g(P), a = a(t)

%grazing intensity 'g'
g = @(P) (omega*P)/beta;

%sin function of 
%a = @(t) abs((0.9*(9*sin(pi*t)+1))/(10));

%-----------------------------------------------------
%solving DFE
% dMdt = a*C*M + phi*T*M - (g*M)/(M+T);
% %M_0 = solve(dMdt == 0, M); % solves for M when dMdt = 0
% %M0 = subs(M0, [C,P,T], [0,0,0]);
% M_0 = 0;
% 
% dTdt = mu1*C + (g(P)*M_0)/(M_0+T) - T*(r*C + phi*M_0);
% T_0 = solve(dTdt == 0, T);
% 
% %dCdt = r*T_0*C + sigma*P_0*C - (a*M_0 + mu1)*C;
% %C_0 = solve(dCdt == 0, C);
% C_0 = 1 - T_0;
% 
% dPdt = q*P*(1-(P/(beta*C))) - P*(h+mu2);
% P_0 = solve(dPdt == 0, P);
% %P_0 = 0;
% %------------------------------------------------------
% 
% %------------------------------------------------------
% % % Solving R0
% 
% %script F
% sF = [a*C*M + phi*T*M];
% F = jacobian(sF, [M]); % jacobian matrix
% F = subs(F, T, T_0); 
% %F = subs(F, M, 0);
% %F = subs(F, C, C_0);
% 
% %script V
% syms G;
% sV = [(g(P)*M)/(M+T)];
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
% dMdt = a*C*M + phi*T*M - (g*M)/(M+T);
% M_E = solve(dMdt == 0, M);
% 
% dTdt = mu1*C + (g*M)/(M+T) - T*(r*C + phi*M);
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

%--------------------------------------------------------------------------
% Sensitivity Analysis
% R0 = -(beta*mu1*(a*(mu1/r - 1) - (mu1*phi)/r))/(P*omega*r);
% param_array = [mu1, mu2, q, omega, sigma, r, phi, beta, h];
% param_values = [0.15, 0.22, 0.47, 1, 0.01, 0.5, 0.8, 1, 0.1];
% 
% 
% for i = 1:length(param_array)
%     sensAns(i) = diff(R0, param_array(i));
%     %subsAns(i) = subs(sensAns[i], (param_array), (param_values));
% end 
%--------------------------------------------------------------------------

% T_E = (mu1 + a*M)/r;
% 
% C_E = 1 - (T+M);
% C_E = subs(C_E, T, T_E);
% 
% P_E = (beta*C*(q-(h+mu2)))/q;
% P_E = subs(P_E, C, C_E);
% %testP = (beta*(1-((mu1+a*M)/r) + M)*(q-(h+mu2)))/(q);
% 
% M_E = (omega*P)/(beta*(a*C+phi*T))-T;
% M_E = subs(M_E, P, P_E);
% M_E = subs(M_E, C, C_E);
% M_E = subs(M_E, T, T_E);
% %testM = (omega*(beta*(1-(mu1+a*M)/r + M)*(q-(h+mu2))))/(q*beta*(a*(1-((mu1+a*M)/r + M))) + phi*((mu1+a*M)/r)) - (mu1+a*M)/r;
% ansM = solve(M_E == 0, M);
% l = latex(ansM);
% 
% % My quadratic equation
% a_quad = a^2*q*phi - a^2*q*r - a^3*q;
% b_quad = a^2*q*r - a*mu1*q*r - 2*a^2*mu1*q + 2*a*mu1*q*phi + omega*a*q - omega*a*h - omega*a*mu2 + omega*r*q - omega*r*h - omega*r*mu2;
% c_quad = a*mu1*q*r - a*mu1^2*q + mu1^2*q*phi - omega*r*q + omega*r*h + omega*r*mu2 + omega*mu1*q - omega*mu1*h - omega*mu1*mu2;
% 
% % Dr. Choi's quadratic equation
% % a_quad = a*q*(a+r-phi)*(a+r)^2;
% % b_quad = (a+r)*(2*a*q*(a+r-phi)*(r-mu1) - r*(omega*(a+r)*(q-h-mu2) - a*q*(a+r) + g*mu1*phi));
% % c_quad = a*q*(a+r-phi)*(r-mu1)^2 + r*(omega*(a+r)*(g-h-mu2) - a*q*(a+r) + q*mu1*phi)*(r-mu1) - q*phi*r^2*(mu1+a);
% % quad_quad = [((-b_quad)+((b_quad^2-4*a_quad*c_quad))^0.5)/(2*a_quad); ((-b_quad)-((b_quad^2-4*a_quad*c_quad))^0.5)/(2*a_quad)];
% 
% quad_ans = solve(M^2*a_quad + M*b_quad + c_quad == 0, M);