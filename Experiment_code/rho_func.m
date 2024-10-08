function [rho,alpha,xi,DE,beta_clip1] = rho_func(delta,K,calT,T_total,d,hat_p,Lmax,N,X_diff,theta_hat,varepsilon)
% confidence radius rho
%% confidence parameters
deltav = delta/(15*K*T_total^3);
deltam = delta/(15*K*N*T_total^3);
deltad = delta/(15*N*T_total^3);
%% VE
alpha = d/T_total*log(2/deltav)+sqrt((d/T_total*log(2/deltav))^2+4*d/T_total*log(2/deltav));
%% DE
phi = min(4*max(hat_p,25/4*Lmax/T_total*log(2/deltad)),1/4);
tilbeta2 = 2.*phi*Lmax/T_total*log(2/deltad);
beta_clip1 = min(1/3*Lmax/T_total*log(2/deltad)+sqrt((1/3*Lmax/T_total*log(2/deltad))^2+tilbeta2),1);
optimal_var = -X_diff'*theta_hat*beta_clip1/sum(beta_clip1);
DE = abs(min(max(X_diff'*theta_hat,-2),2)+optimal_var)*beta_clip1;
%% RE
temp = 25/4*Lmax/T_total*log(2/deltad);
temp1 = d/4*log(2/deltam);
psi = hat_p > temp;
tilpsi = calT < temp1;
%case1
xi_1 = tilpsi'*4 * beta_clip1;
%case2
temp2 = hat_p < temp & calT > temp1;
lemmab1 = d./max(1,calT)*log(2/deltam)+sqrt((d./max(1,calT)*log(2/deltam)).^2+4*d./max(1,calT)*log(2/deltam)); % improved in a similar way as alpha
xi_2 = temp2'*(min(2*lemmab1,4).*beta_clip1);
%case3
xi_31 = sum(psi)*50*sqrt(2)*sqrt(d*Lmax)/T_total*log(2/deltam);
xi_32 = psi'*(min(2*lemmab1,4).*beta_clip1);
xi_3 = min(xi_31,xi_32);

xi = xi_1+xi_2+xi_3;
%% total
rho = 2*alpha+xi+DE;

end

