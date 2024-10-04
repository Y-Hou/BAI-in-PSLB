function [tau_PSBAI,hatx_varepsilon,flag_PSBAI,l_PSBAI,tau_naive,x_dot_varepsilon,flag_naive,l_naive,rho,alpha,xi,DE,beta_clip1] = PSBAI_plus_func(X,Theta,ptheta,Ptheta,Lmin_gt,Lmax,varepsilon,w,gamma,b,delta,C_fixed,seed)
%The PS \epsilon BAI^+ algorithm, intergaretd with N \epsilon BAI
%%
rng(seed)
tic
%default output
tau_PSBAI = -1;
hatx_varepsilon = -1;
l_PSBAI = -1;
l_naive = -1;

rho = [-1;-1];
alpha = -1;
xi = -1;
DE =[-1;-1];
beta_clip1 =[-1;-1];

tau_naive=-1;
x_dot_varepsilon = -1;

flag_PSBAI = 0;
flag_naive = 0;

%input
[d,K] = size(X);
X_index = 1:K;
N = length(Ptheta);
tau_star = 38400*log(80)*N*Lmax/min(varepsilon,1)^2*log(N^2*K*Lmax/(min(varepsilon,1)^2*delta));

%initialization
lambda = minvol(X);     %Compute the G-optimal allocation
Alambda = zeros(d);
Lambda = zeros(K,1);    %cdf of lambda
for i=1:K
    Alambda = Alambda + lambda(i).*(X(:,i)*X(:,i)');
    Lambda(i) = sum(lambda(1:i));
end
Alambda_inv = inv(Alambda);

CDsample = zeros(d,w);
n_cd = 0;
CAid = {};
tCD = inf;
C3 = 1.2021;
tstar = 3*d*log(6*d*K*C3/delta);
%statistics for easier computing
cum_theta = zeros(d,N);
theta_hat = zeros(d,N);
calT = ones(N,1);
T_total = 0;
%statistics to store the empiricals for the reversion step
revert_calT = zeros(N,w/2);
revert_cum_theta = zeros(d,N,w/2);
n_revert =1;

%warm up
[theta_jt,jt,ptheta,Ptheta,Theta] = context_generator(CAid,ptheta,Ptheta,Theta);
[theta_next,jnext,ptheta,Ptheta,Theta] = context_generator(CAid,ptheta,Ptheta,Theta);
l = 2;
t = 1;
c_l = C_fixed(l);
%arm sample block
CA_window = 2*w/2; % double the window for context alignment for more robust performance
[~,~,Ainv_xyt,is_change] = sample_n_arm(X,Lambda,CA_window,Alambda_inv,theta_jt,t,c_l,theta_next);
if is_change % a changepoint has occurred, update the context
    jt = jnext;
    theta_jt = theta_next;
    [theta_next,jnext,ptheta,Ptheta,Theta] = context_generator(CAid,ptheta,Ptheta,Theta);
    l = l+1;
    c_l = C_fixed(l);
end
%initialization
t = w/2;
tCA = t;
CAid{jt} = Ainv_xyt;
%update PSBAI statistics
hatjt = 1;



%main loop
while 1
    t = t + 1;
    %arm sample block
    [~,~,Ainv_xyt,is_change] = sample_n_arm(X,Lambda,1,Alambda_inv,theta_jt,t,c_l,theta_next);
    if is_change % a change point occurred, update the context
        jt = jnext;
        theta_jt = theta_next;
        [theta_next,jnext,ptheta,Ptheta,Theta] = context_generator(CAid,ptheta,Ptheta,Theta);
        l = l+1;
        c_l = C_fixed(l);
    end
    %% Exp phase
    if mod(t-tCA,gamma) >0 
        %update PSBAI statistics
        calT(hatjt) = calT(hatjt) + 1;
        T_total = T_total + 1;
        cum_theta(:,hatjt) = cum_theta(:,hatjt) + Ainv_xyt;
        theta_hat(:,hatjt) = cum_theta(:,hatjt)/calT(hatjt);
        hat_p = calT./T_total;

        % stopping rule check
        [~,xtstar] = max(X'*theta_hat*hat_p);
        subarms = X_index(X_index~=xtstar);
        X_diff = X(:,xtstar)-X(:,subarms);
        rho = rho_func(delta,K,calT,T_total,d,hat_p,Lmax,N,X_diff,theta_hat,varepsilon);
        Delta_hat_diff = X_diff'*theta_hat*hat_p;
        check_procedure_1 = Delta_hat_diff - rho;
        deltad = delta/(15*N*T_total^3);
        if min(check_procedure_1) > -varepsilon && T_total > 2*Lmax/9*log(2/deltad) && isinf(tCD) %procedure one of the stopping rule
            hatx_varepsilon = xtstar;
            tCD = n_cd;
        elseif tCD == n_cd - w/2    %procedure two of the stopping rule
            disp("PSBAI finds an varepsilon best arm!")
            tau_PSBAI = t;
            toc
            flag_PSBAI = 1;
            l_PSBAI = l;
            [rho,alpha,xi,DE,beta_clip1] = rho_func(delta,K,calT,T_total,d,hat_p,Lmax,N,X_diff,theta_hat,varepsilon);
            break;
        end

    else
    %% CD phase
        %store statistics for reversion
        n_revert = mod(n_revert,w/2) + 1;
        revert_calT(:,n_revert) = calT;
        revert_cum_theta(:,:,n_revert) = cum_theta;
        %store CD sample
        n_cd = n_cd + 1;
        n_cd_mod = mod(n_cd-1,w)+1;
        CDsample(:,n_cd_mod) = Ainv_xyt;
        %enough CD sample for changepoint detection
        if n_cd >= w
            if n_cd_mod <= w/2
                cum_theta1 = sum(CDsample(:,n_cd_mod+1:n_cd_mod+w/2),2);
                cum_theta2 = sum(CDsample(:,[1:n_cd_mod,n_cd_mod+w/2+1:w]),2);
            else
                cum_theta1 = sum(CDsample(:,n_cd_mod-w/2:n_cd_mod-1),2);
                cum_theta2 = sum(CDsample(:,[1:n_cd_mod-w/2-1,n_cd_mod:w]),2);
            end
            if LCD(X,w,b,cum_theta1,cum_theta2)
                CDsample = zeros(d,w);
                n_cd = 0;
                %revert the statistics
                n_revert = mod(n_revert,w/2) + 1;
                calT = revert_calT(:,n_revert);
                cum_theta = revert_cum_theta(:,:,n_revert);
                T_total = sum(calT);
                theta_hat= cum_theta./(calT');
                hat_p = calT./T_total;
                %% Context Alignment
                %arm sample block
                [~,~,Ainv_xyt,is_change] = sample_n_arm(X,Lambda,CA_window,Alambda_inv,theta_jt,t,c_l,theta_next);
                t = t + CA_window;
                tCA = t;
                tCD = inf;
                if is_change % a change point has occurred, update the context
                    jt = jnext;
                    theta_jt = theta_next;
                    [theta_next,jnext,ptheta,Ptheta,Theta] = context_generator(CAid,ptheta,Ptheta,Theta);
                    l = l+1;
                    c_l = C_fixed(l);
                end
                [hatjt,CAid] =  LCA(X,2*CA_window,b,CAid,Ainv_xyt);
                if hatjt > N
                    disp("More than N contexts identified!")
                    break
                end
            end
        end
    end
end



end

