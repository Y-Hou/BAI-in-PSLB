%% D\epsilon BAI_\beta algorithm
n_epsilon = 12;
epsilon_range = 0.03*ones(n_epsilon,1);
for i=1:n_epsilon
    epsilon_range(i) = epsilon_range(i)*1.35^(i);
end
n_epsilon_range  = length(epsilon_range);
number_of_exp = 20;
table_DBAIbeta = zeros(3,n_epsilon_range,number_of_exp);
load("change_point_sequence.mat");

for i=1:number_of_exp
    for j=1:n_epsilon_range
        rng(100*i+j);

        % arms
        d=2;
        I=eye(d);
        X = [eye(d),[cos(pi/8);sin(pi/8)]];
        [d,K] =size(X);

        % contexts
        a = cos(pi/8);
        aa = sin(pi/8);
        N = d;
        Theta = zeros(d);
        Theta(:,1) = a*I(:,1) + aa*I(:,2);
        Theta(:,2) = a*I(:,1) - aa*I(:,2);

        % distribution
        varepsilon = epsilon_range(j);
        p = 1/N;
        ptheta = p*ones(N,1);
        Ptheta = ptheta;
        for i1=2:length(ptheta)
            Ptheta(i1) = Ptheta(i1-1) + ptheta(i1);
        end
        
        %others
        delta = 0.05;
        [d,K] = size(X);
        X_index = 1:K;
        N = length(Ptheta);
        cum_p_DBAIbeta = zeros(N,1);
        l = 1;
        a = rand;
        [~,jt] = max(Ptheta>a);

        while 1
            l =l + 1;
            a = rand;
            [~,jt] = max(Ptheta>a);
            % update the statistics
            Lt = C_fixed(l+1)- C_fixed(l);
            cum_p_DBAIbeta(jt) = cum_p_DBAIbeta(jt) + Lt;
            T_total = sum(cum_p_DBAIbeta);
            hatp_DBAIbeta = cum_p_DBAIbeta./T_total;
            deltad = delta/(1.2021*N*T_total^3);
            phi = min(4*max(hatp_DBAIbeta,25/4*Lmax/T_total*log(2/deltad)),1/4);
            tilbeta2 = 2.*phi*Lmax/T_total*log(2/deltad);
            beta_clip1 = min(1/3*Lmax/T_total*log(2/deltad)+sqrt((1/3*Lmax/T_total*log(2/deltad))^2+tilbeta2),1);
            % update the empirically best arm andt the confidence radii
            [~,x_DBAIbeta] = max(X'*Theta*hatp_DBAIbeta);
            subarms_2 = X_index(X_index~=x_DBAIbeta);
            X_diff_DBAIbeta = X(:,x_DBAIbeta)-X(:,subarms_2);
            zeta = -sum(X_diff_DBAIbeta'*Theta,2)/N;
            rho_DBAIbeta = sum(abs(X_diff_DBAIbeta'*Theta+zeta),2).*beta_clip1;
            if min(X_diff_DBAIbeta'*Theta*hatp_DBAIbeta-rho_DBAIbeta) > -varepsilon && T_total > 2*Lmax/9*log(2/deltad)
                l_DBAIbeta = l;
                xtstar_DBAIbeta = x_DBAIbeta;
                break
            end
        end
        table_DBAIbeta (:,j,i) = [l_DBAIbeta;xtstar_DBAIbeta;sum(rho_DBAIbeta)];
    end
end



