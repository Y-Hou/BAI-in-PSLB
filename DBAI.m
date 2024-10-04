%% D\epsilon BAI algorithm
n_epsilon = 12;
epsilon_range = 0.03*ones(n_epsilon,1);
for i=1:n_epsilon
    epsilon_range(i) = epsilon_range(i)*1.35^(i);
end
n_epsilon_range  = length(epsilon_range);
number_of_exp = 20;
table_DBAI = zeros(3,n_epsilon_range,number_of_exp);

for i=1:number_of_exp
    for j=1:n_epsilon_range
        rng(i*100+j)
        %arms
        d=2;
        I=eye(d);
        X = [eye(d),[cos(pi/8);sin(pi/8)]];
        [d,K] =size(X);
        %% contexts
        a = cos(pi/8);
        aa = sin(pi/8);
        N = d;
        Theta = zeros(d);
        Theta(:,1) = a*I(:,1) + aa*I(:,2);
        Theta(:,2) = a*I(:,1) - aa*I(:,2);
        
        %distribution
        varepsilon = epsilon_range(j);
        p = 1/N;
        ptheta = p*ones(N,1);
        Ptheta = ptheta;
        for i1=2:length(ptheta)
            Ptheta(i1) = Ptheta(i1-1) + ptheta(i1);
        end

        %other parameters
        delta =0.05;
        [d,K] = size(X);
        X_index = 1:K;
        N = length(Ptheta);
        cum_p_DBAI = zeros(N,1);
        l = 1;
        a = rand;
        [~,jt] = max(Ptheta>a);
        
        % main loop
        while 1
            l =l + 1;
            a = rand;
            [~,jt] = max(Ptheta>a);
            % update the statistics
            cum_p_DBAI(jt) = cum_p_DBAI(jt)+1;
            hatp_DBAI = cum_p_DBAI./sum(cum_p_DBAI);
            Etheta_hat_DBAI = Theta*hatp_DBAI;
            % compute the confidence radius
            [~,x_DBAI] = max(X'*Etheta_hat_DBAI);
            subarms = X_index(X_index~=x_DBAI);
            X_diff_DBAI = X(:,x_DBAI)-X(:,subarms);
            zeta = -sum(X_diff_DBAI'*Theta,2)/N;
            rho_DBAI = min(sum(abs(X_diff_DBAI'*Theta+zeta),2) * sqrt(1/(2*l)*log(2*1.2021*N*l^3/delta)) , 4);
            % stopping rule check
            if min(X_diff_DBAI'*Theta*hatp_DBAI)-rho_DBAI > -varepsilon
                l_DBAI = l;
                [~,x_DBAI] = max(X'*Etheta_hat_DBAI);
                xtstar_DBAI = x_DBAI;
                break
            end
        end
        table_DBAI (:,j,i) = [l_DBAI;xtstar_DBAI;sum(rho_DBAI)];
    end
end



