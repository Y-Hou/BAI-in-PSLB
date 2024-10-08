%% Setup the epsilon's
n_epsilon = 12;
epsilon_range = 0.03*ones(n_epsilon,1);
for i=1:n_epsilon
    epsilon_range(i) = epsilon_range(i)*1.35^(i);
end
n_epsilon_range  = length(epsilon_range);

number_of_exp = 20; %number of repetitions
table = zeros(16,n_epsilon_range,number_of_exp); % used to store the result

%% arms
d = 2;
I=eye(d);
X = [eye(d),[cos(pi/8);sin(pi/8)]];
K = size(X,2);
%% contexts
a = cos(pi/8);
aa = sin(pi/8);
N = d;
Theta = zeros(d);
Theta(:,1) = a*I(:,1) + aa*I(:,2);
Theta(:,2) = a*I(:,1) - aa*I(:,2);

%% distribution
p = 1/N;
ptheta = p*ones(N,1);
Ptheta = ptheta;
for i1=2:length(ptheta)
    Ptheta(i1) = Ptheta(i1-1) + ptheta(i1);
end
%change points
load("change_point_sequence.mat",'C_fixed');

%others
delta = 0.05;

% groundtruth of Lmin and Lmax
Lmin_gt = 30000;
Lmax_gt = 50000;

for mismatch_index = [1,3]

    % mismatch in Lmin and Lmax
    mismatch = [0.8,1,1.2];
    % mismatch_index = 3;
    Lmin = Lmin_gt*mismatch(mismatch_index); %modify the index for different levels of misspecification
    Lmax = Lmax_gt*mismatch(mismatch_index);
    
    gamma = 6;
    w = floor(Lmin/gamma/3);
    w = w + mod(w,2);
    
    
    for j=n_epsilon_range:-1:1
        varepsilon = epsilon_range(j);
        tau_temp = 38400*log(80)*N*Lmax/min(varepsilon,1)^2*log(N^2*K*Lmax/(min(varepsilon,1)^2*delta));
        delta_FAE = gamma*delta/(4*tau_temp^2*K);
        delta_FA = Lmax*delta/(4*N*tau_temp^2);
        b = d/w*log(2/delta_FAE) + sqrt((d/w*log(2/delta_FAE))^2+ d/w*log(2/delta_FAE));
        for i=1:number_of_exp
            disp(strcat('---current epsilon index:',num2str(j),', current experiment trial:',num2str(i),'---'));
            seed = i*100+j+1000;
            % revoke the PS\epsilon BAI algorithm
            [tau_PSBAI,hatx_varepsilon,flag_PSBAI,l_PSBAI,tau_naive,x_dot_varepsilon,flag_naive,l_naive,rho,alpha,xi,DE,beta_clip1] = PSBAI_func(X,Theta,ptheta,Ptheta,Lmin_gt,Lmax,varepsilon,w,gamma,b,delta,C_fixed,seed);
            % store the results
            table(:,j,i) = [tau_PSBAI,hatx_varepsilon,flag_PSBAI,l_PSBAI,tau_naive,x_dot_varepsilon,flag_naive,l_naive,rho',alpha,xi,DE',beta_clip1']';
        end
        disp(table(2,:,:))
        %save the intermidiate result for easy checking
        filename = strcat('intermidiate_',num2str(j),'.mat');
        save(filename);
    end
    table
    if mismatch_index == 1
        table0_8 =  table;
    elseif mismatch_index == 3
        table1_2 =  table;
    else
        table1_0 =  table;
    end
    filename = strcat('trial_',num2str(mismatch_index),'.mat');
    save(filename)
end    

