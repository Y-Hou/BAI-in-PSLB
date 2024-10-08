%% generate the change points
%It has been excuted once and generates a sequence of changepoints,
%which is stored in "change_point_sequence.mat".
Lmin = 3*10^4;
Lmax = 5*10^4;
rng(1); %fixed the generation process
l_t = 3*10^5;
C_fixed = rand(l_t,1);
C_fixed = Lmax.*(C_fixed>0.8) + Lmin.*(C_fixed<=0.8);
C_fixed = [1;C_fixed];
for i = 2:length(C_fixed)-1
    C_fixed(i) = C_fixed(i-1) + C_fixed(i);
end
C_fixed_small_Lmax = C_fixed;
cp_seq = strcat("change_point_sequence",".mat");
save(cp_seq)

