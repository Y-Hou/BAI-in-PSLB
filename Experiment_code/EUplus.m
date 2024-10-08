function [cum_tiltheta,flag,x_dot_varepsilon] = EUplus(Ainv_xyt,t,cum_tiltheta,tstar,C3,X,Lmax,K,delta,d,varepsilon,X_index)
%Estimates Update: update the statistics for the naive algorithm, and check
%the stopping rule for the naive algorithm. 
%If it stops, the indicator flag is changed and the empirically best arm is
%recommended.
cum_tiltheta = cum_tiltheta + Ainv_xyt;
tiltheta = cum_tiltheta/t;
x_dott = zeros(d,1);
tilrho = 2*varepsilon;
if t > tstar
    [~,x_dott] = max(X'*tiltheta);
    tilrho = sqrt(8*Lmax/t*log(4*K*C3*t^3/delta)) + 5*sqrt(d/t*log(4*K*C3*t^3/delta));
end
subarms = X_index(X_index~=x_dott);
x_dot_varepsilon = 0;
flag = 0;
if X(:,x_dott)'*tiltheta-2*tilrho+varepsilon >= max(X(:,subarms)'*tiltheta)
    x_dot_varepsilon = x_dott;
    flag = 2;
    disp("Naive alg finds an varepsilon-best arm")
end

end

