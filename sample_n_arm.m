function [Xt,Yt,Ainv_XYt,is_change] = sample_n_arm(X,p,n,Alambda_inv,theta_t,t,c_l,theta_next)
%samplel n arms according to the distribution p, p is in cdf.
if n+t-1 < c_l
    a = rand(1,n);
    P = repmat(p,[1,n]);
    [~,Xt] = max(P>a);
    Yt = X(:,Xt)'*theta_t + ClipG(3,n);
    XYt = X(:,Xt)*Yt;
    is_change = 0;
else
    %there is a change point c_l during the sampling
    %sample from the first context 
    m = n;
    n = c_l-t;
    a = rand(1,n);
    P = repmat(p,[1,n]);
    [~,Xt1] = max(P>a);
    Yt1 = X(:,Xt1)'*theta_t + ClipG(3,n);
    XYt1 = X(:,Xt1)*Yt1;
    %sample from the second context
    n = m - n;
    a = rand(1,n);
    P = repmat(p,[1,n]);
    [~,Xt] = max(P>a);
    Yt = X(:,Xt)'*theta_next + ClipG(3,n);
    XYt = X(:,Xt)*Yt;
    %assemble the results
    Xt =[Xt1,Xt];
    Yt = [Yt1;Yt];
    XYt = XYt1 + XYt;
    is_change = 1;
end
Ainv_XYt = Alambda_inv*XYt;
end

