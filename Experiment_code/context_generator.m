function [theta_jt,jt,ptheta,Ptheta,Theta] = context_generator(CAid,ptheta,Ptheta,Theta)
%context_generator
%generate a new context according to Ptheta
a = rand;
[~,jt] = max(Ptheta>a);
theta_jt = Theta(:,jt);
seen = length(CAid);
%if this context did not appear before, label it as the len(seen)+1; and
%update the context labels in CPF Ptheta
if jt > seen
    pjt = ptheta(jt);
    thetajt = Theta(:,jt);
    ptheta(jt) = ptheta(seen+1);
    Theta(:,jt) = Theta(:,seen+1);
    ptheta(seen+1) = pjt;
    Theta(:,seen+1) = thetajt;
    Ptheta = ptheta;
    for i=2:length(ptheta)
        Ptheta(i) = Ptheta(i-1) + ptheta(i);
    end
    jt = seen + 1;
end




end

