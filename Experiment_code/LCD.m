function alarm = LCD(X,w,b,cum_theta1,cum_theta2)
%Linear Change Detection

theta1 = 2/w*cum_theta1;
theta2 = 2/w*cum_theta2;

ind = abs(X'*(theta2-theta1))>b;

if sum(ind)>0
    alarm = 1;
else
    alarm = 0;
end

end