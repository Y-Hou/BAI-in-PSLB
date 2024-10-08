function [noise] = ClipG(sigma,n)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
a = randn(n,1)./sigma;
noise = min(max(a,-1),1);
end

