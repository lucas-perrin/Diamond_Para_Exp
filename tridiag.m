function [M] = tridiag(a,b,c,n)

M = a*diag(ones(n-1,1),-1) + b*diag(ones(n,1),0) + c*diag(ones(n-1,1),+1);