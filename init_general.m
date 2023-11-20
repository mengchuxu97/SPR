function [X, Y, A, supportX] = init_general(n, m, sparsity, isComplex, method)
%% Introduction
% It is used to produce the case.
% Y = AX

% Input
% n,m:        A is a (m x n)-matrix
% sparsity:   The sparsity of the orignal signal x
% isComplex:  Complex or not
% method:     Standard Gaussian initialization or 0-1 initialization.

% Output
% X:          The orignal signal x
% Y:          The measurements signal
% A:          The measurement matrix obeying standard Gaussian distribution
% supportX:   The support set of X.


if isComplex
    A = (randn(m, n)+1i*randn(m, n))/sqrt(2);
else
    A = randn(m,n);
end
switch method
    case "Gaussian"
        x = (randn(n,1)+1i*isComplex*randn(n,1));
    case "0-1"
        x = ones(n,1);
end

loc = randperm(n);
supportX = loc(1:sparsity);
x(loc(sparsity+1:end)) = 0;
X = x;
Y = A*X;
end

