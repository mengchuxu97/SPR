function [x0, supportInit] = SpectralInit(Yt, A, K)
%%
%Yt = |Ax|
%%
[~, n] = size(A);
Zj = Yt.*abs(A);
Z = sum(Zj,1);
[~, index] = sort(Z, 'descend');
x0 = zeros(n,1);
supportInit = index(1:K);
x0(index(1:K)) = randn(K,1);
end