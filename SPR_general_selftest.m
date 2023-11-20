%% Guide
% This file provides a test sample for SPR when using random measurements.

%% Model
% y = |Ax|.^2

%% Settings
clc
clear
close all

m               = 400;              % The number of measurements
n               = 1000;              % The length of x
opNum           = 1;                % The number of optimal subproblem solvers
method          = "Gaussian";       % The generating method of x (0-1 or Gaussian)
sparsity        = 10;                % Sparsity of x
isComplex       = 1;                % If complex signal
K               = sparsity;         % The sparsity level we esitmate.
tol             = 1e-6;             % The threshold.
iterNum         = K*200;            % The maximal iteration number


%% Generate Test Sample
[X,Y,A,supportX] = init_general(n, m, sparsity, isComplex, method);
Yt = abs(Y);
[~,supportInit] = SpectralInit(Yt, A, K);
%% Recovery
[x1,k] = SPsolver_general(Yt, A, K, iterNum, opNum, isComplex, tol);

%% Test Recovery Result
phase = x1(supportX)./X(supportX);
Loss = @(A, x, y)1/m*norm(abs(A*x).^2-y,2);

if norm(x1-X*phase(1)) < tol && Loss(A, x1, abs(Y).^2) <tol
    fprintf("Succeed!\n")
else
    fprintf("Fail!\n")
end
