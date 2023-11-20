function [x1, k] = SPsolver_general(Yt, A, K, iterNum, opNum, isComplex, tol, InitialGuess)
%% Introduction
% Input
% Y = |Ax|.^2
% Y:          The measurement signal
% A:          Random measurement matrix
% K:          Estimated sparisty (not less than half real sparisty)
% iterNum:    The maximal iteration number
% opNum:      The number of subproblem solvers.
% isComplex:  Complex or not
% tol:        The maximal error that we can tolerate.

% Output
% x1:         Solver
% k:          Itration number

% **Remark**
% We strongly recommend that the maximal intensity of the orignal signal x
% does not larger than 1. Or the distribution of x obeys approximately the
% standard normal distribution.

%% Tool function.
[m,n]= size(A);
Loss = @(A, x, y)1/m*norm(abs(A*x).^2-y,2);
Grad = @(A, x, y)2/m*A'*((abs(A*x).^2-y).*(A*x));
Hess = @(A, x, y)4/m*A'*((abs(A*x).^2).*A) + A'*((abs(A*x).^2-y).*A);

%% The algorithm begins


K0 = 3*K;
K1 = K;
K2 = K0-K1;

[~, n] = size(A);

x0 = SpectralInit(Yt, A, K);
Yt = Yt.^2;


grad0 = Grad(A, x0, Yt);
[~, indexes] = sort(abs(grad0));
support = indexes(end-K2+1:end); %

x1 = x0;

%% Loops
repeat = 0; % trick

for k = 1:iterNum
    [~, indexes] = sort(abs(Grad(A, x1, Yt)));
    
    supportcopy = support; %
    
    support0 = indexes(end-K2+1:end);
    supportX = union(support, support0); % update
    Xhat = [];
    R00 = [];
    for i = 1:opNum
        x2 = randn(n,1)+1i*isComplex*randn(n,1);
        %[xhat0,k0] = Newton(Grad, Hess, A, x2, Yt, supportX);
        [xhat0,~] = BB(Loss, Grad, A, x2, Yt, supportX); % Recommended
        %[xhat0,indicator] = PGD(Loss, Grad, A, x2, Yt, supportX);
        Xhat = [Xhat xhat0];
        r00 = norm(abs(A*xhat0).^2-Yt);
        R00 = [R00 r00];
    end
    [~, index] = sort(R00);
    xhat = Xhat(:,index(1)); %
    
    [~, indexes] = sort(abs(xhat));
    support = indexes(end-K1+1:end); %
    
    
    % trick
    %     if isempty(setxor(support,supportcopy))
    %         repeat = repeat + 1;
    %         % p = [p 1];
    %         if repeat >= 5
    %             try
    %                 support = [indexes(end-K); indexes(end-K+2:end)];
    %             catch
    %                 s = randperm(n);
    %                 support = s(1:K);
    %             end
    %             repeat = 0;
    %         end
    %     else
    %         % p = [p 0];
    %         repeat = 0;
    %     end
    x1 = xhat; % different from Cosamp
    
    %x0 = xhat;
    %projM = zeros(n,1);
    %projM(support, 1) = 1;
    %x1 = x0.*projM; %
    %r0 = Yt - LinearOperatorA(A, x1*conj(x1).');
    
    %stopping
    if Loss(A,x1,Yt) < tol
        %         x0 = x1;
        %         projM = zeros(n,1);
        %         projM(support, 1) = 1;
        %         x1 = x0.*projM;
        break
    end
    
end

end

