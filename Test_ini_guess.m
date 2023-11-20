%% Guide
% This file provides a test sample for SPR when using random measurements.

%% Model
% y = |Ax|.^2

%% Settings


m               = 80;               % The number of measurements
n               = 100;              % The length of x
opNum           = 1;                % The number of optimal subproblem solvers
method          = "Gaussian";       % The generating method of x (0-1 or Gaussian)
sparsity        = 5;                % Sparsity of x
isComplex       = 0;                % If complex signal
K               = sparsity;         % The sparsity level we esitmate.
tol             = 1e-6;             % The threshold.
iterNum         = K*200;            % The maximal iteration number


testNum         = 10000;
Record          = zeros(testNum, 3);
for i=1:testNum
    if mod(i, 100) == 0
        fprintf('Now WF: %d\n', i);
    end
InitialGuess    = "WF_Candes";         % "Null"/"GN_XZQ"/"WF_Candes"

%% Generate Test Sample
[X,Y,A,supportX] = init_general(n, m, sparsity, isComplex, method);
Yt = abs(Y).^2;

%% Recovery
[x1,k] = SPsolver_general(Yt, A, K, iterNum, opNum, isComplex, tol, InitialGuess);

%% Test Recovery Result
phase = x1(supportX)./X(supportX);
Loss = @(A, x, y)norm(abs(A*x).^2-y,2);

if norm(x1-X*phase(1)) < tol && Loss(A, x1, abs(Y).^2) <tol
    Record(i,1) = k;
else
    Record(i,1) = -1;
end
end

fprintf('\n\n');

for i=1:testNum
    if mod(i, 100) == 0
        fprintf('Now GN: %d\n', i);
    end
InitialGuess    = "GN_XZQ";         % "Null"/"GN_XZQ"/"WF_Candes"

%% Generate Test Sample
[X,Y,A,supportX] = init_general(n, m, sparsity, isComplex, method);
Yt = abs(Y).^2;

%% Recovery
[x1,k] = SPsolver_general(Yt, A, K, iterNum, opNum, isComplex, tol, InitialGuess);

%% Test Recovery Result
phase = x1(supportX)./X(supportX);
Loss = @(A, x, y)norm(abs(A*x).^2-y,2);

if norm(x1-X*phase(1)) < tol && Loss(A, x1, abs(Y).^2) <tol
    Record(i,2) = k;
else
    Record(i,2) = -1;
end
end

fprintf('\n\n');

for i=1:testNum
    if mod(i, 100) == 0
        fprintf('Now Null: %d\n', i);
    end
InitialGuess    = "Null";         % "Null"/"GN_XZQ"/"WF_Candes"

%% Generate Test Sample
[X,Y,A,supportX] = init_general(n, m, sparsity, isComplex, method);
Yt = abs(Y).^2;

%% Recovery
[x1,k] = SPsolver_general(Yt, A, K, iterNum, opNum, isComplex, tol, InitialGuess);

%% Test Recovery Result
phase = x1(supportX)./X(supportX);
Loss = @(A, x, y)norm(abs(A*x).^2-y,2);

if norm(x1-X*phase(1)) < tol && Loss(A, x1, abs(Y).^2) <tol
    Record(i,3) = k;
else
    Record(i,3) = -1;
end
end

save('Record.mat', 'Record')