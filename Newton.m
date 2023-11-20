function [xhat,k] = Newton(Grad, Hess, A, x, y, support)
%% The Newton method
n = length(x);
maxIter = 1000; %  
epsilon = 1e-10; % 
A0 = A(:, support); % 
x0 = x(support,1); % 
for k = 1:maxIter
    x1 = x0;
    Gk = feval(Grad, A0, x1, y);
    Hk = feval(Hess, A0, x1, y);
    Dk = -Hk\Gk;
    if norm(Dk,1) < epsilon
        break
    end
    x0 = x1+1*Dk;
end

xhat = zeros(n,1);
xhat(support,1) = x0;
end

