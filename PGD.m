function [xhat, indicator] = PGD(Loss, Grad, A, x, y, support)
%% The BB method
n           = length(x);
k           = length(support);

T           = 200; %
eta         = 1;
i_thre      = 10;
g_thre      = 1e-6; %
f_thre      = 1e-7;
i_n         = -i_thre - 1;
r           = 1e-3;

indicator   = 0;
A0          = A(:, support); %
x0          = x(support,1); %
x1          = x0;
fx0 = feval(Loss, A0, x0, y);
fx1 = fx0;

for i = 1:T
    g = feval(Grad, A0, x0, y);
    gnorm = norm(g);
    if gnorm < g_thre && i-i_n > i_thre
        i_n = i;
        xi  = randn(k,1);
        xi  = r*xi/norm(xi);
        x1  = x0;
        fx1 = fx0;
        x0  = xi + x0;
        fx0 = feval(Loss, A0, x0, y);
    end
    
    if i-i_n == i_thre && fx0 - fx1 > -f_thre
        x0 = x1;
        indicator = 1;
        break
    end
    x0 = x0 - eta*g;
end
xhat = zeros(n,1);
xhat(support,1) = x0;
normg = norm(g);
