function [xhat,k, mk] = BB(Loss, Grad, A, x, y, support)
%% The BB method
n           = length(x);
maxIter     = 200; %
alpha       = 0.9; % 
sigma       = 0.8; % 
epsilon     = 1e-10; %
A0          = A(:, support); %
x0          = x(support,1); %


x2 = x0;
x1 = x0;
G2 = 0;

mk = 0;
for k = 1:maxIter
    G1 = G2;
    G2 = feval(Grad, A0, x2, y);
    if norm(G2) < epsilon
        break
    end
    deltaG = G2-G1;
    deltaX = x2-x1;
    if k == 1
        m = 0;
        mk = 0;
        while m<20
            if feval(Loss, A0, x2+sigma^m*(-G2), y)<feval(Loss, A0, x2, y)+ ...
                    alpha*sigma^m*G2'*(-G2)
                mk = m;
                break
            end
            m = m+1;
        end
        BBstep = sigma^mk*(-G2);
    elseif mod(k,4) < 2 
        BBstep = -(deltaX'*deltaX)/(deltaX'*deltaG)*G2;
    else
        BBstep = -(deltaX'*deltaG)/(deltaG'*deltaG)*G2;
    end
    
    x1 = x2;
    x2 = x1+BBstep;
    
    
end

xhat = zeros(n,1);
xhat(support,1) = x2;

end

