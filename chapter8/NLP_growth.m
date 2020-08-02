function [m] = NLP_growth(p, d, L)
% Author: Roman Bauer

% estimate for maximal number of stages
t_max = round((log((L*(d-1))+1)/log(d)));

m = zeros(L,L);
n = 1;

tol_k = 0.00001;

for t = 1:t_max
    k = round( d^t ); % nodes to be added at the current step
    if (abs(d-1)<tol_k)
        k = 1; 
    end
    if (k>(1+tol_k))
        if (n+k>L)
            k=L-n;
            maxlength = n+k;
            t=t_max;
        end
        
        m(n+1:(n+k),1:n) = rand(k, n) < p;
        n = n + k;
    end
end;
end
