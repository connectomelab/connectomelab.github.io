function [m] = NLA_growth(A, d, L)
% Author: Roman Bauer

% estimate for maximal number of stages
t_max = round((log((L*(d-1))+1)/log(d)));

m = zeros(L,L);
n = 1;

for t = 1:t_max
    k = round( d^t ); % nodes to be added at the current step
    if (abs(d-1)<0.00001)
        k = 1;
    end
    
    if (k>=1)  
        if (n+k>L)
            k=L-n;
            t=t_max;
        end
        for i = 1:k
           curr_out = round(A+(0.5*A)*randn);
           m(n+i,randperm(n,min(max(0,curr_out),n))) = 1;
        end
        n = n + k;
    end
end;

end
