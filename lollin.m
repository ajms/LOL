function [ cost, dcost ] = lollin ( phi, I, lambda1, lambda2, nu )
    u = I(:);
    cost = 0;
    [M, N] = size(u);
    dcost = zeros(M*N,1);
    mu1 = mean(u(phi>=0));
    mu2 = mean(u(phi<0));
    mom = 2;
    for i = 1:M*N
        t = phi(i);
        iv = u(i);
        tj = floor(t);
        p = [(t-tj) (tj+1-t)];
        dp = [-1 1];
        for k = 0:1
            if tj+k >= 0
                cost = cost + lambda1*(iv-mu1)^mom * p(k+1) + nu * p(k+1);
                dcost(i) = dcost(i) + lambda1*abs(iv-mu1)^mom * dp(k+1) + nu * dp(k+1);
            elseif tj+k < 0
                cost = cost + lambda2*(iv-mu2)^mom*p(k+1);
                dcost(i) = dcost(i) + lambda2*abs(iv-mu2)^mom * dp(k+1);
            end
        end
    end
end
