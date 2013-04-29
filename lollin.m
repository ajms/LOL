function [ cost, dcost ] = lollin ( phi, I, lambda1, lambda2, nu )
    u = I(:);
    cost = 0;
    [M, N] = size(u);
    dcost = zeros(M*N,1);
    mu1 = mean(u(phi>=0));
    mu2 = mean(u(phi<0));
    mom = 2;
    for i = 2:M*N-1
        phiv = phi(i);
        iv = u(i);
        basis = floor(phiv);
        t = phiv-basis;
        p = [1-t t];
        dp = [-1 1];
        for k = 0:1
            if basis+k >= 0
                cost = cost + (iv-mu1)^2 * p(k+1) + nu * p(k+1);
                dcost(i,1) = dcost(i,1) + lambda1*abs(iv-mu1)^mom * dp(k+1) + nu * dp(k+1);
            elseif basis+k < 0
                cost = cost+ lambda2*(iv-mu2)^2*p(k+1);
                dcost(i,1) = dcost(i,1) + lambda2*abs(iv-mu2)^mom * dp(k+1);
            end
        end
    end
end

