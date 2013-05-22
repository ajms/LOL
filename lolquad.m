function [ F, dF, ddF ] = lolquad ( phi, I, lambda1, lambda2, mu, nu )
    F = 0;
    [M, N] = size(phi);
    dF = zeros(M*N,1);
    mu1 = mean(I(phi>=0));
    mu2 = mean(I(phi<0));
    mom = 2;

    for i = 1:M*N
        phiv = phi(i);
        iv = I(i);
        basis = floor(phiv);
        t = phiv-basis;
        p = 1/2*[t^2-2*t+1 -2*t^2+2*t+1 t^2];
        dp = 1/2*[2*t-2 -4*t+2 2*t];
        for k = -1:1
            if basis+k == 0
                F = F + mu * p(k+2);
                dF(i) = dF(i) + mu*dp(k+2);
            end
            if basis+k >= 0
                F = F + lambda1*abs(iv-mu1)^mom * p(k+2) + nu * p(k+2);
                dF(i) = dF(i) + lambda1*abs(iv-mu1)^mom * dp(k+2) + nu * dp(k+2);
            else
                F = F + lambda2*abs(iv-mu2)^mom * p(k+2);
                dF(i) = dF(i) + lambda2*abs(iv-mu2)^mom * dp(k+2);
            end
        end
    end
end
