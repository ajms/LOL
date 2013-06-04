function [ F, dF, ddF ] = olnarrow ( phi, I, lambda1, lambda2, mu, nu )
% Function to calculate the value of the functional F and the
% derivatives in a narrowband of the 0'th levelset.
    
%      [ F, dF ] = olnarrow ( phi, I, lambda1, lambda2, mu, nu ) returns
%      a scalar value of the functional F and a matrix of
%      derivatives in every point of the image
%
%      Input: 
%      - phi: levelset function
%      - I: image
%      - lambda1,lambda2: parameters to weight the inner and outer
%      segments respectively
%      - mu: parameter to weight the length regularization
%      - nu: parameter to weight the area regularization

    % define variables
    F = 0;
    [M, N] = size(phi);
    dF = zeros(M*N,1);
    mu1 = mean(I(phi>=0));
    mu2 = mean(I(phi<0));
    mom = 2;
    
    % iteration through the image calculating the value of F and
    % the derivatives in a narrowband of the 0'th levelset.
    for i = 1:M*N
        phiv = phi(i);
        iv = I(i);
        basis = floor(phiv);
        t = phiv-basis;
        p = 1/6*[-t^3+3*t^2-3*t+1 3*t^3-6*t^2+4 -3*t^3+3*t^2+3*t+1 t^3];
        dp = 1/6*[-3*t^2+6*t-3 9*t^2-12*t -9*t^2+6*t+3 3*t^2];
        for k = -1:2
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