function [ cost, dcost ] = lol ( phi )

    RGB = imread('test_images/test2.jpg');
    u0 = double(rgb2gray(RGB))/255;
    u = real(ifftn(scalen(fftn(u0),[2,2],[0,0])));
    u = u(:);
    cost = 0;
    [M, N] = size(phi);
    dcost = zeros(M*N,1);
    mu1 = mean(u(phi>=0));
    mu2 = mean(u(phi<0));
    mu = 0;
    nu = 1;
    mom = 3;

    for i = 1:M*N
        phiv = phi(i);
        iv = u(i);
        basis = floor(phiv);
        t = phiv-basis;
        p = [-t^3+3*t^2-3*t+1 3*t^3-6*t^2+4 -3*t^3+3*t^2+3*t+1 t^3];
        dp = [-3*t^2+6*t-3 9*t^2-12*t -9*t^2+6*t+3 3*t^2];
        for k = -1:2
            if basis+k == 0
                cost = cost + mu * p(k+2);
            end
            if basis+k >= 0
                cost = cost + nu*abs(iv-mu1)^mom * p(k+2);
                dcost(i)= nu*abs(iv-mu1)^mom * dp(k+2);
            else
                cost = cost + abs(iv-mu2)^mom * p(k+2);
                dcost(i) = abs(iv-mu2)^mom * dp(k+2);
            end
        end
    end
end
