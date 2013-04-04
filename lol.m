function [ cost, dcost ] = lol ( phi )

    RGB = imread('test_images/test2.jpg');
    u0 = double(rgb2gray(RGB))/255;
    u = u0(:);
    cost = 0;
    dcost = 0;
    [M, N] = size(u);
    pointm = floor(min(phi(:))):-1;
    pointp = 0:ceil(max(phi(:)));
    histm = zeros(size(pointm));
    histp = zeros(size(pointp));
    mu1 = mean(u(phi>=0));
    mu2 = mean(u(phi<0));
    mu = 0;
    nu = 1;

    for i = 1:M*N
        phiv = phi(i);
        iv = u(i);
        basis = floor(phiv);
        t = phiv-basis;
        p = [1/6*(1-t)^3 2/3-1/2*t^2*(2-t) 2/3-1/2*(1-t)^3 1/6*t^3];
        dp = [-3/6*(1-t)^2 3/2*t^2-2*t 3/2*(1-t)^2 3/6*t^2];
        for k = 0:3
            if i+k == 1 || i+k == M*N 
                continue
            end
            if basis+k == 0
                cost = cost + mu*(iv-mu)^2 * p(k+1);
            end
            if basis+k >= 0
                cost = nu*cost+(iv-mu1)^2*p(k+1);
                %histp(basis+k-1) = histp(basis+k-1) + p(k+1);
                dcost = dcost + (iv-mu1)^2 * dp(k+1);
            else
                cost = cost+(iv-mu2)^2*p(k+1);
                %histm(abs(floor(min(phi(:))))+basis+k-1) = histm(abs(floor(min(phi(:))))+basis+k-1) + p(k+1);
                dcost = dcost + (iv-mu2)^2 * dp(k+1);
            end
        end
    end
    %{
    hold on;
    plot(pointm,histm,'r');
    plot(pointp,histp,'b');
    hold off;
    %}
end
