close all;
if ndims(I) == 3
    [M, N, O] = size(I);
else
    [M, N] = size(I);
end

% parameters
no = '10';
maxit = 100;
gamma = 10;
lambda1 = 1;
lambda2 = 1;
mu = 0;%10^(-4);
nu = 0;
sigma = 0.5;
m = 70;
options.Method = 'lbgfs';

% smoothen image using scalen with sigma
if sigma == 0
    Ismooth = I;
else
    if ndims(I) == 3
        Ismooth = real(ifftn(scalen(fftn(I),[sigma,sigma,sigma],[0,0,0])));
    else
        Ismooth = real(ifftn(scalen(fftn(I),[sigma,sigma],[0,0])));
    end
end

%initialize phi
if ndims(I) == 3
    phi = -ones(M,N,O);
    [X Y Z] = meshgrid(1:M);
    phip = (X-floor(M/2)).^2 + (Y-floor(N/2)).^2 + (Z-floor(O/2)).^2;
    phi(phip <= m^2) = 1;
    phi = init(phi);
else
    phi = -ones(M,N);
    [X Y] = meshgrid(1:M);
    phip = (X-50).^2 + (Y-0).^2; 
    phi(phip <= m^2) = 1;
    phi = init(phi);
    initial = figure;
    hold on;
    imagesc(I);
    colorbar();
    colormap('gray');
    contour(reshape(phi,M,N), [0 0], 'Color', [1 0 0],'LineWidth',3);
    axis tight;
    print(initial,'-dpsc',strcat('I',no,'init.eps'));
end

if mu == 0
    phi = phi/(max(phi(:))-min(phi(:)));
end

F = zeros(maxit+1,1);
diffF = zeros(maxit+1,1);
normdF = zeros(maxit+1,1);
dF = 0;
tic;
for i=1:maxit
    fprintf('Iteration: %d, Diff in F: %f, Norm dphi: %f\n', i, diffF(i), normdF(i));
    if mu > 0
        [phi_n, F(i+1)] = minFunc(@lolquad,phi(:),options,Ismooth,lambda1,lambda2,mu,nu);
        phi = init(phi_n);
    else
        [F(i+1), dF] = lolquad(phi(:),Ismooth,lambda1,lambda2,mu,nu);
        phi = phi(:)-gamma*dF;
        normdF(i+1) = norm(gamma*dF);
    end
    diffF(i+1) = abs(F(i+1)-F(i));
    if  diffF(i+1) < 0.05 &&  normdF(i+1) < 0.05
       fprintf(['Done after %d iterations!\n Diff in F: %f, Norm %f\n'],i,diffF(i+1),normdF(i+1));
       break;
    end
end
toc;

if ndims(I) == 2
    seg = figure;
    hold on;
    imagesc(I);
    colorbar();
    colormap('gray');
    contour(reshape(phi,M,N), [0 0], 'Color', [1 0 0],'LineWidth',3);
    axis tight;
    hold off;
    print(seg,'-dpsc',strcat('I',no,'seg.eps'));
end

Fplot = figure;
plot(1:sum(F>0),F(F>0));
set(gca,'xtick',0:1:maxit);
title('Convergence plot','interpreter','latex','FontSize',15);
xlabel('Iteration no.','interpreter','latex','FontSize',15);
ylabel('Value of $F$','interpreter','latex','FontSize',15);
print(Fplot,'-dpsc',strcat('I',no,'con.eps'));

phiplot = figure;
plot(1:sum(normdF>0),gamma*normdF(normdF>0));
set(gca,'xtick',0:1:maxit);
title('Convergence plot','interpreter','latex','FontSize',15);
xlabel('Iteration no.','interpreter','latex','FontSize',15);
ylabel('Difference in $\varphi$','interpreter','latex','FontSize',15);
print(phiplot,'-dpsc',strcat('I',no,'con2.eps'));

dev = sum(abs((phi>=0)-Iref(:)))/(M*N);
fprintf('Deviation in pct. %f\n',dev);