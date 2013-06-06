close all;

% Import image
%I = rgb2gray(imread('/path/to/image.jpg'));

% Check image dimensions.
if ndims(I) == 3
    [M, N, O] = size(I);
else
    [M, N] = size(I);
end

% Parameters
no = '0'; % number of image for figures
maxit = 100; % maximal number of iterations
gamma = 10; % size of steepest descent step
lambda1 = 1; % parameter weigting the inner segment
lambda2 = 1; % parameter weigting the outer segment
mu = 0; % parameter weighting the length regularization
nu = 0; % parameter weighting the area regularization
sigma = 1; % parameter for Gaussian kernel
m =50; % size of circle for initial phi
options.Method = 'lbgfs'; % optimization method for calculations
                          % using the length regularization

% Smoothen image using scalen with sigma
if sigma == 0
    Ismooth = I;
else
    if ndims(I) == 3
        Ismooth = real(ifftn(scalen(fftn(I),[sigma,sigma,sigma],[0,0,0])));
    else
        Ismooth = real(ifftn(scalen(fftn(I),[sigma,sigma],[0,0])));
    end
end

% Initialize phi
if ndims(I) == 3
    phi = -ones(M,N,O);
    [X Y Z] = meshgrid(1:M);
    phip = (X-floor(M/2)).^2 + (Y-floor(N/2)).^2 + (Z-floor(O/2)).^2;
    phi(phip <= m^2) = 1;
else
    phi = -ones(M,N);
    [X Y] = meshgrid(1:M);
    phip = (X-floor(M/2)).^2 + (Y-floor(M/2)).^2; 
    phi(phip <= m^2) = 1;
    phi = init(phi);
    % Plot of initial segmentation.
    initial = figure;
    hold on;
    imagesc(I);
    colorbar();
    colormap('gray');
    contour(reshape(phi,M,N), [0 0], 'Color', [1 0 0],'LineWidth',3);
    axis tight;
    set(gca,'YDir','Reverse');
    print(initial,'-dpsc',strcat('I',no,'init.eps'));
end

% Normalizing the levelset function to be between -1 and 1
if mu == 0
    phi = phi/(max(phi(:))-min(phi(:)));
end

% Initialize variables
F = zeros(maxit+1,1);
diffF = zeros(maxit+1,1);
normdF = zeros(maxit+1,1);
dF = 0;

% Iterating until method has converged
tic;
for i=1:maxit
    fprintf('Iteration: %d, Diff in F: %f, Norm dphi: %f\n', i, diffF(i), normdF(i));
    if mu > 0
        % Using optimization function minFunc in case of length term
        [phi_n, F(i+1)] = minFunc(@olnarrow,phi(:),options,Ismooth,lambda1,lambda2,mu,nu);
        phi = init(reshape(phi_n,M,N));
    else
        % Using step in the steepest descent direction
        [F(i+1), dF] = ol(phi(:),Ismooth,lambda1,lambda2,mu,nu);
        phi = phi(:)-gamma*dF;
        normdF(i+1) = norm(gamma*dF);
    end
    diffF(i+1) = abs(F(i+1)-F(i));
    if  diffF(i+1) < 0.05 &&  normdF(i+1) < 0.1
       fprintf(['Done after %d iterations!\n Diff in F: %f, Norm %f\n'],i,diffF(i+1),normdF(i+1));
       break;
    end
end
toc;

% Plot of final segmentation
if ndims(I) == 2
    seg = figure;
    hold on;
    imagesc(I);
    colorbar();
    colormap('gray');
    contour(reshape(phi,M,N), [0 0], 'Color', [1 0 0],'LineWidth',3);
    axis tight;
    set(gca,'YDir','Reverse');
    hold off;
    print(seg,'-dpsc',strcat('I',no,'seg.eps'));
end

% Convergence plot of F
Fplot = figure;
plot(1:sum(F>0),F(F>0));
set(gca,'xtick',0:1:maxit);
title('Convergence plot','interpreter','latex','FontSize',15);
xlabel('Iteration no.','interpreter','latex','FontSize',15);
ylabel('Value of $F$','interpreter','latex','FontSize',15);
print(Fplot,'-dpsc',strcat('I',no,'con.eps'));

% Convergence plot of phi
phiplot = figure;
plot(1:sum(normdF>0),gamma*normdF(normdF>0));
set(gca,'xtick',0:1:maxit);
title('Convergence plot','interpreter','latex','FontSize',15);
xlabel('Iteration no.','interpreter','latex','FontSize',15);
ylabel('Difference in $\varphi$','interpreter','latex','FontSize',15);
print(phiplot,'-dpsc',strcat('I',no,'con2.eps'));

% For syntetic images: difference to reference image
dev = sum(abs((phi>=0)-Iref(:)))/(M*N);
fprintf('Deviation in pct. %f\n',dev);