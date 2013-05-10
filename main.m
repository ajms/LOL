%clear all;
close all;
%{
[FileName, PathName, FilterIndex] = uigetfile('*', 'Select image');
filename = strcat(PathName, FileName);
[pathstr, name, ext] = fileparts(filename);
if strcmpi(ext, '.mat')
    S = load(filename);
    I = S.I/255;
else
    RGB = imread(filename);
    I = double(rgb2gray(RGB))/255;
end
%}
[M, N] = size(I);

% parameters
plt = 1;
gamma = 100;
lambda1 = 1;
lambda2 = 1;
mu = 0;
nu = 0;
sigma = 2;

Ismooth = real(ifftn(scalen(fftn(I),[sigma,sigma],[0,0])));

phi = -ones(M,N);
[X Y] = meshgrid(1:M);
phip = (X-floor(M/2)).^2 + (Y-floor(N/2)).^2; 
phi(phip <= 50^2) = 1;
phi = init(phi);

a = figure
if plt == 1
    subplot(2,2,1);
    hold on;
    imagesc(I);
    contour(phi, [0 0],'Color',[1 0 0],'LineWidth',3);
    title('Initial contour and smooth image');
    axis tight;
    colormap('gray');
    colorbar();
    hold off;
end
%options.Method = 'lbfgs';
phi = phi/(max(phi(:))-min(phi(:)));
F = zeros(500,1);
diffF = zeros(500,1);
normdF = zeros(500,1);
dF = 0;
tic;
for i=1:100
    [F(i+1,1), dF] = lol(phi(:),Ismooth,lambda1,lambda2,mu,nu);
    %dF = minFunc(@lol,phi(:),options,Ismooth,lambda1,lambda2,mu,nu);
    x = phi(:)-gamma*dF;
    x = dF;
    if plt == 1
        subplot(2,2,2);
        surf(reshape(-gamma*dF,M,N));
        colorbar();
        subplot(2,2,3);
        surf(reshape(x,M,N));
        colorbar();
        title('Final levelset');
        subplot(2,2,4);
        hold on;
        imagesc(Ismooth);
        colormap('gray');
        contour(reshape(phi,M,N), [0 0], 'Color', [1 0 0], 'LineWidth',2);
        title('Final contour and original image');
        axis tight;
        colorbar();
        hold off;
        pause(0.01);
    end
    %phi = init(x);
    phi = x;
    diffF(i+1) = norm(F(i+1)-F(i));
    normdF(i+1) = norm(dF);
    if  diffF(i+1) < 1/gamma &&  normdF(i+1) < 1/gamma
       fprintf('Done after %d iterations!\n\n',i-1);
       break;
    end
end
toc;
if plt == 0
    hold on;
    imagesc(Ismooth);
    colorbar();
    colormap('gray');
    contour(reshape(phi,M,N), [0 0], 'Color', [1 0 0],'LineWidth',3);
    axis tight;
    hold off;
    pause(5);
    surf(reshape(x,M,N));
    colorbar();
end
 