clear all;
close all;

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
%[X Y] = meshgrid(1:M);
%phip = (X-floor(M/2)).^2 + (Y-floor(N/2)).^2; 
%phi(phip <= 20^2) = 1;
phi(30:60,30:60) = 1;
phi = init(phi);

%options.Method = 'sd';
%figure
if plt == 1
    subplot(2,2,1);
    hold on;
    imagesc(Ismooth);
    [G H] = contour(phi, [0 0],'k');
    set(H,'LineWidth',3);
    title('Initial contour and smooth image');
    axis tight;
    colorbar();
    hold off;
end

phi = phi/(max(phi(:))-min(phi(:)));
F = ones(M*N,1);
tic;
for i=1:100
    tempF = F;
    %x = minFunc(@lol,phi(:),options,Ismooth,lambda1,lambda2,nu);
    [F, dF] = lol(phi(:),Ismooth,lambda1,lambda2,mu,nu);
    x = phi(:)-gamma*dF;
    if plt == 1
        subplot(2,2,2);
        surf(reshape(-gamma*dF,M,N));
        subplot(2,2,3);
        surf(reshape(x,M,N));
        colorbar();
        pause(0.01);
        title('Final levelset');
        subplot(2,2,4);
        hold on;
        imagesc(Ismooth);
        [G H] = contour(reshape(phi,M,N), [0 0]);
        set(H,'LineWidth',3);
        title('Final contour and original image');
        axis tight;
        colorbar();
        hold off;
    end
    phi = x;
    diffF = norm(F-tempF);
    NdF = norm(dF);
    fprintf('Iteration: %d, difference: %f, F=%f\n ', i, NdF, F);
    if  diffF < 0.01 &&  NdF < 0.01
       fprintf('Done after %d iterations!\n\n',i-1);
       break;
    end
end
toc;
if plt == 0
    hold on;
    imagesc(Ismooth);
    [H I] = contour(reshape(phi,M,N), [0 0]);
    axis tight;
    hold off;
    pause(5);
    surf(reshape(x,M,N));
end
 