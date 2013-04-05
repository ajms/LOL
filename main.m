clear all;
close all;

[FileName, PathName, FilterIndex] = uigetfile('*', 'Select image');
filename = strcat(PathName, FileName);
[pathstr, name, ext] = fileparts(filename);
if strcmpi(ext, '.mat')
    S = load(filename);
    u0 = S.I/255;
else
    RGB = imread(filename);
    u0 = double(rgb2gray(RGB))/255;
end
%{
RGB = imread('test_images/test2.jpg');
u0 = double(rgb2gray(RGB))/255;
%}
[M, N] = size(u0);

% parameters
circle = 1;
h = 1.0;
dt = 0.1;
di = 0.5;
lambda1 = 1;
lambda2 = 1;
mu = 0;%0.01*255^2;
nu = 200;
sigma = 0.5;
beta = 0.01;
doreinit = 0;

usmooth = real(ifftn(scalen(fftn(u0),[2,2],[0,0])));
umin = floor(min(usmooth(:)));
umax = ceil(max(usmooth(:)));
uinc = (umax-umin)/100;
ulev = umin:uinc:umax;
umean = mean(usmooth(:));

phi = -ones(M,N);
if circle == 1
    [X Y] = meshgrid(1:M);
    phip = (X-150).^2 + (Y-150).^2; 
    phi(phip <= 50^2) = 1;
elseif circle == 0
    phip = zeros(M,N);
    for i=umin:0.005:umean
        phip = phip + PW(i,0.005,usmooth);
    end
    phi(phip >= 0.5) = 1;
end    

phi = init(phi);

options.Method = 'lbfgs';
figure
pause(2);
subplot(2,2,1);
hold on;
imagesc(usmooth);
[G H] = contour(phi, [0 0],'k');
set(H,'LineWidth',3);
title('Initial contour and smooth image');
axis tight;
colorbar();
hold off;

subplot(2,2,2);
surf(phi);
title('Initial levelset');
colorbar();

for i=1:20
    fprintf('Iteration: %d\n',i);
    if i<10
        phi = phi/50;
    end
    x = minFunc(@lol,phi(:),options);
    subplot(2,2,3);
    surf(reshape(x,M,N));
    colorbar();
    phi = -ones(M,N);
    phi(x>=0) = 1;
    phi = init(phi);
    pause(0.1);
end
title('Final levelset');
subplot(2,2,4);
hold on;
imagesc(u0);
[G H] = contour(phi, [0 0]);
set(H,'LineWidth',3);
title('Final contour and original image');
axis tight;
colorbar();
hold off;
