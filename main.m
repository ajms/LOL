clear all;
close all;
%{
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
%}
RGB = imread('test_images/test2.jpg');
u0 = double(rgb2gray(RGB))/255;
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
    phip = (X-50).^2 + (Y-50).^2; 
    phi(phip <= 10^2) = 1;
elseif circle == 0
    phip = zeros(M,N);
    for i=umin:0.005:umean
        phip = phip + PW(i,0.005,usmooth);
    end
    phi(phip >= 0.5) = 1;
end    
phi = init(phi);

options.Method = 'csd';
x = minFunc(@lol,phi(:),options);
hold on;
imagesc(usmooth);
contour(phi, [0 0]);
contour(reshape(x,M,N), [0 0]);
hold off;