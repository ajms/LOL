M = 100;
N = 100;

% Circle + square: small=1, dt=0.5, mu=0.1*255^2, 51 iterations,
% Deviation: 9
[X Y] = meshgrid(1:N,1:M);
I0 = ((X-floor(N/2)).^2 + (Y-floor(M/2)).^2)<=20^2;
I0(floor(M/2)+5:floor(M/2)+30,floor(N/2)+5:floor(N/2)+30) = 1;

I1 = I0 + random('unif',-0.2,0.2,M,N);
I1 = (I1 + abs(min(I1(:))))/(max(I1(:))-min(I1(:)));

% random disjoint squares: small=1, dt=0.5, mu=0.1*255^2, 62
% iterations, Deviation = 27
I2 = zeros(M,N);
randx = sort(randi([1,M],[4,1]));
randy = sort(randi([1,N],[4,1]));
I2(randx(2):randx(3),randy(1):randy(2))=1;
I2(randx(1):randx(4),randy(3):randy(4))=1;

I3 = I2 + random('unif',-0.3,0.3,M,N);
I3 = (I3 + abs(min(I3(:))))/(max(I3(:))-min(I3(:)));

% circles: small=1, dt=0.1, mu = 0.4*255^2, 93 iterations
[X Y] = meshgrid(1:N,1:M);
I4 = (((X-floor(N/2)-10).^2 + (Y-floor(M/2)-10).^2)<=4^2) ... 
     + (((X-floor(N/2)).^2 + (Y-floor(M/2)).^2)<=4^2) ...
     + (((X-floor(N/2)-10).^2 + (Y-floor(M/2)).^2)<=4^2) ...
     + (((X-floor(N/2)).^2 + (Y-floor(M/2)-10).^2)<=4^2) ...
     + (((X-floor(N/2)-5).^2 + (Y-floor(M/2)+10).^2)<=4^2);

I5 = I4 + random('unif',-0.5,0.5,M,N);
I5 = (I5 + abs(min(I5(:))))/(max(I5(:))-min(I5(:)));

% import image
S0 = load('/home/albert/Dropbox/Uni/Bachelorprojekt/LOL/medical_images/images8.mat');
I6 = S0.I;

S1 = load('/home/albert/Dropbox/Uni/Bachelorprojekt/LOL/medical_images/brain.mat');
I7 = S1.img2;


%{
[FileName, PathName, FilterIndex] = uigetfile('*', 'Select image');
filename = strcat(PathName, FileName);
[pathstr, name, ext] = fileparts(filename);
if strcmpi(ext, '.mat')
    S = load(filename);
    I = S.I;
else
    RGB = imread(filename);
    I = double(rgb2gray(RGB));
end

I = real(ifftn(scalen(fftn(I),[10,10],[0,0])));
%}
