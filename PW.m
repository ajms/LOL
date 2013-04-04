function [ isophote ] = PW( i,beta,usmooth )
%PW Summary of this function goes here
%   Detailed explanation goes here
    ui = usmooth-i;
    isophote = exp(-ui.^2/(2*beta^2));
end

