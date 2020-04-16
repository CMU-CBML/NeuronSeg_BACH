function resultat = LoG(x,y,sigma)
% function resultat = LoG(x,y,sigma)
% Laplacien d'une fonction gaussienne
%
% This function computes the LoG filter kernel entries
%
% Inputs :
%    size : size of LoG kernel
%    scale : radius of disk shaped filter
%
% Output :
%    resultat : matrix entries for the LoG kernel

% Fait par JB Fiot pour l'assignement 1 du cours
% de Reconnaissance d'objets et vision artificielle

% Date : Oct. 2008
resultat = exp(-(x^2+y^2)/(2*sigma^2))*(x^2+y^2-2*sigma^2);

end
