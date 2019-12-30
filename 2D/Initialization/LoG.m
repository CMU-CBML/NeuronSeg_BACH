function resultat = LoG(x,y,sigma)
% function resultat = LoG(x,y,sigma)
% Laplacien d'une fonction gaussienne
% 
% Ce Laplacien est normalisé pour avoir un 
% scale space normalisé.
%
% Inputs :
%    x, y : points où l'on calcule la fonction
%    sigma : largeur de la gaussienne
% 
% Output :
%    resultat : valeur de la fonction en ce point

% Fait par JB Fiot pour l'assignement 1 du cours 
% de Reconnaissance d'objets et vision artificielle

% Date : Oct. 2008
    resultat = exp(-(x^2+y^2)/(2*sigma^2))*(x^2+y^2-2*sigma^2);
