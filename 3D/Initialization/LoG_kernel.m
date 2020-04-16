function resultat = LoG_kernel(size,scale)
% function resultat = LoG_kernel(size,scale)
%
% This function computes the LoG filter kernel
%
% Inputs :
%    size : size of LoG kernel
%    scale : radius of disk shaped filter
%
% Output :
%    resultat : LoG kernel
%
% Dépendance : utilise la fonction LoG.m

% Fait par JB Fiot pour l'assignement 1 du cours
% de Reconnaissance d'objets et vision artificielle

% Date : Oct. 2008

resultat =zeros(size,size);
for i=1:size
    for j=1:size
        resultat(i,j)=LoG((size+1)/2-i,(size+1)/2-j,scale);
    end
end

end