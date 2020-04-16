function BBV = bcondition1(BBvector,D, bff, nbu)

%% This function applies Dirichlet zero deformation boundary condition on the B-spline grid

%--Input Variable:
%BBvector: the RHS vector for updating the control points
%D: THB-spline data structure element
%bff: B-spline indices of active splines
%nbu: number of B-splines

%--Output Variable:
%BBV: output RHS vector for updating the control points after applying
%Dirichlet BC

nlev = size(BBvector,1);

for i=1:nlev
    BB = D{bff(i,3),1};
    bbc = bff(i,1:2);
    bbin = nbu(bff(i,3),1)*(bbc(1,2)-1) + bbc(1,1);
    if(BB{bbin,9}==1)
        BBvector(i,1) = 0;
        BBvector(i,2) = 0;
    end
end
BBV = BBvector;
end
