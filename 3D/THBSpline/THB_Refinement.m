function [Elem_final,Basis_final,P_final,ActiveNodes] = THB_Refinement(Elem,Basis,Cpt,knotvectorU, knotvectorV,knotvectorW,ActiveBasis,CellGrad,meanGrad,param,multilev,ActiveNodes)

%% This function performs THBSpline refinement

%--Input Variables:
%Elem: Data struct variable storing the control grid element information
%Basis: Basis data struct variable storing the B-spline basis function information
%Cpt: THB-spline control points
%knotvectorU: knot vector in u direction
%knotvectorV: knot vector in v direction
%knotvectorW: knot vector in w direction
%ActiveBasis: Basis function to be refined
%CellGrad: Refinement criterion based on image gradient for the active
%elements
%meanGrad: mean image gradient value

%param: struct variable storing the parameters for the neuron
%segmentation consisting of the following fields
%rho: threshold parameter for the refinement of B-splines
%orderGauss: Gaussian quadrature order
%maxlevel: number of refinement levels
%lambda1, lambda2, lambda3: regularization parameters
%maxiteration: number of maximum iterations for each refinement level
%timestep: time step value for each refinement level
%nelemx, nelemy: number of elements in each direction at the first
%refinement level
%pU, pV: B-spline degree
%sigma: G_sigma: the key parameter which needs to be tuned properly for local image fitting
%sigma_phi:%G_sigma_phi: the variance of regularized Gaussian kernel
%epsilon:width of level set function

%multilev: current refinement level
%ActiveNodes: Nodes of active elements


%--Output Variables:

%Elem_final: Final data struct variable storing the control grid element
%information after refinement
%Basis_final: Final basis data struct variable storing the B-spline basis function information
%after refinement
%P_final: THB-spline control points after refinement
%ActiveNodes: Nodes of active elements after refinement


%convert struct aray to cell array
% the reason is refinement function is faster when it works on cell arrays
% than struct data structure

Em = struct2cell(Elem)';
Dm = struct2cell(Basis)';
Pm = struct2cell(Cpt)';
bf_ct = size(ActiveBasis,1);
bf = ActiveBasis;

%degree of B-splines
pU = param.pU;
pV = param.pV;
pW = param.pW;

%loop over the active B-splines
for j =1:bf_ct
    
    %B-spline index, B-spline level
    bbc = bf(j,1);
    bf_lev = bf(j,2);
    
    %refinement parameter
    rho = param.rho(multilev);
    
    %loop over the support cells of the active splines
    supp_cells = Basis(bf_lev).suppCell(bbc,:);
    
    %compute the average image difference value over the supporting
    %elements of the active B-spline
    grad = 0;
    supp_ct = 0;
    for i =1:size(supp_cells,2)
        if(supp_cells(1,i)~=0)
            supp_ct = supp_ct + 1;
            ac_ind = Elem(bf_lev).actE(supp_cells(1,i),1);
            grad  = grad + CellGrad(ac_ind,1);
        end
    end
    grad = grad/supp_ct;
    
    %Refinement to create next level
    %if the average value of Ig over the support elements is greater than
    %mean image difference times rho, REFINE spline
    if(grad>=rho*meanGrad)
        [Dm,Em,Pm,ActiveNodes] =  Refine3D_truncation(bbc,bf_lev,Dm,Em,Pm,knotvectorU,knotvectorV,knotvectorW,pU,pV,pW,ActiveNodes);
    end
end


%convert cell array to struct array
fieldElem = {'knot_ind','flag_active','IEN','chdElem','cell_centre','node','parElem','actE','nodes'};
Elem_final = cell2struct(Em,fieldElem,2);
fieldBasis = {'basis_ind','flag_active','chdBasis','coeffBasis','suppCell','flag_trunc','flag_ref','flag_bdy','actB'};
Basis_final = cell2struct(Dm,fieldBasis,2);
P_final = cell2struct(Pm,'pts',2);
end