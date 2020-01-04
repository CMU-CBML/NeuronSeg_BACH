function [phi,VF] = MultipleResolution3D_heliximg(F,phi,param,Img1)

[X,Y,Z] = meshgrid(0.5:size(F,1)-0.5,0.5:size(F,2)-0.5,0.5:size(F,3)-0.5);

img_pU = 3;
img_pV = 3;
img_pW = 3;
orderGauss = param.orderGauss;

[FDX,FDY,FDZ] = gradient(F);
[DphiX,DphiY,DphiZ] = gradient(phi);

F_orig = F;
phi_orig = phi;

phi00 = phi;
F00 = F;

%Basis value computed by substituting corresponding midpoint
%Bspline kernel
[Bsplinekernel] = BsplineKernel3D;

% Display intial config
iterct = 0;
xlen = 4;

nx = size(F,1);
ny = size(F,2);
nz = size(F,3);

%% Multilevel framework
VXLF=X(:)';
VYLF=Y(:)';
VZLF=Z(:)';

coef_xf = img2coef3D_mex(VXLF,nx,ny,nz);   %coefficients for position x
coef_yf = img2coef3D_mex(VYLF,nx,ny,nz);   %coefficients for position y
coef_zf = img2coef3D_mex(VZLF,nx,ny,nz);   %coefficients for position z

coef_xf_fun = reshape(coef_xf, 1, (nx+img_pU)*(ny+img_pV)*(nz+img_pW));
coef_yf_fun = reshape(coef_yf, 1, (nx+img_pU)*(ny+img_pV)*(nz+img_pW));
coef_zf_fun = reshape(coef_zf, 1, (nx+img_pU)*(ny+img_pV)*(nz+img_pW));

% Start multilevel
VXf = X;
VYf = Y;
VZf = Z;

ori_size_x=size(F,2);
ori_size_y=size(F,1);
ori_size_z=size(F,3);

disp('coefficients of initial vectors');
toc

for level = param.maxlevel:-1:1
    tic
    disp(['Registration level:' num2str(param.maxlevel-level+1)]);
    
    %% downsample image
    scale = 2^-(level-1);  % image scale
    
    new_size_x=round(ori_size_x*scale);
    new_size_y=round(ori_size_y*scale);
    new_size_z=round(ori_size_z*scale);
    npx = new_size_x; %size of scaled image
    npy = new_size_y;
    npz = new_size_z;
    
    F = resize3Dmatrix(new_size_x,new_size_y,new_size_z,F);
    vecF = F(:)';
    coef_F= img2coef3D_mex(vecF,npx,npy,npz);
    CIF = coef_F(:)';
    coef_matF = reshape(CIF, npy+img_pU, npx+img_pV, npz+img_pW); % this is an unnecessary step, check later!!!!

    phi = resize3Dmatrix(new_size_x,new_size_y,new_size_z,phi); % does not matter if phi00 or phi, because this only counts
    % for the first level which is anyway phi00.
    vecP = phi(:)';
    coef_P = img2coef3D_mex(vecP,npx,npy,npz);
    CIP = coef_P(:)';
    coef_matP = reshape(CIP,npy+img_pU, npx+img_pV, npz+img_pW);

    F_orig = resize3Dmatrix(new_size_x,new_size_y,new_size_z,F00);
    vecForig = F_orig(:)';
    coef_Forig= img2coef3D_mex(vecForig,npx,npy,npz);
    CIForig = coef_Forig(:)';
    coef_matForig = reshape(CIForig,npy+img_pU, npx+img_pV, npz+img_pW);
    
    phi_orig = resize3Dmatrix(new_size_x,new_size_y,new_size_z,phi00);
    vecPorig = phi_orig(:)';
    coef_Porig= img2coef3D_mex(vecPorig,npx,npy,npz);
    CIPorig = coef_Porig(:)';
    coef_matPorig = reshape(CIPorig,npy+img_pU, npx+img_pV, npz+img_pW);

    [FDX,FDY,FDZ] = gradient(F);
    [DphiX,DphiY,DphiZ] = gradient(phi);
    
    disp('coefficient of scaled images');
    toc
    tic

    if(level==param.maxlevel)
        
        nx = size(F,2);
        ny = size(F,1);
        nz = size(F,3);
        
        [Vx_ident,Vy_ident,Vz_ident] = meshgrid((0.5:nx-0.5),(0.5:ny-0.5),(0.5:nz-0.5));
        [X_vtk,Y_vtk,Z_vtk] = meshgrid(1:nx,1:ny,1:nz);
        
        Vxlf = Vx_ident;
        Vylf = Vy_ident;
        Vzlf = Vz_ident;

        Vxlf_inv = Vx_ident;
        Vylf_inv = Vy_ident;
        Vzlf_inv = Vz_ident;

    else
        
        nx = size(F,2);
        ny = size(F,1);
        nz = size(F,3);
        
        [Vx_ident, Vy_ident,Vz_ident] = meshgrid((0.5:nx-0.5),(0.5:ny-0.5),(0.5:nz-0.5));
        [X_vtk,Y_vtk,Z_vtk] = meshgrid(1:nx,1:ny,1:nz);
        
        Vxlf = resize3Dmatrix(new_size_x,new_size_y,new_size_z,Vxf*scale);
        Vylf = resize3Dmatrix(new_size_x,new_size_y,new_size_z,Vyf*scale);
        Vzlf = resize3Dmatrix(new_size_x,new_size_y,new_size_z,Vzf*scale);

        Vxlf_inv = resize3Dmatrix(new_size_x,new_size_y,new_size_z,Vxf_inv*scale);
        Vylf_inv = resize3Dmatrix(new_size_x,new_size_y,new_size_z,Vyf_inv*scale);
        Vzlf_inv = resize3Dmatrix(new_size_x,new_size_y,new_size_z,Vzf_inv*scale);

        Vxlf = Vxlf+Vx_ident;
        Vylf = -Vylf+Vy_ident;
        Vzlf = Vzlf+Vz_ident;

        Vxlf_inv = Vxlf_inv+Vx_ident;
        Vylf_inv = -Vylf_inv+Vy_ident;
        Vzlf_inv = Vzlf_inv + Vz_ident;

        II0F_ori = F_orig(:)';
        II0phi_ori = phi_orig(:)';
        
        %% coef is in the vector form, see cdmffd code

        coef_IF = img2coef3D_mex(II0F_ori, nx,ny, nz);
        coef_Iphi = img2coef3D_mex(II0phi_ori,nx,ny,nz);

        CI0phi = coef_Iphi(:)';
        
        VXLf = Vxlf(:)';
        VYLf = Vylf(:)';
        VZLf = Vzlf(:)';

        phi = BsplineComposeImage3D_mex(nx,ny,nz,VXLf, VYLf, VZLf, CI0phi);
        vecP = phi(:)';
        coef_P = img2coef3D_mex(vecP, nx, ny,nz);
        CIP = coef_P(:)';
        coef_matP = reshape(CIP,npy+img_pU, npx+img_pV, npz+img_pW);

        [FDX,FDY,FDZ] = gradient(F);
        [DphiX,DphiY,DphiZ] = gradient(phi);

    end
    
    %removed the grid plotting functions
    VXLf = Vxlf(:)';
    VYLf = Vylf(:)';
    VZLf = Vzlf(:)';

    coef_xf = img2coef3D_mex(VXLf,nx,ny,nz);
    coef_yf = img2coef3D_mex(VYLf,nx,ny,nz);
    coef_zf = img2coef3D_mex(VZLf,nx,ny,nz);

    VXLf_inv = Vxlf_inv(:)';
    VYLf_inv = Vylf_inv(:)';
    VZLf_inv = Vzlf_inv(:)';

    coef_xf_inv = img2coef3D_mex(VXLf_inv,nx,ny,nz);
    coef_yf_inv = img2coef3D_mex(VYLf_inv,nx,ny,nz);
    coef_zf_inv = img2coef3D_mex(VZLf_inv,nx,ny,nz);

    disp('set initial vectors for the level');
    toc
    tic
    
    %% Construct B-spline grid
    maxlev = param.maxlevel-level+1;
    
    setBsplineGrid3D
    
    ActiveNodes = [];
    Node = [];
    
    for levels = 1:maxlev,
        [nx1,ny1,nz1] = meshgrid(uknotvectorV{levels,1},uknotvectorU{levels,1},uknotvectorW{levels,1});
        Node = [Node;ny1(:),nx1(:),nz1(:)];
    end
    
    disp('set bspline grid');
    toc
    tic
    %Loop over each refinement level
    disp('Loop over each refinement level...');
    for multil = 0:1:maxlev-1
        
        iteration = 0;
        
        fprintf('Refinement at level %i...\n',multil+1);
        if(multil>0)
            [Em,Dm,Pm,ActiveNodes] = THB_Refinement(Em,Dm,Pm,knotvectorU, knotvectorV,knotvectorW,bf,CellGrad,meanGrad,param,multil,ActiveNodes);
        end
        
        disp('Collecting active elements, control points and basis functions...');
        [ac, bf, ACP, RHS,Em,Dm,ActiveNodes] = storeActiveElem(Em,Dm,Pm,multil,ActiveNodes);
        
        ac_ct = size(ac,1);
        bf_ct = size(bf,1);
        
        ActiveNodes = unique(ActiveNodes);
        ACCf = ACP;
        ACPf = ACP;
        ACPold = ACP;
        
        cell_co = zeros(ac_ct,3);
        for i = 1:ac_ct
            cell_id = ac(i,1);
            cell_le = ac(i,2);
            cell_co(i,:) = Em(cell_le).cell_centre(cell_id,:);
        end
        
        Idiff = sqrt((FDX).^2 + (FDY).^2 + (FDZ).^2);
        CellGrad = interp3(Vx_ident,Vy_ident,Vz_ident,Idiff,cell_co(:,2),cell_co(:,1),cell_co(:,3));
        meanGrad = mean2(Idiff);
    end
    
    Pmold = Pm;
    
    [Jm, Coeff] = computeNonZeroSplines(ac, param, Em, Dm, multil);
    
    disp('Computing the basis functions at pixel coordinates...');
    numPixels = int64(prod(sizeImage));
    pix =  [Vy_ident(:),Vx_ident(:),Vz_ident(:)];
    [Pixel, Pix2] = storePixelPhi(numPixels, multil,pix, knotvectorU, knotvectorV, knotvectorW, Em, Coeff, param);
    for p_ind = 1:numPixels
        Pixel(p_ind).phi = Pix2{p_ind};
    end
    clear Pix2
    
    disp('Computing the basis functions at gaussian points...');
    %compute the gaussian points and weights of the given gauss order
    [Gu,Wu] = ggquad(param.orderGauss);
    [Gv,Wv] = ggquad(param.orderGauss);
    [Gw,Ww] = ggquad(param.orderGauss);
    
    [PHI,PHIU,PHIV,PHIW,BIGX,BIGY,BIGZ,H] = GaussPhi(ac,Em,knotvectorU,knotvectorV,knotvectorW,Coeff,param,maxlev);
    % interpolate the intensity values of the target image at the gauss
    % points stored in BIGX, BIGY, BIGZ
    BIGXF = BIGX;
    BIGYF = BIGY;
    BIGZF = BIGZ;
    
    %cI0fd = interp3(X_vtk,Y_vtk,Z_vtk,F_orig,BIGY,BIGX,BIGZ,'*linear',min(M(:)));
    [cI0f, tempu1, tempu2, tempu3] =  BsplineComposeImage3D_single_mex(BIGY, BIGX, BIGZ, coef_matForig, size(BIGX,1), size(BIGX,2), size(BIGX,3));

    Node(:,4) = 0;
    for i=1:size(ActiveNodes,1)
        Node(ActiveNodes(i,1),4) = i;
    end
    
    PlotGrid
    
    disp('refinement done')
    toc
    tic
   
    %THBcheck
    
    %% Update the iteration loop here
    timestep = param.timestep(param.maxlevel-level+1,1);
    
    RHS_init = RHS;
    
    figure
    
    PHI1 = cell2struct(PHI,'mat',2);
    PHIU1 = cell2struct(PHIU,'mat',2);
    PHIV1 = cell2struct(PHIV,'mat',2);
    PHIW1 = cell2struct(PHIW,'mat',2);
    
    disp('set iteration parameters');
    toc;

    
    VXf_new = Vx_ident;
    VYf_new = Vy_ident;
    VZf_new = Vz_ident;
    ACCf = ACPold;
    
    iterationloop_heliximg
    
    tic
    Vxf = resize3Dmatrix(size(F00,1), size(F00,2), size(F00,3),Vxf/scale);
    Vyf = resize3Dmatrix(size(F00,1), size(F00,2), size(F00,3),Vyf/scale);
    Vzf = resize3Dmatrix(size(F00,1), size(F00,2), size(F00,3),Vzf/scale);

    Vxf_inv = resize3Dmatrix(size(F00,1), size(F00,2), size(F00,3),Vxf_inv/scale);
    Vyf_inv = resize3Dmatrix(size(F00,1), size(F00,2), size(F00,3),Vyf_inv/scale);
    Vzf_inv = resize3Dmatrix(size(F00,1), size(F00,2), size(F00,3),Vzf_inv/scale);

    F = resize3Dmatrix(size(F00,1), size(F00,2), size(F00,3),F);
    phi = resize3Dmatrix(size(F00,1), size(F00,2), size(F00,3),phi);
    
    disp('scale up images')
    toc
    
end
VF{1,1} =  VXf_new;
VF{1,2} =  VYf_new;
VF{1,3} =  VZf_new;
end
