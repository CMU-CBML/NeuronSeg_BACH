function [phi,VF] = MultipleResolution3D_heliximg(F,phi,param,Img1)

%-- In this function, the neuron segmentation for multiple refinement
%levels is carried out.

%--Input Variables
%F: input image to be segmented
%phi: level set function describing the initial guess for the segmented image


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

% Img1: The original image without noise to check accuracy

%--Output Variables

%phi: level set function capturing the final segmented image result
%VF: deformation field of the evolution of the level set function

img_pU = 3;
img_pV = 3;
img_pW = 3;
orderGauss = param.orderGauss;

phi00 = phi;
F00 = F;

%Basis value computed by substituting corresponding midpoint
%Bspline kernel
[Bsplinekernel] = BsplineKernel3D;

% Display intial config
iterct = 0;

%% Multilevel framework


% Start multilevel

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
    
    phi = resize3Dmatrix(new_size_x,new_size_y,new_size_z,phi); % does not matter if phi00 or phi, because this only counts
    % for the first level which is anyway phi00.
    vecP = phi(:)';
    coef_P = img2coef3D(vecP,npx,npy,npz);
    CIP = coef_P(:)';
    coef_matP = reshape(CIP,npy+img_pU, npx+img_pV, npz+img_pW);
    
    F_orig = resize3Dmatrix(new_size_x,new_size_y,new_size_z,F00);
    vecForig = F_orig(:)';
    coef_Forig= img2coef3D(vecForig,npx,npy,npz);
    CIForig = coef_Forig(:)';
    coef_matForig = reshape(CIForig,npy+img_pU, npx+img_pV, npz+img_pW);
    
    phi_orig = resize3Dmatrix(new_size_x,new_size_y,new_size_z,phi00);
    
    [FDX,FDY,FDZ] = gradient(F);
    
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
        
        II0phi_ori = phi_orig(:)';
        
        %% coef is in the vector form, see cdmffd code
        
        coef_Iphi = img2coef3D(II0phi_ori,nx,ny,nz);
        
        CI0phi = coef_Iphi(:)';
        
        VXLf = Vxlf(:)';
        VYLf = Vylf(:)';
        VZLf = Vzlf(:)';
        
        phi = BsplineComposeImage3D(nx,ny,nz,VXLf, VYLf, VZLf, CI0phi);
        vecP = phi(:)';
        coef_P = img2coef3D(vecP, nx, ny,nz);
        CIP = coef_P(:)';
        coef_matP = reshape(CIP,npy+img_pU, npx+img_pV, npz+img_pW);
        
        [FDX,FDY,FDZ] = gradient(F);
        
        
    end
    
    %removed the grid plotting functions
    VXLf = Vxlf(:)';
    VYLf = Vylf(:)';
    VZLf = Vzlf(:)';
    
    coef_xf = img2coef3D(VXLf,nx,ny,nz);
    coef_yf = img2coef3D(VYLf,nx,ny,nz);
    coef_zf = img2coef3D(VZLf,nx,ny,nz);
    
    VXLf_inv = Vxlf_inv(:)';
    VYLf_inv = Vylf_inv(:)';
    VZLf_inv = Vzlf_inv(:)';
    
    coef_xf_inv = img2coef3D(VXLf_inv,nx,ny,nz);
    coef_yf_inv = img2coef3D(VYLf_inv,nx,ny,nz);
    coef_zf_inv = img2coef3D(VZLf_inv,nx,ny,nz);
    
    disp('set initial vectors for the level');
    toc
    tic
    
    %% Construct B-spline grid
    maxlev = param.maxlevel-level+1;
    
    setBsplineGrid3D
    
    ActiveNodes = [];
    Node = [];
    
    for levels = 1:maxlev
        [nx1,ny1,nz1] = meshgrid(uknotvectorV{levels,1},uknotvectorU{levels,1},uknotvectorW{levels,1});
        Node = [Node;ny1(:),nx1(:),nz1(:)];
    end
    
    disp('set bspline grid');
    toc
    tic
    %Loop over each refinement level
    disp('Loop over each refinement level...');
    for multil = 0:1:maxlev-1
        
        
        fprintf('Refinement at level %i...\n',multil+1);
        if(multil>0)
            [Em,Dm,Pm,ActiveNodes] = THB_Refinement(Em,Dm,Pm,knotvectorU, knotvectorV,knotvectorW,bf,CellGrad,meanGrad,param,multil,ActiveNodes);
        end
        
        disp('Collecting active elements, control points and basis functions...');
        [ac, bf, ACP, RHS,Em,Dm,ActiveNodes] = storeActiveElem(Em,Dm,Pm,multil,ActiveNodes);
        
        ac_ct = size(ac,1);
        ActiveNodes = unique(ActiveNodes);
        
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
    [~,Wu] = ggquad(param.orderGauss);
    [~,Wv] = ggquad(param.orderGauss);
    [~,Ww] = ggquad(param.orderGauss);
    
    [PHI,PHIU,PHIV,PHIW,BIGX,BIGY,BIGZ,H] = GaussPhi(ac,Em,knotvectorU,knotvectorV,knotvectorW,Coeff,param,maxlev);
    % interpolate the intensity values of the target image at the gauss
    % points stored in BIGX, BIGY, BIGZ
    
    
    %cI0fd = interp3(X_vtk,Y_vtk,Z_vtk,F_orig,BIGY,BIGX,BIGZ,'*linear',min(M(:)));
    [cI0f, ~, ~, ~] =  BsplineComposeImage3D_single(BIGY, BIGX, BIGZ, coef_matForig, size(BIGX,1), size(BIGX,2), size(BIGX,3));
    
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
    
    K = fspecial3('gaussian',2*round(2*param.sigma)+1,param.sigma);
    
    for iteration = 1:param.maxiteration+1
        
        tic
        Pm = Pmold;
        RHS = RHS_init;
        Bvectf = RHS_init;
        
        D1E = Delta(phi,param.epsilon);
        vecD = D1E(:)';
        coef_D = img2coef3D(vecD,nx,ny,nz);
        CID = coef_D(:)';
        coef_matD = reshape(CID, nx+img_pU,ny+img_pV, nz+img_pW);
        
        Hphi = Heaviside(phi,param.epsilon);
        vecH = Hphi(:)';
        coef_H = img2coef3D(vecH,nx,ny,nz);
        CIH = coef_H(:)';
        coef_matH = reshape(CIH, nx+img_pU,ny+img_pV, nz+img_pW);
        
        [f1,f2] = Local_Avr(F_orig,Hphi,K);
        vec_f1 = f1(:)';
        coef_f1 = img2coef3D(vec_f1,nx,ny,nz);
        CI_f1 = coef_f1(:)';
        coef_matf1 = reshape(CI_f1, nx+img_pU,ny+img_pV, nz+img_pW);
        
        vec_f2 = f2(:)';
        coef_f2 = img2coef3D(vec_f2,nx,ny,nz);
        CI_f2 = coef_f2(:)';
        coef_matf2 = reshape(CI_f2, nx+img_pU,ny+img_pV, nz+img_pW);
        
        iterct = iterct +1;
        
        disp('iterct:');
        disp(iterct);
        toc
        tic
        [~,cDphiY, cDphiX,cDphiZ] =  BsplineComposeImage3D_single(BIGY, BIGX, BIGZ, coef_matP, size(BIGX,1), size(BIGX,2), size(BIGX,3));
        [dDelta, ~, ~, ~] =  BsplineComposeImage3D_single(BIGY, BIGX, BIGZ, coef_matD, size(BIGX,1), size(BIGX,2), size(BIGX,3));
        [dHphi, ~, ~, ~] =  BsplineComposeImage3D_single(BIGY, BIGX, BIGZ, coef_matH, size(BIGX,1), size(BIGX,2), size(BIGX,3));
        [df1, ~, ~, ~] =  BsplineComposeImage3D_single(BIGY, BIGX, BIGZ, coef_matf1, size(BIGX,1), size(BIGX,2), size(BIGX,3));
        [df2, ~, ~, ~] =  BsplineComposeImage3D_single(BIGY, BIGX, BIGZ, coef_matf2, size(BIGX,1), size(BIGX,2), size(BIGX,3));
        
        disp('cphi');
        toc
        
        tic
        Bseg = -dDelta.*(cI0f -df1.*dHphi -df2.*(1-dHphi)).*(df1-df2);
        Bsegx = cDphiY;
        Bsegy  = cDphiX;
        Bsegz = cDphiZ;
        
        disp('Bterm1');
        toc
        tic
        
        RHS = compute_Integ_Domain_fidelity_mex(Jm,Bseg, Bsegx,Bsegy,Bsegz,RHS,PHI1,Wu,Wv,Ww,H);
        
        toc
        tic
        RHS = bcondition3D(RHS);
        
        disp('integration done')
        toc
        tic
        
        ACCf(:,1:3) = ACCf(:,1:3) - timestep.*RHS(:,1:3);
        
        tic
        [~, ~, ~, BIGMUXf, BIGMUYf, BIGMUZf, BIGMVXf, BIGMVYf, BIGMVZf, BIGMWXf, BIGMWYf, BIGMWZf] =  computenewPoints(Jm, ACCf, PHI1, PHIU1, PHIV1, PHIW1, orderGauss);
        
        disp('BIGXX computed');
        toc
        
        Bvectf = compute_Integ_Domain_regularization(Jm,BIGMUXf, BIGMUYf, BIGMUZf, BIGMVXf, BIGMVYf, BIGMVZf, BIGMWXf, BIGMWYf, BIGMWZf, Bvectf,PHIU1, PHIV1, PHIW1,param.lambda_1,param.lambda_2,param.lambda_3,Wu,Wv,Ww,H);
        
        clear('BIGXXf','BIGYYf','BIGZZf','BIGMUXf','BIGMUYf','BIGMUZf','BIGMVXf','BIGMVYf','BIGMVZf','BIGMWXf','BIGMWYf','BIGMWZf');
        toc
        
        tic
        Bvectf = bcondition3D(Bvectf);
        disp('integration done')
        toc
        tic
        ACPf(:,1:3) = ACCf(:,1:3) - timestep.*Bvectf(:,1:3);
        tic
        
        [pxxf,pyyf,pzzf] = tripleIterLoop_mex(sizeImage, Pixel, Jm, ACPf);
        
        Vxf_temp = pxxf-Vx_ident;
        Vyf_temp = pyyf-Vy_ident;
        Vzf_temp = pzzf-Vz_ident;
        
        Vxf_old = Vxf_temp + Vx_ident;
        Vyf_old = Vyf_temp + Vy_ident;
        Vzf_old = Vzf_temp + Vz_ident;
        
        Vxf_old = reshape(Vxf_old,1,nx*ny*nz);
        Vyf_old = reshape(Vyf_old,1,nx*ny*nz);
        Vzf_old = reshape(Vzf_old,1,nx*ny*nz);
        
        coef_xf_fun = reshape(coef_xf, 1, (nx+img_pU)*(ny+img_pV)*(nz+img_pW));
        coef_yf_fun = reshape(coef_yf, 1, (nx+img_pU)*(ny+img_pV)*(nz+img_pW));
        coef_zf_fun = reshape(coef_zf, 1, (nx+img_pU)*(ny+img_pV)*(nz+img_pW));
        
        Vxf_old_inv = -Vxf_temp + Vx_ident;
        Vyf_old_inv = -Vyf_temp + Vy_ident;
        Vzf_old_inv = -Vzf_temp + Vz_ident;
        
        Vxf_old_inv = reshape(Vxf_old_inv,1,nx*ny*nz);
        Vyf_old_inv = reshape(Vyf_old_inv,1,nx*ny*nz);
        Vzf_old_inv = reshape(Vzf_old_inv,1,nx*ny*nz);
        
        coef_xf_inv_fun = reshape(coef_xf_inv, 1, (nx+img_pU)*(ny+img_pV)*(nz+img_pW));
        coef_yf_inv_fun = reshape(coef_yf_inv, 1, (nx+img_pU)*(ny+img_pV)*(nz+img_pW));
        coef_zf_inv_fun = reshape(coef_zf_inv, 1, (nx+img_pU)*(ny+img_pV)*(nz+img_pW));
        
        [VXf_new, VYf_new, VZf_new] = BsplineCompose3D_mex(nx, ny, nz, Vxf_old, Vyf_old, Vzf_old, coef_xf_fun, coef_yf_fun, coef_zf_fun);
        [VXf_new_inv, VYf_new_inv, VZf_new_inv] = BsplineCompose3D_mex( nx, ny, nz, Vxf_old_inv, Vyf_old_inv, Vzf_old_inv, coef_xf_inv_fun, coef_yf_inv_fun, coef_zf_inv_fun);
        
        temp_coef_xf = imfilter(VXf_new, Bsplinekernel); %Vx
        temp_coef_yf = imfilter(VYf_new, Bsplinekernel); %Vy
        temp_coef_zf = imfilter(VZf_new, Bsplinekernel);
        
        coef_xf(4:ny,4:nx, 4:nz) = temp_coef_xf(2:end-2,2:end-2,2:end-2);
        coef_yf(4:ny,4:nx, 4:nz) = temp_coef_yf(2:end-2,2:end-2,2:end-2);
        coef_zf(4:ny,4:nx, 4:nz) = temp_coef_zf(2:end-2,2:end-2,2:end-2);
        
        temp_coef_xf_inv = imfilter(VXf_new_inv, Bsplinekernel); %Vx
        temp_coef_yf_inv = imfilter(VYf_new_inv, Bsplinekernel); %Vy
        temp_coef_zf_inv = imfilter(VZf_new_inv, Bsplinekernel);
        
        coef_xf_inv(4:ny,4:nx, 4:nz) = temp_coef_xf_inv(2:end-2,2:end-2,2:end-2);
        coef_yf_inv(4:ny,4:nx, 4:nz) = temp_coef_yf_inv(2:end-2,2:end-2,2:end-2);
        coef_zf_inv(4:ny,4:nx, 4:nz) = temp_coef_zf_inv(2:end-2,2:end-2,2:end-2);
        
        Vxf = VXf_new - Vx_ident;
        Vyf = VYf_new - Vy_ident;
        Vzf = VZf_new - Vz_ident;
        
        Vxf_inv = VXf_new_inv - Vx_ident;
        Vyf_inv = VYf_new_inv - Vy_ident;
        Vzf_inv = VZf_new_inv - Vz_ident;
        
        Vyf=-Vyf;
        Vyf_inv = -Vyf_inv;
        
        disp('update transformations')
        toc
        
        tic
        
        phi = interp3(Vx_ident,Vy_ident,Vz_ident,phi_orig,VXf_new,VYf_new, VZf_new);
        phi(isnan(phi)) = 0;
        vecP = phi(:)';
        coef_P = img2coef3D(vecP,nx,ny,nz);
        CIP = coef_P(:)';
        coef_matP = reshape(CIP, ny+img_pU, nx+img_pV, nz+img_pW);
        
        save('phi_helix_final.mat','phi');
        Img1_orig = resize3Dmatrix(nx,ny,nz,Img1);
        Mimg = zeros(size(F));
        Fimg = zeros(size(F));
        Mimg(phi<=0) = 1;
        Fimg(Img1_orig>=128) = 1;
        
        dice_metric = evaluate_dicesimilarity3D(nx,ny,nz,Mimg,Fimg,2);
        fprintf("Dice metric = %f\n",dice_metric);
        
        disp('new images')
        toc
        tic
        
        path = pwd;
        
        filename_vtk = sprintf('/postprocessing/evolve%d.vtk',iterct);
        str =strcat(path,filename_vtk);
        vtkwrite(str, 'structured_grid',X_vtk,Y_vtk,Z_vtk,'scalars','Intensity',phi);
        
        p = patch(isosurface(X_vtk,Y_vtk,Z_vtk,F_orig,128));
        isonormals(X_vtk,Y_vtk,Z_vtk,F_orig,p);
        p.FaceColor = 'red';
        p.EdgeColor = 'none';
        daspect([1 1 1])
        view(3);
        axis tight
        camlight
        lighting gouraud
        hold on
        p1 = patch(isosurface(X_vtk,Y_vtk,Z_vtk,phi,0));
        isonormals(X_vtk,Y_vtk,Z_vtk,phi,p1);
        p1.FaceColor = 'green';
        p1.EdgeColor = 'none';
        daspect([1 1 1])
        view(3);
        axis tight
        hold on
        drawnow
        
        if(level==3)
            if (iteration==31)
                break;
            end
        end
        
        if(level==2)
            if (iteration==20)
                break;
            end
        end
        
        disp('display images')
        toc
        
    end
    
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