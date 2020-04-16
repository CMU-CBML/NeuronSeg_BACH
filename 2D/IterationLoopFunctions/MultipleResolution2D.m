function [phi,VF] = MultipleResolution2D(F,phi,param)

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

%--Output Variables
%phi: level set function capturing the final segmented image result
%VF: deformation field of the evolution of the level set function

%initialize phi, F
phi00 = phi;
F00 = F;

%% Basis value computed by substituting corresponding midpoint
Mb = [1/24, 11/24, 11/24, 1/24];

%% Bspline kernel
M_filt = Mb'*Mb;

%% Start multilevel
for level = param.maxlevel:-1:1
    
    disp(['Registration level:' num2str(param.maxlevel-level+1)]);
    
    %% downsample image
    scale = 2^-(level-1);  % image scale
    F = imresize(F,scale);
    npx = size(F,1); %size of scaled image
    npy = size(F,2);
    %vecF = F(:)';
    %coef_F= img2coef2D(npx,npy, vecF);
    %CIF = coef_F(:)';
    %coef_matF = reshape(CIF, npy+param.pU, npx+param.pV);
    
    phi = imresize(phi,scale);
    vecP = phi(:)';
    coef_P= img2coef2D(npx,npy, vecP);
    CIP = coef_P(:)';
    coef_matP = reshape(CIP, npy+param.pU, npx+param.pV);
    
    F_orig = imresize(F00,scale);
    vecForig = F_orig(:)';
    coef_Forig= img2coef2D(npx,npy, vecForig);
    CIForig = coef_Forig(:)';
    coef_matForig = reshape(CIForig, npy+param.pU, npx+param.pV);
    
    phi_orig = imresize(phi00,scale);
    vecphiorig = phi_orig(:)';
    coef_phiorig= img2coef2D(npx,npy, vecphiorig);
    CIPorig = coef_phiorig(:)';
    %coef_matPorig = reshape(CIPorig, npy+param.pU, npx+param.pV);
    
    if(level==param.maxlevel)
        
        nx = size(F,2);
        ny = size(F,1);
        
        [Vx_ident, Vy_ident] = meshgrid((0.5:nx-0.5),(0.5:ny-0.5));
        
        Vxlf = Vx_ident;
        Vylf = Vy_ident;
        
        Vxlf_inv = Vx_ident;
        Vylf_inv = Vy_ident;
        
    else
        
        nx = size(F,2);
        ny = size(F,1);
        
        [Vx_ident, Vy_ident] = meshgrid((0.5:nx-0.5),(0.5:ny-0.5));
        
        Vxlf = imresize(Vxf*scale,scale);
        Vylf = imresize(Vyf*scale,scale);
        
        Vxlf_inv = imresize(Vxf_inv*scale,scale);
        Vylf_inv = imresize(Vyf_inv*scale,scale);
        
        Vxlf = Vxlf+Vx_ident;
        Vylf = -Vylf+Vy_ident;
        
        Vxlf_inv = Vxlf_inv+Vx_ident;
        Vylf_inv = -Vylf_inv+Vy_ident;
        
        %II0F_ori = F_orig(:)';
        II0phi_ori = phi_orig(:)';
        
        %coef_IF = img2coef2D(nx,ny, II0F_ori);
        coef_Iphi = img2coef2D(nx,ny, II0phi_ori);
        
        CI0phi = coef_Iphi(:)';
        
        VXLf=Vxlf(:)';
        VYLf=Vylf(:)';
        
        phi = BsplineComposeImage2D(VXLf, VYLf, CI0phi, nx, ny);
        vecP = phi(:)';
        coef_P = img2coef2D(nx,ny, vecP);
        CIP = coef_P(:)';
        coef_matP = reshape(CIP, ny+param.pU, nx+param.pV);
    end
    
    VXLf = Vxlf(:)';
    VYLf = Vylf(:)';
    
    coef_xf = img2coef2D(nx,ny,VXLf);
    coef_yf = img2coef2D(nx,ny,VYLf);
    
    VXLf_inv = Vxlf_inv(:)';
    VYLf_inv = Vylf_inv(:)';
    
    coef_xf_inv = img2coef2D(nx,ny,VXLf_inv);
    coef_yf_inv = img2coef2D(nx,ny,VYLf_inv);
    
    
    %% Construct B-spline grid
    maxlev = param.maxlevel-level+1;
    
    [Dm,Pm,Em,~,knotvectorU,knotvectorV,nobU,nobV,nelemU] = setBsplineGrid(maxlev,param,F);
    Pmold = Pm;
    for multilev = 0:1:maxlev-1
        if(multilev>0)
            for j =1:bf_ct
                bbc = bf(j,1:2);
                bf_lev = bf(j,3);
                EEM = Em{bf_lev,1};
                BEM = Dm{bf_lev,1};
                bind = nobU(bf_lev,1)*(bbc(1,2)-1)+bbc(1,1);
                supp_cells = BEM{bind,6};
                grad = 0;
                supp_ct = 0;
                for i =1:size(supp_cells,2)
                    if(supp_cells(1,i)~=0)
                        supp_ct = supp_ct + 1;
                        ac_ind = EEM{supp_cells(1,i),11};
                        grad  = grad + Cell_grad(ac_ind,1);
                    end
                end
                
                grad = grad/supp_ct;
                
                %Refinement to create next level
                rho = param.rho(multilev+1);
                if(grad>=(rho*meangrad))
                    [~,~,Pmold] =  Refine2D(bbc(1,1),bbc(1,2),bf_lev,Dm,Em,Pmold,knotvectorU,knotvectorV,pU,pV);
                    [Dm,Em,Pm] =  Refine2D(bbc(1,1),bbc(1,2),bf_lev,Dm,Em,Pm,knotvectorU,knotvectorV,pU,pV);
                end
            end
        end
        
        Pmold = Pm;
        
        ac_ct = 0;
        bf_ct = 0;
        ac = zeros(1,2);
        bf = zeros(1,3);
        
        for lev = 1:(multilev+1)
            
            EE = Em{lev,1};
            BE = Dm{lev,1};
            sizee = size(EE,1);
            sizeb = size(BE,1);
            
            for i = 1:sizee
                if(EE{i,4}==1)
                    ac_ct = ac_ct+1;
                    ac(ac_ct,1) = EE{i,1};
                    ac(ac_ct,2) = lev;
                    EE{i,11} = ac_ct;
                end
            end
            
            for j = 1:sizeb
                if(BE{j,3}==1)
                    bf_ct = bf_ct + 1;
                    bf(bf_ct,1:2) = BE{j,1};
                    bf(bf_ct,3) = lev;
                    BE{j,10} = bf_ct;
                end
            end
            
            Em{lev,1} = EE;
            Dm{lev,1} = BE;
            
            
        end
        
        [Jm,Coeff,Pixel,HH,PHI,PHIU,PHIV,BIGX,BIGY] = constructAdaptiveGrid(ac,param,Dm,Em,F,knotvectorU,knotvectorV,multilev,nobU,nobV,nelemU);
        
        cell_co = zeros(ac_ct,2);
        for i = 1:ac_ct
            cell_id = ac(i,1);
            cell_le = ac(i,2);
            EEM = Em{cell_le,1};
            cell_co(i,1) = EEM{cell_id,8};
            cell_co(i,2) = EEM{cell_id,9};
        end
        
        F_255 = F.*255;
        [FDX,FDY] = gradient(F_255);
        Idiff = sqrt((FDX).^2 + (FDY).^2);
        Cell_grad = interp2(Vx_ident,Vy_ident,Idiff,cell_co(:,2),cell_co(:,1));
        meangrad = mean2(Idiff);
        displayAdaptiveGrid(ac,Coeff,Em,knotvectorU,knotvectorV,Jm,Pmold,param,nx,ny);
    end
    
    timestep = param.timestep(param.maxlevel-level+1);
    pU = param.pU;
    pV = param.pV;
    
    orderGauss=param.orderGauss;
    term_seg = zeros(ac_ct*orderGauss,orderGauss,3);
    
    K = fspecial('gaussian',2*round(2*param.sigma)+1,param.sigma);
    K_phi = fspecial('gaussian',5,param.sigma_phi);
    
    ngridX = 150;
    ngridY = 150;
    
    VGx = zeros(ngridY,ngridX);
    VGy = zeros(ngridY,ngridX);
    
    [FgridX, FgridY, ~, ~] = makeGrid(ngridX,ngridY,nx,ny);
    
    gridX0 = FgridX;
    gridY0 = FgridY;
    
    figure
    [cI0,~,~] = BsplineComposeImage2D_single(BIGY, BIGX, coef_matForig, size(BIGX,1), size(BIGX,2));
    
    for iteration = 1:param.maxiteration
        
        Pm = Pmold;
        Cm = Pmold;
        Bvect = zeros(bf_ct,2);
        
        Hphi = Heaviside(phi,param.epsilon);
        D1E = Delta(phi,param.epsilon);
        [f1,f2] = Local_Avr(F_orig,Hphi,K);
        
        vecD = D1E(:)';
        coef_D = img2coef2D(nx,ny, vecD);
        CID = coef_D(:)';
        coef_matD = reshape(CID, ny+3, nx+3);
        
        vecHphi = Hphi(:)';
        coef_Hphi = img2coef2D(nx,ny, vecHphi);
        CIHphi = coef_Hphi(:)';
        coef_matHphi = reshape(CIHphi, ny+3, nx+3);
        
        vec_f1 = f1(:)';
        coef_f1 = img2coef2D(nx,ny, vec_f1);
        CI_f1 = coef_f1(:)';
        coef_mat_f1 = reshape(CI_f1, ny+3, nx+3);
        
        vec_f2 = f2(:)';
        coef_f2 = img2coef2D(nx,ny, vec_f2);
        CI_f2 = coef_f2(:)';
        coef_mat_f2 = reshape(CI_f2, ny+3, nx+3);
        
        [~,cDphiX, cDphiY] = BsplineComposeImage2D_single(BIGY, BIGX, coef_matP, size(BIGX,1), size(BIGX,2));
        [dDelta, ~, ~] = BsplineComposeImage2D_single(BIGY, BIGX, coef_matD, size(BIGX,1), size(BIGX,2));
        [dHphi, ~, ~] = BsplineComposeImage2D_single(BIGY, BIGX, coef_matHphi, size(BIGX,1), size(BIGX,2));
        [df1, ~, ~] = BsplineComposeImage2D_single(BIGY, BIGX, coef_mat_f1, size(BIGX,1), size(BIGX,2));
        [df2, ~, ~] = BsplineComposeImage2D_single(BIGY, BIGX, coef_mat_f2, size(BIGX,1), size(BIGX,2));
        
        term_seg1 = -dDelta.*(cI0 -df1.*dHphi -df2.*(1-dHphi)).*(df1-df2);
        term_seg_x = cDphiX;
        term_seg_y = cDphiY;
        term_seg(:,:,1) = term_seg1;
        term_seg(:,:,2) = term_seg_x;
        term_seg(:,:,3) = term_seg_y;
        
        %Step1: Compute energy functional for segmentation
        [Bvect] = computeIntegrationFidelity(term_seg,Jm,PHI,Dm,HH,Bvect);
        
        Bvect = bcondition1(Bvect,Dm,bf,nobU);
        
        %Update the control points
        C1f = zeros(bf_ct,2);
        
        for i=1:bf_ct
            bbc = bf(i,1:2);
            bf_lev = bf(i,3);
            bi = nobU(bf_lev,1)*(bbc(1,2)-1)+bbc(1,1);
            Cif = Cm{bf_lev,1};
            
            C1f(i,1) = Cif(bi,1);
            C1f(i,2) = Cif(bi,2);
        end
        
        fprintf('timestep = %f, iteration = %d\n',timestep,iteration);
        
        C1f = C1f - timestep*Bvect; %update set of first control points
        
        for i = 1:bf_ct
            
            bbc = bf(i,1:2);
            bf_lev = bf(i,3);
            bi = nobU(bf_lev,1)*(bbc(1,2)-1)+bbc(1,1);
            
            Cif = Cm{bf_lev,1};
            Cif(bi,1) = C1f(i,1);
            Cif(bi,2) = C1f(i,2);
            Cm{bf_lev,1} = Cif;
        end
        
        %% Step 2: Regularization of the deformation field
        RHSf = zeros(bf_ct,2);
        
        [BIGXf,BIGYf] = computedeformation(Jm,PHI,PHIU,PHIV,Cm);
        
        PHII = cell(1,3);
        PHII{1,1} = PHI;
        PHII{1,2} = PHIU;
        PHII{1,3} = PHIV;
        
        [RHSff] = computeIntegrationRegularization(param,HH,BIGXf,BIGYf, Jm,Dm, PHII,RHSf);
        
        RHSff = bcondition1(RHSff,Dm,bf,nobU);
        
        %Update the control points
        P1f = zeros(bf_ct,2);
        
        for i=1:bf_ct
            bbc = bf(i,1:2);
            bf_lev = bf(i,3);
            bi = nobU(bf_lev,1)*(bbc(1,2)-1)+bbc(1,1);
            Pif = Cm{bf_lev,1};
            
            P1f(i,1) = Pif(bi,1);
            P1f(i,2) = Pif(bi,2);
            
        end
        
        P1fd = P1f - timestep*RHSff;
        
        for i = 1:bf_ct
            bbc = bf(i,1:2);
            bf_lev = bf(i,3);
            bi = nobU(bf_lev,1)*(bbc(1,2)-1)+bbc(1,1);
            
            Pif1 = Pm{bf_lev,1};
            Pif1(bi,1) = P1fd(i,1);
            Pif1(bi,2) = P1fd(i,2);
            Pm{bf_lev,1} = Pif1;
        end
        
        %Perform composition to update transformation function
        [PXF] = computenewpoints(F,Pixel,Jm,Pm);
        pxxf = PXF(:,:,1);
        pyyf = PXF(:,:,2);
        
        Vxf_temp = pxxf-Vx_ident;
        Vyf_temp = pyyf-Vy_ident;
        
        Vxf_old = Vxf_temp + Vx_ident;
        Vyf_old = Vyf_temp + Vy_ident;
        
        Vxf_old = reshape(Vxf_old,1,nx*ny);
        Vyf_old = reshape(Vyf_old,1,nx*ny);
        
        coef_xf_fun = reshape(coef_xf, 1, (nx+pU)*(ny+pV));
        coef_yf_fun = reshape(coef_yf, 1, (nx+pU)*(ny+pV));
        
        Vxf_old_inv = -Vxf_temp + Vx_ident;
        Vyf_old_inv = -Vyf_temp + Vy_ident;
        
        Vxf_old_inv = reshape(Vxf_old_inv,1,nx*ny);
        Vyf_old_inv = reshape(Vyf_old_inv,1,nx*ny);
        
        coef_xf_inv_fun = reshape(coef_xf_inv, 1, (nx+pU)*(ny+pV));
        coef_yf_inv_fun = reshape(coef_yf_inv, 1, (nx+pU)*(ny+pV));
        
        [VXf_new, VYf_new] = BsplineCompose2D( Vxf_old, Vyf_old, coef_xf_fun, coef_yf_fun, nx, ny);
        [VXf_new_inv, VYf_new_inv] = BsplineCompose2D( Vxf_old_inv, Vyf_old_inv, coef_xf_inv_fun, coef_yf_inv_fun, nx, ny);
        
        temp_coef_xf = imfilter(VXf_new, M_filt); %Vx  %doubt
        temp_coef_yf = imfilter(VYf_new, M_filt); %Vy
        
        coef_xf(4:ny,4:nx) = temp_coef_xf(2:end-2,2:end-2);
        coef_yf(4:ny,4:nx) = temp_coef_yf(2:end-2,2:end-2);
        
        temp_coef_xf_inv = imfilter(VXf_new_inv, M_filt); %Vx  %doubt
        temp_coef_yf_inv = imfilter(VYf_new_inv, M_filt); %Vy
        
        coef_xf_inv(4:ny,4:nx) = temp_coef_xf_inv(2:end-2,2:end-2);
        coef_yf_inv(4:ny,4:nx) = temp_coef_yf_inv(2:end-2,2:end-2);
        
        Vxf = VXf_new - Vx_ident;
        Vyf = VYf_new - Vy_ident;
        
        Vxf_inv = VXf_new_inv - Vx_ident;
        Vyf_inv = VYf_new_inv - Vy_ident;
        
        Vyf=-Vyf;
        Vyf_inv = -Vyf_inv;
        
        %Evaluate new phi
        phi = BsplineComposeImage2D(VXf_new(:)',VYf_new(:)', CIPorig, nx, ny);
        phi(isnan(phi)) = 0;
        vecP = phi(:)';
        coef_P = img2coef2D(nx,ny, vecP);
        CIP = coef_P(:)';
        coef_matP = reshape(CIP, ny+3, nx+3);
        phi = conv2(phi,K_phi,'same');
        
        for ii=1:ngridX
            for jj=1:ngridY
                
                j = gridY0(jj,ii);
                i = gridX0(jj,ii);
                
                VGx(jj,ii) = Vxf_inv(j,i);
                VGy(jj,ii) = Vyf_inv(j,i);
                
                FgridX(jj,ii) = round(i+VGx(jj,ii));
                FgridY(jj,ii) = round(j-VGy(jj,ii));
                
            end
        end
        
        subplot(1,2,1)
        imagesc(F);
        colormap gray;
        hold on
        contour(phi,[0,0],'g','Linewidth',3.24);
        hold off
        title('Segmented Image');
        
        th1 = subplot(1,2,2);
        cla(th1);
        plotGrid(FgridX, FgridY)
        axis ([1 nx 1 ny])
        colormap('gray')
        %set(gca,'position',[0 0 1 1],'units','normalized')
        set(gca,'YDir','reverse');
        hold off
        drawnow
        
    end
    
    figure
    imagesc(F);
    phi_img = 255.*ones(nx,ny);
    phi_img(phi<=0) = -255;
    hold on
    contour(phi_img,[0,0],'g','Linewidth',3.24);
    hold off
    set(gca,'position',[0 0 1 1],'units','normalized')
    colormap gray

    figure
    plotGrid(FgridX, FgridY)
    axis ([1 nx 1 ny])
    colormap('gray')
    set(gca,'position',[0 0 1 1],'units','normalized')
    set(gca,'YDir','reverse');
    hold off

    %upsample the image to original size
    Vxf = imresize(Vxf/scale,size(F00));
    Vyf = imresize(Vyf/scale,size(F00));
    
    Vxf_inv = imresize(Vxf_inv/scale,size(F00));
    Vyf_inv = imresize(Vyf_inv/scale,size(F00));
    
    F = imresize(F,size(F00));
    phi = imresize(phi,size(phi00));
end
VF = struct('VXf_new',VXf_new,'VXf_new_inv',VXf_new_inv,'VYf_new',VYf_new,'VYf_new_inv',VYf_new_inv);
end
