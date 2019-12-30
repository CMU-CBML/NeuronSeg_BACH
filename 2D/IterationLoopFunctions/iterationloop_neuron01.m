for iteration = 1:param.maxiteration
    
    Pm = Pmold;
    Cm = Pmold;
    
    rs_initial = rs_final;
    
    Hphi = Heaviside(phi,parameters.epsilon);
    D1E = Delta(phi,parameters.epsilon);
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
    
    Bvectf = zeros(bf_ct,2);
    
    [cI1,cDI1X, cDI1Y] = BsplineComposeImage2D_single(BIGY, BIGX, coef_matM, size(BIGX,1), size(BIGX,2));
    [ctempp,cDphiX, cDphiY] = BsplineComposeImage2D_single(BIGY, BIGX, coef_matP, size(BIGX,1), size(BIGX,2));
    [dDelta, ~, ~] = BsplineComposeImage2D_single(BIGY, BIGX, coef_matD, size(BIGX,1), size(BIGX,2));
    [dHphi, ~, ~] = BsplineComposeImage2D_single(BIGY, BIGX, coef_matHphi, size(BIGX,1), size(BIGX,2));
    [df1, ~, ~] = BsplineComposeImage2D_single(BIGY, BIGX, coef_mat_f1, size(BIGX,1), size(BIGX,2));
    [df2, ~, ~] = BsplineComposeImage2D_single(BIGY, BIGX, coef_mat_f2, size(BIGX,1), size(BIGX,2));
    
    denominator = sqrt(cDI1X.^2+cDI1Y.^2 + smallNumber);
    
    term_reg_x = (cI1 - cI0).*2.*cDI1X./denominator;
    term_reg_y = (cI1 - cI0).*2.*cDI1Y./denominator;
    term_reg(:,:,1) = term_reg_x;
    term_reg(:,:,2) = term_reg_y;
    
    term_seg1 = -dDelta.*(cI0 -df1.*dHphi -df2.*(1-dHphi)).*(df1-df2);
    term_seg_x = cDphiX;
    term_seg_y = cDphiY;
    term_seg(:,:,1) = term_seg1;
    term_seg(:,:,2) = term_seg_x;
    term_seg(:,:,3) = term_seg_y;
    
    [Bvectf] = computeIntegrationFidelity(term_seg,Jm,PHI,Dm,HH,Bvectf);
    
    Bvectf = bcondition1(Bvectf,Dm,bf,nobU);
    
    %Update the control points
    C1f = zeros(bf_ct,2);
    
    for i=1:bf_ct,
        bbc = bf(i,1:2);
        bf_lev = bf(i,3);
        bi = nobU(bf_lev,1)*(bbc(1,2)-1)+bbc(1,1);
        Cif = Cm{bf_lev,1};
        
        C1f(i,1) = Cif(bi,1);
        C1f(i,2) = Cif(bi,2);
    end
    
    timestep
    C1f = C1f - timestep*Bvectf; %update set of first control points
    
    for i = 1:bf_ct
        
        bbc = bf(i,1:2);
        bf_lev = bf(i,3);
        bi = nobU(bf_lev,1)*(bbc(1,2)-1)+bbc(1,1);
        
        Cif = Cm{bf_lev,1};
        Cif(bi,1) = C1f(i,1);
        Cif(bi,2) = C1f(i,2);
        Cm{bf_lev,1} = Cif;
    end
    
    %% now do the regularization
    RHSf = zeros(bf_ct,2);
    
    [BIGXf,BIGYf] = computedeformation(Jm,PHI,PHIU,PHIV,Cm);
    
    PHII = cell(1,3);
    PHII{1,1} = PHI;
    PHII{1,2} = PHIU;
    PHII{1,3} = PHIV;
    
    BIGXXf = BIGXf(:,:,1);
    BIGMUXf = BIGXf(:,:,2);
    BIGMUYf = BIGXf(:,:,3);
    
    BIGYYf = BIGYf(:,:,1);
    BIGMVXf = BIGYf(:,:,2);
    BIGMVYf = BIGYf(:,:,3);
    
    [RHSff] = computeIntegrationRegularization(parameters,HH,BIGXf,BIGYf, Jm,Dm, PHII,RHSf);
    
    RHSff = bcondition1(RHSff,Dm,bf,nobU);
    
    %Update the control points
    P1f = zeros(bf_ct,2);
    
    for i=1:bf_ct,
        bbc = bf(i,1:2);
        bf_lev = bf(i,3);
        bi = nobU(bf_lev,1)*(bbc(1,2)-1)+bbc(1,1);
        Pif = Cm{bf_lev,1};
        
        P1f(i,1) = Pif(bi,1);
        P1f(i,2) = Pif(bi,2);
        
    end
    
    P1fd = P1f - timestep*RHSff;
    
    for i = 1:bf_ct,
        bbc = bf(i,1:2);
        bf_lev = bf(i,3);
        bi = nobU(bf_lev,1)*(bbc(1,2)-1)+bbc(1,1);
        
        Pif1 = Pm{bf_lev,1};
        Pif1(bi,1) = P1fd(i,1);
        Pif1(bi,2) = P1fd(i,2);
        Pm{bf_lev,1} = Pif1;
    end
    
    [PXF] = computenewpoints(M,Pixel,Jm,Pm);
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
    
    M = BsplineComposeImage2D(VXf_new(:)',VYf_new(:)', CIMorig, nx, ny);
    M(isnan(M)) = 0;
    vecM = M(:)';
    coef_M= img2coef2D(nx,ny, vecM);
    CIM = coef_M(:)';
    coef_matM = reshape(CIM, ny+3, nx+3);
    
    phi = BsplineComposeImage2D(VXf_new(:)',VYf_new(:)', CIPorig, nx, ny);
    phi(isnan(phi)) = 0;
    vecP = phi(:)';
    coef_P = img2coef2D(nx,ny, vecP);
    CIP = coef_P(:)';
    coef_matP = reshape(CIP, ny+3, nx+3);
    
    M = imgaussfilt(M,param.sigma_smooth);
    phi = imgaussfilt(phi,param.sigma_smooth);
    
    [MDX,MDY] = gradient(M);
    [DphiX,DphiY] = gradient(phi);
    
    phi_img = 255.*ones(nx,ny);
    phi_img(phi<=5) = -255;
    subplot(1,2,1);
    imagesc(phi);
    hold on
    contour(phi_img,[0,0],'g');
    hold off
    title('Source Image');
    axis([1 nx 1 ny])
    colormap gray
    
    colormap gray;
    subplot(1,2,2);
    imagesc(F);
    hold on
    contour(phi_img,[0,0],'g');
    hold off
    title('Segmented Image');
    axis([1 nx 1 ny])
    drawnow
    
end