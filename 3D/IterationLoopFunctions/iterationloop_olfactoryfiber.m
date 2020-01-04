K = fspecial3('gaussian',2*round(2*param.sigma)+1,param.sigma);
K_phi = fspecial3('gaussian',5,param.sigma_phi);

%gaussian filter
sigma_par = round(3*scale);
smooth_par = 2*sigma_par;
Hsmooth=fspecial3('gaussian',smooth_par,sigma_par);

for iteration = 1:param.maxiteration
    
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
    iteration = iteration+1;
    disp('iterct:');
    disp(iterct);
    toc
    tic
    [ctempp,cDphiY, cDphiX,cDphiZ] =  BsplineComposeImage3D_single(BIGY, BIGX, BIGZ, coef_matP, size(BIGX,1), size(BIGX,2), size(BIGX,3));
    [dDelta, dt, dt1, dt2] =  BsplineComposeImage3D_single(BIGY, BIGX, BIGZ, coef_matD, size(BIGX,1), size(BIGX,2), size(BIGX,3));
    [dHphi, dt, dt1, dt2] =  BsplineComposeImage3D_single(BIGY, BIGX, BIGZ, coef_matH, size(BIGX,1), size(BIGX,2), size(BIGX,3));
    [df1, dt, dt1, dt2] =  BsplineComposeImage3D_single(BIGY, BIGX, BIGZ, coef_matf1, size(BIGX,1), size(BIGX,2), size(BIGX,3));
    [df2, dt, dt1, dt2] =  BsplineComposeImage3D_single(BIGY, BIGX, BIGZ, coef_matf2, size(BIGX,1), size(BIGX,2), size(BIGX,3));
    
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
    
    RHS = compute_Integ_Domain_fidelity(Jm,Bseg, Bsegx,Bsegy,Bsegz,RHS,PHI1,Wu,Wv,Ww,H);
    
    toc
    tic
    RHS = bcondition3D(RHS);
    
    disp('integration done')
    toc
    tic
    timestep
    ACCf(:,1:3) = ACCf(:,1:3) - timestep.*RHS(:,1:3);
    
    tic
    [BIGXXf, BIGYYf, BIGZZf, BIGMUXf, BIGMUYf, BIGMUZf, BIGMVXf, BIGMVYf, BIGMVZf, BIGMWXf, BIGMWYf, BIGMWZf] =  computenewPoints(Jm, ACCf, PHI1, PHIU1, PHIV1, PHIW1, orderGauss);
    
    disp('BIGXX computed');
    toc
    
    [s1,w1] = ggquad(orderGauss);
    [s2,w2] = ggquad(orderGauss);
    [s3,w3] = ggquad(orderGauss);

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

    VXFF_new = VXf_new(:)';
    VYFF_new = VYf_new(:)';
    VZFF_new = VZf_new(:)';

    phi = interp3(Vx_ident,Vy_ident,Vz_ident,phi_orig,VXf_new,VYf_new, VZf_new);
    phi(isnan(phi)) = 0;
    vecP = phi(:)';
    coef_P = img2coef3D(vecP,nx,ny,nz);
    CIP = coef_P(:)';
    coef_matP = reshape(CIP, ny+img_pU, nx+img_pV, nz+img_pW);

    filename_phi= sprintf('phi%d.mat',iterct);
    save(filename_phi,'phi');
    
    disp('new images')
    toc
    tic
     
    p = patch(isosurface(X_vtk,Y_vtk,Z_vtk,F_orig,64));
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
    camlight
    lighting gouraud
    hold off
    drawnow

    [FDX,FDY,FDZ] = gradient(F);
    [DphiX,DphiY,DphiZ] = gradient(phi);
    
    if(level==1)
    [xy_accuracy,z_accuracy] = evaluate_accuracy_olfactoryfiber(phi,tree);
    disp('xy_accuracy');
    disp(xy_accuracy);
    disp('z_accuracy');
    disp(z_accuracy);
    end
    
    disp('display images')
    toc
    
end