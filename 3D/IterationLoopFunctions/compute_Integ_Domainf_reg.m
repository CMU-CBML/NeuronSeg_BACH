function RHS = compute_Integ_Domainf_reg(Jm,BIGMUXf, BIGMUYf, BIGMUZf, BIGMVXf, BIGMVYf, BIGMVZf, BIGMWXf, BIGMWYf, BIGMWZf,RHSf,PHIU1, PHIV1, PHIW1, lambda_1,lambda_2, lambda_3, w1,w2,w3,H)

ac_ct = size(Jm,1);
bf_ct = size(RHSf,1);
xlen = 4;

parfor i = 1:ac_ct
    
    RHSf1 = zeros(bf_ct,4);
    
    SB = Jm(i).nzsplines;

    supp_phiu = PHIU1(i).mat;
    supp_phiv = PHIV1(i).mat;
    supp_phiw = PHIW1(i).mat;
    supp_size = size(SB,1);

    hu = H(i,1);
    hv = H(i,2);
    hw = H(i,3);

    term_fxx = BIGMUXf(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    term_fxy = BIGMUYf(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    term_fxz = BIGMUZf(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    
    term_fyx = BIGMVXf(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    term_fyy = BIGMVYf(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    term_fyz = BIGMVZf(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    
    term_fzx = BIGMWXf(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    term_fzy = BIGMWYf(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    term_fzz = BIGMWZf(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    
    val1f = zeros(supp_size,1);
    val2f = zeros(supp_size,1);
    val3f = zeros(supp_size,1);
    
    reg_length_x = zeros(supp_size,xlen,xlen,xlen);
    reg_length_y = zeros(supp_size,xlen,xlen,xlen);
    reg_length_z = zeros(supp_size,xlen,xlen,xlen);
    
    reg_area_x = zeros(supp_size,xlen,xlen,xlen);
    reg_area_y = zeros(supp_size,xlen,xlen,xlen);
    reg_area_z = zeros(supp_size,xlen,xlen,xlen);
    
    reg_volume_x = zeros(supp_size,xlen,xlen,xlen);
    reg_volume_y = zeros(supp_size,xlen,xlen,xlen);
    reg_volume_z = zeros(supp_size,xlen,xlen,xlen);
    
    for gg1 = 1:xlen
        for gg2 = 1:xlen
            for gg3 = 1:xlen

                phi_ui = supp_phiu(:,gg1,gg2,gg3);
                phi_vi = supp_phiv(:,gg1,gg2,gg3);
                phi_wi = supp_phiw(:,gg1,gg2,gg3);

                reg_length_x(:,gg1,gg2,gg3) = 2.*lambda_1.*(term_fxx(gg1,gg2,gg3)-1).*phi_ui+term_fxy(gg1,gg2,gg3).*phi_vi+term_fxz(gg1,gg2,gg3).*phi_wi;
                reg_length_y(:,gg1,gg2,gg3) = 2.*lambda_1.*term_fyx(gg1,gg2,gg3).*phi_ui+(term_fyy(gg1,gg2,gg3)-1).*phi_vi+term_fyz(gg1,gg2,gg3).*phi_wi;
                reg_length_z(:,gg1,gg2,gg3) = 2.*lambda_1.*term_fzx(gg1,gg2,gg3).*phi_ui+term_fzy(gg1,gg2,gg3).*phi_vi+(term_fzz(gg1,gg2,gg3)-1).*phi_wi;
         
                surface_area = 2*((term_fxx(gg1,gg2,gg3)^2+term_fyx(gg1,gg2,gg3)^2+term_fzx(gg1,gg2,gg3)^2)*(term_fxy(gg1,gg2,gg3)^2+term_fyy(gg1,gg2,gg3)^2+term_fzy(gg1,gg2,gg3)^2)-(term_fxx(gg1,gg2,gg3)*term_fxy(gg1,gg2,gg3)+term_fyx(gg1,gg2,gg3)*term_fyy(gg1,gg2,gg3)+term_fzx(gg1,gg2,gg3)*term_fzy(gg1,gg2,gg3)))+...
                    2*((term_fxy(gg1,gg2,gg3)^2+term_fyy(gg1,gg2,gg3)^2+term_fzy(gg1,gg2,gg3)^2)*(term_fxz(gg1,gg2,gg3)^2+term_fyz(gg1,gg2,gg3)^2+term_fzz(gg1,gg2,gg3)^2)-(term_fxy(gg1,gg2,gg3)*term_fxz(gg1,gg2,gg3)+term_fyy(gg1,gg2,gg3)*term_fyz(gg1,gg2,gg3)+term_fzy(gg1,gg2,gg3)*term_fzz(gg1,gg2,gg3)))+...
                    2*((term_fxx(gg1,gg2,gg3)^2+term_fyx(gg1,gg2,gg3)^2+term_fzx(gg1,gg2,gg3)^2)*(term_fxz(gg1,gg2,gg3)^2+term_fyz(gg1,gg2,gg3)^2+term_fzz(gg1,gg2,gg3)^2)-(term_fxx(gg1,gg2,gg3)*term_fxz(gg1,gg2,gg3)+term_fyx(gg1,gg2,gg3)*term_fyz(gg1,gg2,gg3)+term_fzx(gg1,gg2,gg3)*term_fzz(gg1,gg2,gg3)));
                
                reg_area_x(:,gg1,gg2,gg3) = (4.*lambda_2.*max(surface_area-1,0)).*(term_fxx(gg1,gg2,gg3).*(term_fxy(gg1,gg2,gg3).^2+term_fyy(gg1,gg2,gg3).^2+term_fzy(gg1,gg2,gg3).^2)).*phi_ui+(4*lambda_2*max(surface_area-1,0)).*(term_fxy(gg1,gg2,gg3).*(term_fxx(gg1,gg2,gg3).^2+term_fyx(gg1,gg2,gg3).^2+term_fzx(gg1,gg2,gg3).^2)).*phi_vi-(4*lambda_2*max(surface_area-1,0)).*((term_fxx(gg1,gg2,gg3).*term_fxy(gg1,gg2,gg3)+term_fyx(gg1,gg2,gg3).*term_fyy(gg1,gg2,gg3)+term_fzx(gg1,gg2,gg3).*term_fzy(gg1,gg2,gg3)).*(term_fxx(gg1,gg2,gg3).*phi_ui+term_fxy(gg1,gg2,gg3).*phi_vi))+...
                     (4.*lambda_2.*max(surface_area-1,0)).*(term_fxy(gg1,gg2,gg3).*(term_fxz(gg1,gg2,gg3).^2+term_fyz(gg1,gg2,gg3).^2+term_fzz(gg1,gg2,gg3).^2)).*phi_vi+(4.*lambda_2.*max(surface_area-1,0)).*(term_fxz(gg1,gg2,gg3).*(term_fxy(gg1,gg2,gg3).^2+term_fyy(gg1,gg2,gg3).^2+term_fzy(gg1,gg2,gg3).^2)).*phi_wi-(4*lambda_2*max(surface_area-1,0)).*((term_fxy(gg1,gg2,gg3).*term_fxz(gg1,gg2,gg3)+term_fyy(gg1,gg2,gg3).*term_fyz(gg1,gg2,gg3)+term_fzy(gg1,gg2,gg3).*term_fzz(gg1,gg2,gg3)).*(term_fxy(gg1,gg2,gg3).*phi_vi+term_fxz(gg1,gg2,gg3).*phi_wi))+...
                     (4.*lambda_2.*max(surface_area-1,0)).*(term_fxz(gg1,gg2,gg3).*(term_fxx(gg1,gg2,gg3).^2+term_fyx(gg1,gg2,gg3).^2+term_fzx(gg1,gg2,gg3).^2)).*phi_wi+(4.*lambda_2.*max(surface_area-1,0)).*(term_fxx(gg1,gg2,gg3).*(term_fxz(gg1,gg2,gg3).^2+term_fyz(gg1,gg2,gg3).^2+term_fzz(gg1,gg2,gg3).^2)).*phi_ui-(4*lambda_2*max(surface_area-1,0)).*((term_fxz(gg1,gg2,gg3).*term_fxx(gg1,gg2,gg3)+term_fyz(gg1,gg2,gg3).*term_fyx(gg1,gg2,gg3)+term_fzz(gg1,gg2,gg3).*term_fzx(gg1,gg2,gg3)).*(term_fxz(gg1,gg2,gg3).*phi_wi+term_fxx(gg1,gg2,gg3)*phi_ui));
                 
                reg_area_y(:,gg1,gg2,gg3) = (4.*lambda_2.*max(surface_area-1,0)).*(term_fyx(gg1,gg2,gg3).*(term_fxy(gg1,gg2,gg3).^2+term_fyy(gg1,gg2,gg3).^2+term_fzy(gg1,gg2,gg3).^2)).*phi_ui+(4*lambda_2*max(surface_area-1,0)).*(term_fyy(gg1,gg2,gg3).*(term_fxx(gg1,gg2,gg3).^2+term_fyx(gg1,gg2,gg3).^2+term_fzx(gg1,gg2,gg3).^2)).*phi_vi-(4*lambda_2*max(surface_area-1,0)).*((term_fxx(gg1,gg2,gg3).*term_fxy(gg1,gg2,gg3)+term_fyx(gg1,gg2,gg3).*term_fyy(gg1,gg2,gg3)+term_fzx(gg1,gg2,gg3).*term_fzy(gg1,gg2,gg3)).*(term_fyx(gg1,gg2,gg3).*phi_ui+term_fyy(gg1,gg2,gg3).*phi_vi))+...
                     (4.*lambda_2.*max(surface_area-1,0)).*(term_fyy(gg1,gg2,gg3).*(term_fxz(gg1,gg2,gg3).^2+term_fyz(gg1,gg2,gg3).^2+term_fzz(gg1,gg2,gg3).^2)).*phi_vi+(4.*lambda_2.*max(surface_area-1,0)).*(term_fyz(gg1,gg2,gg3).*(term_fxy(gg1,gg2,gg3).^2+term_fyy(gg1,gg2,gg3).^2+term_fzy(gg1,gg2,gg3).^2)).*phi_wi-(4*lambda_2*max(surface_area-1,0)).*((term_fxy(gg1,gg2,gg3).*term_fxz(gg1,gg2,gg3)+term_fyy(gg1,gg2,gg3).*term_fyz(gg1,gg2,gg3)+term_fzy(gg1,gg2,gg3).*term_fzz(gg1,gg2,gg3)).*(term_fyy(gg1,gg2,gg3).*phi_vi+term_fyz(gg1,gg2,gg3).*phi_wi))+...
                     (4.*lambda_2.*max(surface_area-1,0)).*(term_fyz(gg1,gg2,gg3).*(term_fxx(gg1,gg2,gg3).^2+term_fyx(gg1,gg2,gg3).^2+term_fzx(gg1,gg2,gg3).^2)).*phi_wi+(4.*lambda_2.*max(surface_area-1,0)).*(term_fyx(gg1,gg2,gg3).*(term_fxz(gg1,gg2,gg3).^2+term_fyz(gg1,gg2,gg3).^2+term_fzz(gg1,gg2,gg3).^2)).*phi_ui-(4*lambda_2*max(surface_area-1,0)).*((term_fxz(gg1,gg2,gg3).*term_fxx(gg1,gg2,gg3)+term_fyz(gg1,gg2,gg3).*term_fyx(gg1,gg2,gg3)+term_fzz(gg1,gg2,gg3).*term_fzx(gg1,gg2,gg3)).*(term_fyz(gg1,gg2,gg3).*phi_wi+term_fyx(gg1,gg2,gg3)*phi_ui));
                
                reg_area_z(:,gg1,gg2,gg3) = (4.*lambda_2.*max(surface_area-1,0)).*(term_fzx(gg1,gg2,gg3).*(term_fxy(gg1,gg2,gg3).^2+term_fyy(gg1,gg2,gg3).^2+term_fzy(gg1,gg2,gg3).^2)).*phi_ui+(4*lambda_2*max(surface_area-1,0)).*(term_fzy(gg1,gg2,gg3).*(term_fxx(gg1,gg2,gg3).^2+term_fyx(gg1,gg2,gg3).^2+term_fzx(gg1,gg2,gg3).^2)).*phi_vi-(4*lambda_2*max(surface_area-1,0)).*((term_fxx(gg1,gg2,gg3).*term_fxy(gg1,gg2,gg3)+term_fyx(gg1,gg2,gg3).*term_fyy(gg1,gg2,gg3)+term_fzx(gg1,gg2,gg3).*term_fzy(gg1,gg2,gg3)).*(term_fzx(gg1,gg2,gg3).*phi_ui+term_fzy(gg1,gg2,gg3).*phi_vi))+...
                     (4.*lambda_2.*max(surface_area-1,0)).*(term_fzy(gg1,gg2,gg3).*(term_fxz(gg1,gg2,gg3).^2+term_fyz(gg1,gg2,gg3).^2+term_fzz(gg1,gg2,gg3).^2)).*phi_vi+(4.*lambda_2.*max(surface_area-1,0)).*(term_fzz(gg1,gg2,gg3).*(term_fxy(gg1,gg2,gg3).^2+term_fyy(gg1,gg2,gg3).^2+term_fzy(gg1,gg2,gg3).^2)).*phi_wi-(4*lambda_2*max(surface_area-1,0)).*((term_fxy(gg1,gg2,gg3).*term_fxz(gg1,gg2,gg3)+term_fyy(gg1,gg2,gg3).*term_fyz(gg1,gg2,gg3)+term_fzy(gg1,gg2,gg3).*term_fzz(gg1,gg2,gg3)).*(term_fzy(gg1,gg2,gg3).*phi_vi+term_fzz(gg1,gg2,gg3).*phi_wi))+...
                     (4.*lambda_2.*max(surface_area-1,0)).*(term_fzz(gg1,gg2,gg3).*(term_fxx(gg1,gg2,gg3).^2+term_fyx(gg1,gg2,gg3).^2+term_fzx(gg1,gg2,gg3).^2)).*phi_wi+(4.*lambda_2.*max(surface_area-1,0)).*(term_fzx(gg1,gg2,gg3).*(term_fxz(gg1,gg2,gg3).^2+term_fyz(gg1,gg2,gg3).^2+term_fzz(gg1,gg2,gg3).^2)).*phi_ui-(4*lambda_2*max(surface_area-1,0)).*((term_fxz(gg1,gg2,gg3).*term_fxx(gg1,gg2,gg3)+term_fyz(gg1,gg2,gg3).*term_fyx(gg1,gg2,gg3)+term_fzz(gg1,gg2,gg3).*term_fzx(gg1,gg2,gg3)).*(term_fzz(gg1,gg2,gg3).*phi_wi+term_fzx(gg1,gg2,gg3)*phi_ui));
                
                determinant = term_fxx(gg1,gg2,gg3)*(term_fyy(gg1,gg2,gg3)*term_fzz(gg1,gg2,gg3)-term_fyz(gg1,gg2,gg3)*term_fzy(gg1,gg2,gg3))-term_fyy(gg1,gg2,gg3)*(term_fyx(gg1,gg2,gg3)*term_fzz(gg1,gg2,gg3)-term_fyz(gg1,gg2,gg3)*term_fzx(gg1,gg2,gg3))+term_fxz(gg1,gg2,gg3)*(term_fyx(gg1,gg2,gg3)*term_fzy(gg1,gg2,gg3)+term_fyy(gg1,gg2,gg3)*term_fzx(gg1,gg2,gg3));
                
                reg_volume_x(:,gg1,gg2,gg3) = (2*lambda_3*(determinant+1)*((determinant-1)/determinant)^3).*((term_fxx(gg1,gg2,gg3)*term_fzz(gg1,gg2,gg3)-term_fyz(gg1,gg2,gg3)*term_fzy(gg1,gg2,gg3)).*phi_ui-(term_fyx(gg1,gg2,gg3)*term_fzz(gg1,gg2,gg3)-term_fyz(gg1,gg2,gg3)*term_fzx(gg1,gg2,gg3)).*phi_vi+(term_fyx(gg1,gg2,gg3)*term_fzy(gg1,gg2,gg3)-term_fyy(gg1,gg2,gg3)*term_fzx(gg1,gg2,gg3)).*phi_wi);
                reg_volume_y(:,gg1,gg2,gg3) = (2*lambda_3*(determinant+1)*((determinant-1)/determinant)^3).*((term_fxz(gg1,gg2,gg3)*term_fzy(gg1,gg2,gg3)-term_fxy(gg1,gg2,gg3)*term_fzz(gg1,gg2,gg3)).*phi_ui-(term_fxx(gg1,gg2,gg3)*term_fzz(gg1,gg2,gg3)-term_fxz(gg1,gg2,gg3)*term_fzx(gg1,gg2,gg3)).*phi_vi+(term_fzx(gg1,gg2,gg3)*term_fxy(gg1,gg2,gg3)-term_fzy(gg1,gg2,gg3)*term_fxx(gg1,gg2,gg3)).*phi_wi);
                reg_volume_z(:,gg1,gg2,gg3) = (2*lambda_3*(determinant+1)*((determinant-1)/determinant)^3).*((term_fxy(gg1,gg2,gg3)*term_fyz(gg1,gg2,gg3)-term_fyz(gg1,gg2,gg3)*term_fzy(gg1,gg2,gg3)).*phi_ui-(term_fxz(gg1,gg2,gg3)*term_fyx(gg1,gg2,gg3)-term_fxx(gg1,gg2,gg3)*term_fyz(gg1,gg2,gg3)).*phi_vi+(term_fxx(gg1,gg2,gg3)*term_fyy(gg1,gg2,gg3)-term_fxy(gg1,gg2,gg3)*term_fyx(gg1,gg2,gg3)).*phi_wi);
                
                val1f(:,1) = val1f(:,1) + w1(gg1,1).*w2(gg2,1).*w3(gg3,1).*reg_length_x(:,gg1,gg2,gg3).*hu.*hv.*hw + w1(gg1,1).*w2(gg2,1).*w3(gg3,1).*reg_area_x(:,gg1,gg2,gg3).*hu.*hv.*hw + w1(gg1,1).*w2(gg2,1).*w3(gg3,1).*reg_volume_x(:,gg1,gg2,gg3).*hu.*hv.*hw;
                val2f(:,1) = val2f(:,1) + w1(gg1,1).*w2(gg2,1).*w3(gg3,1).*reg_length_y(:,gg1,gg2,gg3).*hu.*hv.*hw + w1(gg1,1).*w2(gg2,1).*w3(gg3,1).*reg_area_y(:,gg1,gg2,gg3).*hu.*hv.*hw + w1(gg1,1).*w2(gg2,1).*w3(gg3,1).*reg_volume_y(:,gg1,gg2,gg3).*hu.*hv.*hw;
                val3f(:,1) = val3f(:,1) + w1(gg1,1).*w2(gg2,1).*w3(gg3,1).*reg_length_z(:,gg1,gg2,gg3).*hu.*hv.*hw + w1(gg1,1).*w2(gg2,1).*w3(gg3,1).*reg_area_z(:,gg1,gg2,gg3).*hu.*hv.*hw + w1(gg1,1).*w2(gg2,1).*w3(gg3,1).*reg_volume_z(:,gg1,gg2,gg3).*hu.*hv.*hw;
                
            end
        end
    end
    
    RHSf1(SB,1) = val1f;
    RHSf1(SB,2) = val2f;
    RHSf1(SB,3) = val3f;
    
    RHSf = RHSf  + RHSf1;
    
end

RHS = RHSf;
end


