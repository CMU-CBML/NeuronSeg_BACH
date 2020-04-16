function [BIGXf,BIGYf] = computedeformation(Jm,PHI,PHIU,PHIV,Pmf)

%% This function evaluates the deformed positions of the integration points 

%--Input Variable:
%Jm: Struct variable storing supoort B-splines over each active cell
%PHI: Struct variable storing the basis function at each Gaussian point (B(x))
%PHIU: Struct variable storing the x-derivative of basis function at each Gaussian point (B_u(x))
%PHIV: Struct variable storing the y-derivative of basis function at each Gaussian point (B_v(x))
%Pmf: Control points stored at each refinement level

%--Output Variable:
%BIGXf: Variable storing deformation and its derivatives (x-coordinate) at
%each Gaussian point
%BIGYf: Variable storing deformation and its derivatives (y-coordinate) at
%each Gaussian point

ac_ct = size(Jm,1);
xlen = 6;

BIGGXf = zeros(ac_ct*xlen,xlen,3);
BIGGYf = zeros(ac_ct*xlen,xlen,3);

BIGXXf = zeros(ac_ct*xlen,xlen);
BIGYYf = zeros(ac_ct*xlen,xlen);
BIGMUXf = zeros(ac_ct*xlen,xlen);
BIGMUYf = zeros(ac_ct*xlen,xlen);
BIGMVXf = zeros(ac_ct*xlen,xlen);
BIGMVYf = zeros(ac_ct*xlen,xlen);

for i = 1:ac_ct
    
    SB = Jm{i,1};
    supp_phi = PHI{i,1};
    supp_phiu = PHIU{i,1};
    supp_phiv = PHIV{i,1};
    supp_size = size(SB,1);
    
    %Now to compute the vectors, fx, fux, fvx
    SBXf = zeros(xlen,xlen);
    SBYf = zeros(xlen,xlen);
    SBUXf = zeros(xlen,xlen);
    SBUYf = zeros(xlen,xlen);
    SBVXf = zeros(xlen,xlen);
    SBVYf = zeros(xlen,xlen);

    for gg1 = 1:xlen
        for gg2 = 1:xlen
            
            sumbxf = 0;
            sumbyf = 0;
            sumbuxf = 0;
            sumbuyf = 0;
            sumbvxf = 0;
            sumbvyf = 0;

            for kg = 1:supp_size
                CEbf = Pmf{SB(kg,2),1};
                
                pif = CEbf(SB(kg,1),1);
                pjf = CEbf(SB(kg,1),2);

                sumbxf = sumbxf + pif*supp_phi(kg,gg1,gg2);
                sumbyf = sumbyf + pjf*supp_phi(kg,gg1,gg2);
                sumbuxf = sumbuxf + pif*supp_phiu(kg,gg1,gg2);
                sumbuyf = sumbuyf + pjf*supp_phiu(kg,gg1,gg2);
                sumbvxf = sumbvxf + pif*supp_phiv(kg,gg1,gg2);
                sumbvyf = sumbvyf + pjf*supp_phiv(kg,gg1,gg2);

            end
            
            SBXf(gg1,gg2) = sumbxf;
            SBYf(gg1,gg2) = sumbyf;
            SBUXf(gg1,gg2) = sumbuxf;
            SBUYf(gg1,gg2) = sumbuyf;
            SBVXf(gg1,gg2) = sumbvxf;
            SBVYf(gg1,gg2) = sumbvyf;

        end
    end
    
    BIGXXf(1+(i-1)*xlen:i*xlen,1:6) = SBXf;
    BIGYYf(1+(i-1)*xlen:i*xlen,1:6) = SBYf;
    BIGMUXf(1+(i-1)*xlen:i*xlen,1:6) = SBUXf;
    BIGMUYf(1+(i-1)*xlen:i*xlen,1:6) = SBUYf;
    BIGMVXf(1+(i-1)*xlen:i*xlen,1:6) = SBVXf;
    BIGMVYf(1+(i-1)*xlen:i*xlen,1:6) = SBVYf;
    
end

BIGGXf(:,:,1) = BIGXXf;
BIGGXf(:,:,2) = BIGMUXf;
BIGGXf(:,:,3) = BIGMUYf;

BIGGYf(:,:,1) = BIGYYf;
BIGGYf(:,:,2) = BIGMVXf;
BIGGYf(:,:,3) = BIGMVYf;


BIGXf = BIGGXf;
BIGYf = BIGGYf;


end
