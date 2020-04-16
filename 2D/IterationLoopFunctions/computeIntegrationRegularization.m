function [RHSff] = computeIntegrationRegularization(parameters,H,BIGXf,BIGYf,Jm,Dm, PHII,RHSf)

%% This function evaluates the first order variation of the enrgy functional corresponding to hyperelastic regularization E_REG

%--Input variable:
%parameters,
%H: Variable storing element size of each active cell
%BIGXf: Variable storing gauss points (x-coordinate) for each active cell
%BIGYf:  Variable storing gauss points (y-coordinate) for each active cell
%Jm:  Struct variable storing supoort B-splines over each active cell
%Dm: THB-spline data structure element
%PHII: Sturct variable storing the basis function and its derivatives at
%each Gaussian point (phi, phi_x, phi_y)
%RHSf: Input RHS vector before regularization for update of control points

%--Output variable:
%RHSff: Output RHS vector after regularization for update of control points


ac_ct = size(Jm,1);
xlen = 6;
[~,w1] = ggquad(xlen);
w2 = w1;
BIGMUXf = BIGXf(:,:,2);
BIGMUYf = BIGXf(:,:,3);
BIGMVXf = BIGYf(:,:,2);
BIGMVYf = BIGYf(:,:,3);

PHUI = PHII{1,2};
PHVI = PHII{1,3};
for i = 1:ac_ct
    
    %length term
    term3f = BIGMUXf(1+(i-1)*xlen:i*xlen,1:6);
    term4f = BIGMUYf(1+(i-1)*xlen:i*xlen,1:6);
    term5f = BIGMVXf(1+(i-1)*xlen:i*xlen,1:6);
    term6f = BIGMVYf(1+(i-1)*xlen:i*xlen,1:6);

    SB = Jm{i,1};
    supp_phiu = PHUI{i,1};
    supp_phiv = PHVI{i,1};
    supp_size = size(SB,1);

    for bg = 1:supp_size
        valm1f = zeros(xlen,xlen);
        valm2f = zeros(xlen,xlen);

        for gg1 = 1:xlen
            for gg2 = 1:xlen

                phi_ui = supp_phiu(bg,gg1,gg2);
                phi_vi = supp_phiv(bg,gg1,gg2);

                %area term
                area = (term3f(gg1,gg2).^2 + term4f(gg1,gg2).^2)*(term5f(gg1,gg2).^2 + term6f(gg1,gg2).^2)-(term3f(gg1,gg2)*term5f(gg1,gg2) + term4f(gg1,gg2)*term6f(gg1,gg2));
                
                arxf = term3f(gg1,gg2)*(term5f(gg1,gg2).^2 + term6f(gg1,gg2).^2) - (term3f(gg1,gg2)*term5f(gg1,gg2) + term4f(gg1,gg2)*term6f(gg1,gg2))*term5f(gg1,gg2);
                brxf = term5f(gg1,gg2)*(term3f(gg1,gg2).^2 + term4f(gg1,gg2).^2) - (term3f(gg1,gg2)*term5f(gg1,gg2) + term4f(gg1,gg2)*term6f(gg1,gg2))*term3f(gg1,gg2);
                aryf = term4f(gg1,gg2)*(term5f(gg1,gg2).^2 + term6f(gg1,gg2).^2) - (term3f(gg1,gg2)*term5f(gg1,gg2) + term4f(gg1,gg2)*term6f(gg1,gg2))*term6f(gg1,gg2);
                bryf = term6f(gg1,gg2)*(term3f(gg1,gg2).^2 + term4f(gg1,gg2).^2) - (term3f(gg1,gg2)*term5f(gg1,gg2) + term4f(gg1,gg2)*term6f(gg1,gg2))*term4f(gg1,gg2);

                
                valm1f(gg1,gg2) = 2*parameters.lambda1*((term3f(gg1,gg2)-1)*phi_ui+term5f(gg1,gg2)*phi_vi)+2*max(area-1,0)*parameters.lambda2*(arxf*phi_ui+brxf*phi_vi);
                valm2f(gg1,gg2) = 2*parameters.lambda1*(term4f(gg1,gg2)*phi_ui+(term6f(gg1,gg2)-1)*phi_vi)+2*max(area-1,0)*parameters.lambda2*(aryf*phi_ui+bryf*phi_vi);
               
            end
        end
        
        h1 = H(i,1);
        h2 = H(i,2);
        val1f = w1' * valm1f * w2 * h1 * h2;
        val2f = w1' * valm2f * w2 * h1 * h2;

        btempp = Dm{SB(bg,2),1};
        bact_i = btempp{SB(bg,1),10};
        RHSf(bact_i,1) = RHSf(bact_i,1) + val1f;
        RHSf(bact_i,2) = RHSf(bact_i,2) + val2f;

    end
end

RHSff = RHSf;
end
