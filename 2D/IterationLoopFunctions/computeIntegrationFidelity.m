function [RHS] = computeIntegrationFidelity(Bterm_seg, Jm, PHI, Dm, H,Bvect)

%% This function evaluates the first order variation of the energy functional delta_E_LIF

%--Input Variable:
%Bterm_seg(:,:,1): The segmentation term at each Gaussian point 
%Bterm_seg(:,:,2) and Bterm_seg(:,:,3): x and y derivatives of phi. 
%Jm: Struct variable storing supoort B-splines over each active cell
%PHI:  Variable storing the basis function at each gaussian point in each
%active cell
%Dm: THB-spline data structure element
%H: Variable storing element size of each active cell
%Bvect: RHS vector to update control points before applying fidelity term

%--Output Variable:

%RHS: RHS vector to update control points after applying fidelity term


[~,w1] = ggquad(6);
xlen =6;
w2 = w1;
ac_ct = size(Jm,1);

for i = 1:ac_ct

    term_seg1 = Bterm_seg(1+(i-1)*xlen:i*xlen,1:6,1);
    term_segx = Bterm_seg(1+(i-1)*xlen:i*xlen,1:6,2);
    term_segy = Bterm_seg(1+(i-1)*xlen:i*xlen,1:6,3);
    
    SB = Jm{i,1};
    supp_phi = PHI{i,1};
    supp_size = size(SB,1);

    for bg = 1:supp_size

        valm1f = zeros(xlen,xlen);
        valm2f = zeros(xlen,xlen);

        for gg1 = 1:xlen
            for gg2 = 1:xlen
                
                phi_i = supp_phi(bg,gg1,gg2);
                
                valm1f(gg1,gg2) = phi_i*term_seg1(gg1,gg2)*term_segx(gg1,gg2);
                valm2f(gg1,gg2) = phi_i*term_seg1(gg1,gg2)*term_segy(gg1,gg2);

            end
        end
        
        h1 = H(i,1);
        h2 = H(i,2);
        
        val1f = w1' * valm1f * w2 * h1 * h2;
        val2f = w1' * valm2f * w2 * h1 * h2;

        btempp = Dm{SB(bg,2),1};
        bact_i = btempp{SB(bg,1),10};
        Bvect(bact_i,1) = Bvect(bact_i,1) + val1f;
        Bvect(bact_i,2) = Bvect(bact_i,2) + val2f;

    end
end

RHS = Bvect;

end
