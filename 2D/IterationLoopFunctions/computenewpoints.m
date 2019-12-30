function [PXF] = computenewpoints(M,Pixel,Jm,Pm)

px = 0;
pxxf= zeros(size(M,1),size(M,2));
pyyf= zeros(size(M,1),size(M,2));

for i = 1:size(M,1),
    for j = 1:size(M,2),
        px = px+1;
        ac_ind = Pixel{px,1};
        supp = Pixel{px,2};
        SB = Jm{ac_ind,1};
        ss = size(SB,1);
        fxxf = 0;
        fyyf = 0;

        for k = 1:ss,
            
            CEbf = Pm{SB(k,2),1};
            pif = CEbf(SB(k,1),1);
            pjf = CEbf(SB(k,1),2);
            
            fxxf = fxxf + pif*supp(k,1);
            fyyf = fyyf + pjf*supp(k,1);

        end
        
        pxxf(i,j) = fxxf;
        pyyf(i,j) = fyyf;

    end
end
PXF(:,:,1) = pxxf;
PXF(:,:,2) = pyyf;

end