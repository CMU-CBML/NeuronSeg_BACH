function dice_metric = evaluate_dicesimilarity3D(nx, ny, nz, II0, II1 ,LabelSize)
%DiceSimilarity Compute DiceSimilarity between given segmented image (seg_result) and 
%expertly segmented image (expert_seg)

%Function uses the code: https://github.com/stellaccl/cdmffd-image-registration
%Paper: Chan, C. L., Anitescu, C., Zhang, Y., & Rabczuk, T. (2017). 
%Two and three dimensional image registration based on B-spline composition and level sets. 
%Communications in Computational Physics, 21(2), 600-622.

%Input: 
%        seg_result = segmented image result
%        expert_seg = expertly segmented image

%Output  dice_metric =  Dice similarity between image seg_result and expert_seg 

I0=reshape(II0, [ny nx nz]);
I1=reshape(II1, [ny nx nz]);


Dice=zeros(1,LabelSize);

for Label=0:LabelSize-1

  C1=numel(find(I0==Label));
            C2=numel(find(I1==Label));
            count=0;
            for i=1:nx
                for j=1:ny
                    for k =1:nz
                        if I0(j,i,k)==(Label) && I1(j,i,k)==(Label)
                            count=count+1;
                        end 
                    end
                end   
            end 
            Dice(Label+1)=(2*count)/(C1+C2);
end

dice_metric = Dice;

end