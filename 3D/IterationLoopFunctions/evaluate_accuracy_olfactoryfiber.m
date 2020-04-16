function [xy_Accuracy, z_Accuracy] = evaluate_accuracy_olfactoryfiber(phi,tree)

%% Compare accuracy with gold standard

%--Input Variable:
%phi: final segmentation result
%tree: ground truth tracing tree (SWC file) for checking accuracy

%--Output Variable:
%xy_Accuracy: the accuracy of matched points with ground truth in XY plane
%z_Accuracy: the accuracy of matched points with ground truth in Z axis

phi_img = zeros(size(phi));
phi_img(phi<=0) = 1;
OPVolLogical = logical(phi_img);
OPVolSkel = bwskel(OPVolLogical);

mytree_x = [];
mytree_y = [];
mytree_z = [];

for i=1:size(phi,1)
    for j=1:size(phi,2)
        for k=1:size(phi,3)
            if(OPVolSkel(i,j,k)==1)
                mytree_x = [mytree_x;i];
                mytree_y = [mytree_y;j];
                mytree_z = [mytree_z;k];
            end
        end
    end
end

%% Compute Accuracy
xy_threshold = 3.94;
z_threshold = 5;

target = [tree.Y,tree.X,tree.Z];
source = [mytree_x,mytree_y,mytree_z];

[~,Idist] = pdist2(source,target,'euclidean','Smallest',1);

xycount = 0;
zcount = 0;
final_treex = [];
final_treey = [];
final_treez = [];
for points = 1:length(tree.X)
    xy1 = [mytree_x(Idist(points)),mytree_y(Idist(points))];
    xy2 = [tree.Y(points),tree.X(points)];
    zz1 = mytree_z(Idist(points));
    zz2 = tree.Z(points);
    Dxy = pdist2(xy1,xy2,'euclidean');
    Dz = pdist2(zz1,zz2,'euclidean');
    if(Dxy<=xy_threshold)
        xycount = xycount+1;
    end
    
    if(Dz<=z_threshold)
        zcount = zcount+1;
    end
    
    
    if(Dxy<=xy_threshold && Dz<=z_threshold)
        final_treex = [final_treex;mytree_x(Idist(points))];
        final_treey = [final_treey;mytree_y(Idist(points))];
        final_treez = [final_treez;mytree_z(Idist(points))];
    end
end

xy_Accuracy = xycount/length(tree.X)*100;
z_Accuracy = zcount/length(tree.X)*100;

end
