function [x, y] = find_nucleus_center(inputImg, kernelSize, kernelScale, threRatio, minSizeCell)

% Source code for soma detection taken from: https://github.com/kilho/NIA
% Paper: Kim, K. M., Son, K., & Palmore, G. T. R. (2015). Neuron image analyzer:
% Automated and accurate extraction of neuronal data from low quality images. Scientific reports, 5, 17062.
% Paper: Li, C., Xu, C., Gui, C., & Fox, M. D. (2010). Distance regularized level set evolution and its
% application to image segmentation. IEEE transactions on image processing, 19(12), 3243-3254.

%--Input Variable:
%inputImg: input image to be segmented
%kernelSize: kernel size for LoG filter
%kernelScale: readius of disk shaped filter
%threRatio: threshold for detecting cell
%minSizeCell: minimum size of the cell

%--Output Variable:
%x: x-coordinate of cell centers of detected neurons
%y: y-coordinate of cell centers of detected neurons

[height, width] = size(inputImg);
% make kernel and convolution
LoG = LoG_kernel(kernelSize,kernelScale);
LoGImg = conv2(inputImg, LoG, 'same');

% threshold
threshold = min(min(LoGImg));
threRatioImg = LoGImg/threshold;
threImg = (threRatioImg >= threRatio);

% getting rid of detection on boundary of image
threImg(1:minSizeCell,:) = 0;
threImg(height-minSizeCell:end,:) = 0;
threImg(:,1:minSizeCell) = 0;
threImg(:,width-minSizeCell:end,:) = 0;

% non maximum supression
for j = minSizeCell+1 : height-minSizeCell
    for k = minSizeCell+1 : width-minSizeCell
        if (threImg(j,k))
            tempImg = LoGImg(j-minSizeCell:j+minSizeCell,k-minSizeCell:k+minSizeCell);
            minValue = min(tempImg(:));
            if minValue ~= LoGImg(j,k)
                threImg(j,k) = 0;
            end
        end
    end
end

[y, x] = find(threImg == 1);

end
