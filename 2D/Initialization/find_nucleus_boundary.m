function [x, y, boundaryNecleus] = find_nucleus_boundary(inputImg, inputImg1, kernelSize, kernelScale, threRatio, minSizeCell, levelsetR)
%% necleus detection (blob detector, LOG + Ratio threshold + NMS + Level set) 

[height, width] = size(inputImg);
% make kernel and convolution
LoG = LoG_kernel(kernelSize,kernelScale);
LoGImg = conv2(inputImg, LoG, 'same');
% figure(2)
% imagesc(LoGImg)
% title('LoG image')
% colorbar
% axis off

%     figure(21)
%     [tempX tempY] = meshgrid(1:kernelSize, 1:kernelSize);
%     surf(tempX, tempY, LoG)

% threshold 
threshold = min(min(LoGImg));
threRatioImg = LoGImg/threshold;
threImg = (threRatioImg >= threRatio);

% saveThreImg = threImg;
%     figure(3)
%     imagesc(threImg);
%     title('threshold image')

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
            if minValue ~= LoGImg(j,k);
                threImg(j,k) = 0;
            end
        end
    end
end
%     figure(4)
%     imagesc(threImg);
%     title('threshold image with non maximum supression')

% center of cell
%tempIdx = find(threImg == 1);
%x = ceil((tempIdx-1)/height);
%y = mod(tempIdx-1,height)+1;
[y, x] = find(threImg == 1);
numNeuron = length(x);

% figure(3)
% imshow(saveThreImg)
% hold on
% plot(x, y, 'g.', 'MarkerSize', 20)
% hold off

%% necleus boundary detection
boundaryNecleus = 2*ones(height,width);
for j = 1:numNeuron
    centerPoint = [x(j) y(j)];
    tempBoundaryNecleus = levelset_nucleus(centerPoint, levelsetR , inputImg1, height, width);
    boundaryNecleus = boundaryNecleus.*tempBoundaryNecleus/2;
end
%     figure(100)
%     imagesc(inputImg,[0, 1]); axis off; axis equal; colormap(gray); 
%     hold on;
%     contour(boundaryNecleus, [0,0], 'r');
%     hold off;
%     pause