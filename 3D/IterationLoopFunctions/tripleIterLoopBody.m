function [ pxx, pyy, pzz ] = tripleIterLoopBody(i, sizeImage, Pixel, Jm, ACP  )

%This function efficiently computes the value of deformed position of the
%pixels

%--Input Variable:
%i:index of Pixel location
%sizeImage: image size
%Pixel: struct variable storing the active cell index and support B-splines of each pixel coordinates
%Jm: Struct variables storing support B-splines over each active cell
%ACP: active control points

%--Output Variable:
%pxx:deformed positions in x direction
%pyy:deformed positions in y direction
%pzz:deformed positions in z direction

pxx = zeros(sizeImage(1,1),sizeImage(1,2));
pyy = zeros(sizeImage(1,1),sizeImage(1,2));
pzz = zeros(sizeImage(1,1),sizeImage(1,2));
for j = 1:sizeImage(1,2)
    for k = 1:sizeImage(1,1)
        px = sizeImage(1,1)*sizeImage(1,2)*(i-1)+sizeImage(1,1)*(j-1)+k;
        ac_ind = Pixel(px,1).active_cell;
        supp = Pixel(px,1).phi;
        SB = Jm(ac_ind,1).nzsplines;
        
        pts = ACP(SB,1:3);
        FXX = pts'*supp;
        pxx(k,j) = FXX(1,1);
        pyy(k,j) = FXX(2,1);
        pzz(k,j) = FXX(3,1);
    end
end
end

