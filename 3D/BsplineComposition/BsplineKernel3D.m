function [ kernel ] = BsplineKernel3D( )
%source code: https://github.com/stellaccl/cdmffd-image-registration
%paper: Chan, C. L., Anitescu, C., Zhang, Y., & Rabczuk, T. (2017). Two and three dimensional image registration based on B-spline composition and level sets. 
%Communications in Computational Physics, 21(2), 600-622.

%This function creates a 4x4x4 Bspline kernel 

%--Output Variable:
%kernel: 4x4x4 Bspline kernel 

b = [1/24, 11/24, 11/24, 1/24];

kernel=zeros(4,4,4);

for i=1:4
    for j=1:4
        for k=1:4
            kernel(j,i,k)=b(i)*b(j)*b(k);
        end
    end
end
            
end
