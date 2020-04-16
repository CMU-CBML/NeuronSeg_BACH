function coef = img2coef2D(nx, ny, I)

%source code: https://github.com/stellaccl/cdmffd-image-registration
%paper: Chan, C. L., Anitescu, C., Zhang, Y., & Rabczuk, T. (2017). Two and three dimensional image registration based on B-spline composition and level sets. 
%Communications in Computational Physics, 21(2), 600-622.

%img2coe2D  Compute coeffident for B-spline level set representation of image I

%--Input Variable:     
%nx: number of row of Image I 
%ny: number of column of Image I
%I: input image

%--Output Variable
%coef: coefficient of B-spline level set function

%degree of level set function
p=3;
q=3;

I=reshape(I,[ny nx]);
coef = zeros(ny+q, nx+p);

%filter in the x-direction
parfor i=1:ny
    coef(i,:) = bspline1dfilt(I(i,:));    
end

%filter in the y-direction
for i=1:nx+3
    coef(:,i) = bspline1dfilt(coef(1:end-3,i)');
end
