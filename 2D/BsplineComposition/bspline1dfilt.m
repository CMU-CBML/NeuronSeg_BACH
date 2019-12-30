function coef=bspline1dfilt(I)
%source code: https://github.com/stellaccl/cdmffd-image-registration
%paper: Chan, C. L., Anitescu, C., Zhang, Y., & Rabczuk, T. (2017). Two and three dimensional image registration based on B-spline composition and level sets. 
%Communications in Computational Physics, 21(2), 600-622.

%This function computes the 1D coefficients 

% Input: image
% coef: coefficients

N=length(I);
p=3;
num_basis=N+p;
M = [ 1/24, 11/24, 11/24, 1/24];

numerator=zeros(1,num_basis);
denominator=ones(1,num_basis);
denominator(1:p)=[1/24 0.5 23/24];
denominator(num_basis-2:num_basis)=[23/24 0.5 1/24];

for i=1:N
    numerator_temp=zeros(1,num_basis);
    numerator_temp(i:i+p)=M*I(i);
    numerator=numerator+numerator_temp;
end

coef=numerator./denominator;



