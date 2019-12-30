function [f1,f2] = Local_Avr(I,H,K)
%source obtained from: http://kaihuazhang.net/J_papers/PR_10.rar
%Paper: K. Zhang, H. Song, and L. Zhang., Active Contours Driven by Local Image Fitting Energy., Pattern recognition, vol.43, issue 4, pp. 1199-1206, April 2010. 

%This function computes f1 and f2 constants in local image fitting term
%INPUT: I: input image
%H: Heaviside function
%K: Gaussian kernel
%OUTPUT: f1, f2: constants for local image fitting term

f1 = conv2(I.*H,K,'same');
c1 = conv2(H,K,'same');
f1 = f1./c1;
f2 = conv2(I.*(1-H),K,'same');
c2 = conv2(1-H,K,'same');
f2 = f2./c2;

end