function H = Heaviside(phi,epsilon)
%This function computes regularized Heaviside function
% Input: phi: level set function
% epsilon: width of level set function
%Ouput:regularized Heaviside function

H = 0.5*(1+(2/pi)*atan(phi./epsilon));
end
