function Delta_h = Delta(phi,epsilon)

%This function computes dirac delta function
% Input: phi: level set function
% epsilon: width of level set function
%Ouput: dirac delta function

Delta_h = (epsilon/pi)./(epsilon^2+phi.^2);
end
