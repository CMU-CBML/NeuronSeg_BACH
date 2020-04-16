function T1= Initial(knv1,knv2)

% This function intiializes the coefficient matrix before knot insertion

%--Input Variable:
%knv1: coarse level knot vector
%knv2: fine level knot vector

%--Output Variable:
%T1: initialized coefficient matrix


T1 = zeros(length(knv2)-1,length(knv1)-1);
for i = 1: length(knv2)-1
    for j = 1:length(knv1)-1
        if (knv2(i) >= knv1(j) && knv2(i) < knv1(j+1))
            T1(i,j) = 1;
        else
            T1(i,j) = 0;
        end
    end
end

end