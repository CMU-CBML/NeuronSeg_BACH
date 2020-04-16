function Tnew = Tmatrix(knv1,knv2,p)

%Function to calculate new control points using Oslo's Knot Insertion
%algorithm

%--Input Variables:
%knv1: coarse level knot vector
%knv2: fine level knot vector
%p: B-spline degree

%--Output Variable:
%Tnew: new coefficient matrix after knot insertion

T1 = Initial(knv1,knv2);
%Tnew = T1;
for i = 1:p
    T2 = KnotInsert(knv1,knv2,T1,i);
    %Tnew = Tnew(1:size(T2,1),1:size(T2,2));
    T1 = T2;
end
Tnew= T1;


end
