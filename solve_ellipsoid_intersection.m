function t = solve_ellipsoid_intersection(prev_pos,vec,A)
%created 4/1/16
%JH
%find intersection for reflection in an ellipsoid
%of form r'Ar=1
%prev_pos is previous position
%vec is vector direction of movement
%A is matrix containing ellipse scalings 1/a^2,1/b^2,1/c^2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alpha1 = diag(A)'*vec.^2;
alpha2 = prev_pos'*A*vec;
alpha3 = diag(A)'*prev_pos.^2 - 1;

t = [(-alpha2+sqrt(alpha2^2-alpha1*alpha3))/alpha1;
    (-alpha2-sqrt(alpha2^2-alpha1*alpha3))/alpha1];

