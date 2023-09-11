function a = vecangle360(v1,v2,n)
%% calculate angle between two vectors with a shared normal vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
% v1: vector 1; v2: vector 2; n: normal vector
%
% Output
% a: angle in degrees
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = cross(v1,v2);
c = sign(dot(x,n))*norm(x);
a = atan2d(c,dot(v1,v2));
end