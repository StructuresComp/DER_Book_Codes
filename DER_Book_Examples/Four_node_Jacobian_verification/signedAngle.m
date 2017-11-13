function angle = signedAngle( u, v, n )
w = cross(u,v);
angle = atan2( norm(w), dot(u,v) );
if (dot(n,w) < 0) 
    angle = -angle;
end
end

% C++ code
% inline Scalar signedAngle(const Vec3d& u, const Vec3d& v, const Vec3d& n)
% {
%   Vec3d w = u.cross(v);
%   Scalar angle = atan2(w.norm(), u.dot(v));
%   if (n.dot(w) < 0) return -angle;
%   return angle;
% }