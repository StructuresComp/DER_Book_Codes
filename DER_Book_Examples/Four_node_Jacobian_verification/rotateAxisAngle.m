function vNew = rotateAxisAngle( v, z, theta )

if (theta == 0) 
    vNew = v;
else
    c = cos(theta);
    s = sin(theta);
    vNew = c*v + s*cross(z,v) + dot(z,v) * (1.0-c) * z;
end
end
