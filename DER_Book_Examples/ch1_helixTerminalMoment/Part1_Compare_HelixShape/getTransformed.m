function ptTransformed = getTransformed(M1, n, x, y, z)

zAxis = n;
zAxis = zAxis / norm(zAxis);

M0 = [x(1); y(1); z(1)];
yAxis = (M1 - M0) - dot(M1-M0, zAxis) * zAxis;
yAxis = yAxis / norm(yAxis);
xAxis = cross(yAxis, zAxis);
xAxis = xAxis / norm(xAxis);

CoordMatrix = [xAxis yAxis zAxis];

ptTransformed = zeros(length(x), 3);
for k=1:length(x)
    pt = [x(k); y(k); z(k)] - M1;
    ptTransformed(k,:) = CoordMatrix \ pt;
end

end
