function residual = helixDistance(A)
global currentX currentY currentZ

M1 = [A(1) A(2) A(3)]';
n1 = A(4);
n2 = A(5);
n3_comp = n1^2 + n2^2;
if (n3_comp >= 1.0)
    n3 = 0;
else
    n3 = sqrt(1 - n3_comp);
end

s = [n1 n2 n3]';
distance = currentX*0;
for k = 1:length(currentX)
    M0 = [currentX(k) currentY(k) currentZ(k)]';
    M0M1 = M1 - M0;
    distance(k) = norm(cross(M0M1, s)) / norm(s);
end

% distance

residual = sum((mean(distance) - distance).^2);

end
