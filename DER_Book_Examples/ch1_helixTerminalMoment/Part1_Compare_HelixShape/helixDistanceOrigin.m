function residual = helixDistanceOrigin(A)

global N1

A_new = [A; N1];
residual = helixDistance(A_new);

end
