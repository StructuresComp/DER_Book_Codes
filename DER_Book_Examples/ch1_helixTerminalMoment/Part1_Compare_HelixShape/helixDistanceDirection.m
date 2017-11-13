function residual = helixDistanceDirection(A)

global M1

A_new = [M1; A(1); A(2)];
residual = helixDistance(A_new);

end
