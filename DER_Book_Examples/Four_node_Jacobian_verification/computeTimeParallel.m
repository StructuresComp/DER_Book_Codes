function  [d1, d2] = computeTimeParallel(d1_old, x0, x)

ne = (length(x) + 1) / 4 - 1;
tangent = zeros(ne, 3);
for c=1:ne
    dt = x0(4*c+1:4*c+3) - x0( 4*(c-1) + 1: 4*(c-1) + 3);
    dt = dt / norm(dt);
    tangent(c, :) = dt;
end

d1 = zeros(ne,3);
d2 = zeros(ne,3);
for c=1:ne
    edge = x(4*c+1:4*c+3) - x( 4*(c-1) + 1: 4*(c-1) + 3);
    t = edge/norm(edge);
    t = t';
    d1_l = parallel_transport(d1_old(c,:), tangent(c,:), t);
    d1_l = (d1_l - dot(d1_l, t) );
    d1_l = d1_l/ norm(d1_l);
    d1(c,:) = d1_l;
    d2(c,:) = cross(t, d1_l);    
end
