function [m1, m2] = computeMaterialDirectors(d1, d2, theta0)
ne = length(theta0);
m1 = zeros(ne,3);
m2 = m1;
for c=1:ne
    cs = cos(theta0(c));
    ss = sin(theta0(c));
    d1_l = d1(c,:);
    d2_l = d2(c,:);
    m1(c,:) = cs * d1_l + ss * d2_l;
    m2(c,:) = - ss * d1_l + cs * d2_l;
end
end
