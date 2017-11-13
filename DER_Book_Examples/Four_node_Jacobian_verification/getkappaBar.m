function kappaBar = getkappaBar( x, m1, m2 )
nv = (length(x)+1)/4;
ne = nv-1;

% compute tangent
tangentL = zeros(ne,3); % tangent, computed locally
for c=1:ne
    dx = x(4*c+1:4*c+3) - x( 4*(c-1) + 1: 4*(c-1)+3);
    tangentL(c,1:3) = dx / norm(dx);
end

% Curvature binormal
% Computes the curvature binormal given two unit vectors.
% t0 The first tangent
% t1 The second tangent
kb = zeros(nv, 3);
for c=1:nv
    if c==1 || c==nv
        kb(c,:) = [0 0 0];
    else
        t0 = tangentL(c-1,:);
        t1 = tangentL(c,:);
        kb(c,:) = 2.0 * cross(t0, t1) / (1.0 + dot(t0, t1));
    end
end

% Computer KappaBar
kappaBar = zeros(nv, 2);
for c=1:nv
    if c==1
        m1e = m1(1,:);
        m2e = m2(1,:);
        m1f = m1(1,:);
        m2f = m2(1,:);
    elseif c==nv
        m1e = m1(ne,:);
        m2e = m2(ne,:);
        m1f = m1(ne,:);
        m2f = m2(ne,:);
    else
        m1e = m1(c-1,:);
        m2e = m2(c-1,:);
        m1f = m1(c,:);
        m2f = m2(c,:);
    end
    kappa1 = 0.5 * dot( kb(c,:), m2e + m2f); % CHECKED
    kappa2 = -0.5 * dot( kb(c,:), m1e + m1f); % CHECKED
    kappaBar(c,1) = kappa1;
    kappaBar(c,2) = kappa2;
end

end

