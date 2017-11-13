function [Ft, Jt] = getFt(x, refTwist)
%%
global GJ nv ne voronoiRefLen undeformedTwist

%% compute tangent
tangentL = zeros(ne,3); % tangent, computed locally
for c=1:ne
    dx = x(4*c+1:4*c+3) - x( 4*(c-1) + 1: 4*(c-1)+3);
    tangentL(c,1:3) = dx / norm(dx);
end

%% Edge length
edgeLen = zeros(ne, 1);
for c=1:ne
    dx = x(4*c+1:4*c+3) - x( 4*(c-1) + 1: 4*(c-1)+3);
    edgeLen(c) = norm(dx);
end

%% Curvature binormal
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

%% Compute twist
theta_f = x(4:4:end);
theta_e = [0; theta_f(1:end-1)];
deltam = theta_f - theta_e; % x(4:4:end); % m_i - m_{i-1}

%% Compute grad twist
gradTwist = zeros(nv, 11);
for c=2:ne
    if c==1
        norm_e = edgeLen(c);
        norm_f = edgeLen(c);
    elseif c==nv
        norm_e = edgeLen(c-1);
        norm_f = edgeLen(c-1);
    else
        norm_e = edgeLen(c-1);
        norm_f = edgeLen(c);
    end
    
    gradTwist(c, 1:3) = -0.5 / norm_e * kb(c,:);
    gradTwist(c, 9:11) = 0.5 / norm_f * kb(c,:);
    gradTwist(c, 5:7) = -(gradTwist(c, 1:3)+gradTwist(c, 9:11));
    gradTwist(c, 4) = -1;
    gradTwist(c, 8) = 1;
end

%%
Ft = x * 0;
for c=2:ne
    % Figure out twisting stiffness
    value = GJ / voronoiRefLen(c) * ( deltam(c) ...
        + refTwist(c) - undeformedTwist(c));
    
    ci = 4*(c-1) + 1 - 4;
    cig = 1; % initial location in gradTwist matrix
    if ci < 1
        cig = 1 - (ci-1); % initial location in gradKappa matrix
        ci = 1;
    end
    cf = 4*(c-1) + 1 + 6;
    cfg = 11;
    if cf > 3*nv + ne
        cfg = 11 - (cf - (3*nv+ne)); % final location in gradTwist matrix
        cf = 3*nv+ne;
    end
    
    force = - value * gradTwist(c,:);
    Ft( ci: cf) = Ft( ci: cf) + force(cig:cfg)';
end

%% Jacobian
DDtwist = zeros(nv, 11, 11);
for c=2:ne
    if c==1
        norm_e = edgeLen(c);
        norm_f = edgeLen(c);
    elseif c==nv
        norm_e = edgeLen(c-1);
        norm_f = edgeLen(c-1);
    else
        norm_e = edgeLen(c-1);
        norm_f = edgeLen(c);
    end
    
    norm2_e = norm_e^2;
    norm2_f = norm_f^2;
    
    if c==1
        te = tangentL(c,:);
        tf = te;
    elseif c==nv
        te = tangentL(c-1,:);
        tf = te;
    else
        te = tangentL(c-1,:);
        tf = tangentL(c,:);
    end
    
    kbLocal = kb(c,:);
    chi = 1.0 + dot(te, tf);
    tilde_t = (te + tf) / chi;
    
    D2mDe2 = -0.25 / norm2_e * ( kbLocal' * (te + tilde_t) ...
        + (te + tilde_t)' * kbLocal);
    D2mDf2 = -0.25 / norm2_f  * ( kbLocal' * (tf + tilde_t) ...
        + (tf + tilde_t)' * kbLocal );
    D2mDeDf = 0.5 / ( norm_e * norm_f ) * ( 2.0 / chi * crossMat( te ) - ...
        kbLocal' * tilde_t );
    D2mDfDe = D2mDeDf';
    
    DDtwist(c, 1:3,1:3) = D2mDe2;
    DDtwist(c, 1:3, 5:7) = -D2mDe2 + D2mDeDf;
    DDtwist(c, 5:7, 1:3) = -D2mDe2 + D2mDfDe;
    DDtwist(c, 5:7, 5:7) = D2mDe2 - ( D2mDeDf + D2mDfDe ) + D2mDf2;
    DDtwist(c, 1:3, 9:11) = -D2mDeDf;
    DDtwist(c, 9:11, 1:3) = -D2mDfDe;
    DDtwist(c, 9:11, 5:7) = D2mDfDe - D2mDf2;
    DDtwist(c, 5:7,9:11) = D2mDeDf - D2mDf2;
    DDtwist(c, 9:11,9:11) = D2mDf2;
end

%%
Jt = zeros(length(x), length(x));
for c=2:ne
    % Figure out twisting stiffness   
    ci = 4*(c-1) + 1 - 4;
    cig = 1; % initial location in gradKappa matrix
    if ci < 1
        cig = 1 - (ci-1); % initial location in gradKappa matrix
        ci = 1;
    end
    cf = 4*(c-1) + 1 + 6;
    cfg = 11;
    if cf > 3*nv+ne
        cfg = 11 - (cf - (3*nv+ne)); % final location in gradKappa matrix
        cf = 3*nv+ne;
    end
    
    milen = -1/voronoiRefLen(c);
    hessTwist = DDtwist(c,:,:);
    hessTwist = squeeze(hessTwist);
    gradTwistLocal = gradTwist(c,:);
    J = GJ*milen*( (deltam(c) + refTwist(c)-undeformedTwist(c)) * hessTwist ...
        + gradTwistLocal'*gradTwistLocal);
    
    J = J(cig:cfg,cig:cfg);
    Jt( ci: cf, ci: cf) = Jt( ci: cf, ci: cf) + J;
end

end
