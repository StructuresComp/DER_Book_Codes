function [Fb, Jb] = getFb(x, m1, m2)
%%
global EI nv ne voronoiRefLen kappaBar

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

%% Computer Kappa
kappa = zeros(nv, 2);
for c=2:ne
    
    m1e = m1(c-1,:);
    m2e = m2(c-1,:);
    m1f = m1(c,:);
    m2f = m2(c,:);
    
    kappa1 = 0.5 * dot( kb(c,:), m2e + m2f); % CHECKED
    kappa2 = -0.5 * dot( kb(c,:), m1e + m1f); % CHECKED
    kappa(c,1) = kappa1;
    kappa(c,2) = kappa2;
end

%% Compute grad kappa
gradKappa = zeros(nv, 11, 2);
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
    
    if c==1
        d1e = m1(1,:);
        d2e = m2(1,:);
        d1f = m1(1,:);
        d2f = m2(1,:);
    elseif c==nv
        d1e = m1(ne,:);
        d2e = m2(ne,:);
        d1f = m1(ne,:);
        d2f = m2(ne,:);
    else
        d1e = m1(c-1,:);
        d2e = m2(c-1,:);
        d1f = m1(c,:);
        d2f = m2(c,:);
    end
    chi = 1.0 + dot(te, tf);
    tilde_t = (te + tf) / chi;
    tilde_d1 = (d1e + d1f) / chi;
    tilde_d2 = (d2e + d2f) / chi;
    
    kappa1 = kappa(c, 1);
    kappa2 = kappa(c, 2);
    
    Dkappa1De = 1.0 / norm_e * (-kappa1 * tilde_t + cross(tf,tilde_d2));
    Dkappa1Df = 1.0 / norm_f * (-kappa1 * tilde_t - cross(te,tilde_d2));
    
    Dkappa2De = 1.0 / norm_e * (-kappa2 * tilde_t - cross(tf,tilde_d1));
    Dkappa2Df = 1.0 / norm_f * (-kappa2 * tilde_t + cross(te,tilde_d1));
    
    gradKappa(c, 1:3, 1) = -Dkappa1De;
    gradKappa(c, 5:7, 1) = Dkappa1De - Dkappa1Df;
    gradKappa(c, 9:11, 1) = Dkappa1Df;
    
    gradKappa(c, 1:3, 2) = -Dkappa2De;
    gradKappa(c, 5:7, 2) = Dkappa2De - Dkappa2Df;
    gradKappa(c, 9:11, 2) = Dkappa2Df;
    
    kbLocal = kb(c,:);
    gradKappa(c, 4, 1) = -0.5 * dot(kbLocal, d1e);
    gradKappa(c, 8, 1) = -0.5 * dot(kbLocal, d1f);
    gradKappa(c, 4, 2) = -0.5 * dot(kbLocal, d2e);
    gradKappa(c, 8, 2) = -0.5 * dot(kbLocal, d2f);
end
%%
Fb = x * 0;
for c=2:ne
    % Figure out bending stiffness
    EIMat = [ EI 0; ...
        0 EI];
    
    % Now the real business
    ci = 4*(c-1) + 1 - 4;
    cig = 1; % initial location in gradKappa matrix
    if ci < 1
        cig = 1 - (ci-1); % initial location in gradKappa matrix
        ci = 1;
    end
    cf = 4*(c-1) + 1 + 6;
    cfg = 11;
    if cf > 3*nv + ne
        cfg = 11 - (cf - (3*nv+ne)); % final location in gradKappa matrix
        cf = 3*nv+ne;
    end
    
    relevantPart = gradKappa(c, cig:cfg, :); %units of 1/length
    relevantPart = squeeze(relevantPart); % converts a 3D into 2D matrix
    kappaL = kappa(c,:) - kappaBar(c,:);
    f = relevantPart * EIMat * kappaL';
    Fb( ci: cf) = Fb( ci: cf) - f/voronoiRefLen(c);
end

%% Jacobian
DDkappa1All = zeros(nv, 11, 11);
DDkappa2All = zeros(nv, 11, 11);
DDkappa1 = zeros(11, 11);
DDkappa2 = zeros(11, 11);

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
    
    if c==1
        d1e = m1(1,:);
        d2e = m2(1,:);
        d1f = m1(1,:);
        d2f = m2(1,:);
    elseif c==nv
        d1e = m1(ne,:);
        d2e = m2(ne,:);
        d1f = m1(ne,:);
        d2f = m2(ne,:);
    else
        d1e = m1(c-1,:);
        d2e = m2(c-1,:);
        d1f = m1(c,:);
        d2f = m2(c,:);
    end
    
    chi = 1.0 + dot(te, tf);
    tilde_t = (te + tf) / chi;
    tilde_d1 = (d1e + d1f) / chi;
    tilde_d2 = (d2e + d2f) / chi;
    
    kappa1 = kappa(c, 1);
    kappa2 = kappa(c, 2);
    
    kbLocal = kb(c,:);
    
    tt_o_tt = tilde_t' * tilde_t; % must be 3x3. tilde_t is 1x3
    tmp = cross(tf, tilde_d2);
    tf_c_d2t_o_tt = tmp' * tilde_t; % must be 3x3
    tt_o_tf_c_d2t = tf_c_d2t_o_tt'; % must be 3x3
    kb_o_d2e = kbLocal' * d2e; % must be 3x3
    d2e_o_kb = kb_o_d2e'; % must be 3x3
    
    Id3 = eye(3);
    D2kappa1De2 ...
        = 1.0 / norm2_e * (2 * kappa1 * tt_o_tt - tf_c_d2t_o_tt - tt_o_tf_c_d2t) ...
        - kappa1 / (chi * norm2_e) * (Id3 - te'*te) ...
        + 1.0 / (4.0 * norm2_e) * (kb_o_d2e + d2e_o_kb);
    
    tmp = cross(te, tilde_d2);
    te_c_d2t_o_tt = tmp' * tilde_t;
    tt_o_te_c_d2t = te_c_d2t_o_tt';
    kb_o_d2f = kbLocal' * d2f;
    d2f_o_kb = kb_o_d2f';
    
    D2kappa1Df2 ...
        = 1.0 / norm2_f * (2 * kappa1 * tt_o_tt + te_c_d2t_o_tt + tt_o_te_c_d2t) ...
        - kappa1 / (chi * norm2_f) * (Id3 - tf'*tf) ...
        + 1.0 / (4.0 * norm2_f) * (kb_o_d2f + d2f_o_kb);
    
    D2kappa1DeDf ...
        = -kappa1/(chi * norm_e * norm_f) * (Id3 + te'*tf) ...
        + 1.0 / (norm_e*norm_f) * (2 * kappa1 * tt_o_tt - tf_c_d2t_o_tt + ...
        tt_o_te_c_d2t - crossMat(tilde_d2));
    D2kappa1DfDe = D2kappa1DeDf';
    
    tmp = cross(tf, tilde_d1);
    tf_c_d1t_o_tt = tmp'*tilde_t; % must be 3x3
    tt_o_tf_c_d1t = tf_c_d1t_o_tt'; % must be 3x3
    kb_o_d1e = kbLocal'*d1e; % must be 3x3
    d1e_o_kb = kb_o_d1e'; % must be 3x3
    
    D2kappa2De2 ...
        = 1.0 / norm2_e * (2 * kappa2 * tt_o_tt + tf_c_d1t_o_tt + tt_o_tf_c_d1t) ...
        - kappa2 / (chi * norm2_e) * (Id3 - te'*te) ...
        - 1.0 / (4.0 * norm2_e) * (kb_o_d1e + d1e_o_kb);
    
    tmp = cross(te, tilde_d1);
    te_c_d1t_o_tt = tmp'*tilde_t; % must be 3x3
    tt_o_te_c_d1t = te_c_d1t_o_tt'; % must be 3x3
    kb_o_d1f = kbLocal'*d1f; % must be 3x3
    d1f_o_kb =  kb_o_d1f'; % must be 3x3
    
    D2kappa2Df2 ...
        = 1.0 / norm2_f * (2 * kappa2 * tt_o_tt - te_c_d1t_o_tt - tt_o_te_c_d1t) ...
        - kappa2 / (chi * norm2_f) * (Id3 - tf'*tf) ...
        - 1.0 / (4.0 * norm2_f) * (kb_o_d1f + d1f_o_kb); % must be 3x3
    
    D2kappa2DeDf ...
        = -kappa2/(chi * norm_e * norm_f) * (Id3 + te'*tf) ...
        + 1.0 / (norm_e*norm_f) * (2 * kappa2 * tt_o_tt + tf_c_d1t_o_tt - tt_o_te_c_d1t + crossMat(tilde_d1));
    % must be 3x3
    D2kappa2DfDe = D2kappa2DeDf'; % must be 3x3
    
    D2kappa1Dthetae2 = -0.5 * dot(kbLocal, d2e);
    D2kappa1Dthetaf2 = -0.5 * dot(kbLocal, d2f);
    D2kappa2Dthetae2 =  0.5 * dot(kbLocal, d1e);
    D2kappa2Dthetaf2 =  0.5 * dot(kbLocal, d1f);
    
    D2kappa1DeDthetae ...
        = 1.0 / norm_e * (0.5 * dot(kbLocal, d1e) * tilde_t - 1.0 / chi * cross(tf, d1e));
    D2kappa1DeDthetaf ...
        = 1.0 / norm_e * (0.5 * dot(kbLocal, d1f) * tilde_t - 1.0 / chi * cross(tf, d1f));
    D2kappa1DfDthetae ...
        = 1.0 / norm_f * (0.5 * dot(kbLocal, d1e) * tilde_t + 1.0 / chi * cross(te, d1e));
    D2kappa1DfDthetaf ...
        = 1.0 / norm_f * (0.5 * dot(kbLocal, d1f) * tilde_t + 1.0 / chi * cross(te, d1f));
    D2kappa2DeDthetae ...
        = 1.0 / norm_e * (0.5 * dot(kbLocal, d2e) * tilde_t - 1.0 / chi * cross(tf, d2e));
    D2kappa2DeDthetaf ...
        = 1.0 / norm_e * (0.5 * dot(kbLocal, d2f) * tilde_t - 1.0 / chi * cross(tf, d2f));
    D2kappa2DfDthetae ...
        = 1.0 / norm_f * (0.5 * dot(kbLocal, d2e) * tilde_t + 1.0 / chi * cross(te, d2e));
    D2kappa2DfDthetaf ...
        = 1.0 / norm_f * (0.5 * dot(kbLocal, d2f) * tilde_t + 1.0 / chi * cross(te, d2f));
    
    % Curvature terms
    DDkappa1(1:3, 1:3)  =   D2kappa1De2;
    DDkappa1(1:3, 5:7)  = - D2kappa1De2 + D2kappa1DeDf;
    DDkappa1(1:3, 9:11) =               - D2kappa1DeDf;
    DDkappa1(5:7, 1:3)  = - D2kappa1De2                + D2kappa1DfDe;
    DDkappa1(5:7, 5:7)  =   D2kappa1De2 - D2kappa1DeDf - D2kappa1DfDe + D2kappa1Df2;
    DDkappa1(5:7, 9:11) =                 D2kappa1DeDf                - D2kappa1Df2;
    DDkappa1(9:11, 1:3)  =                              - D2kappa1DfDe;
    DDkappa1(9:11, 5:7)  =                                D2kappa1DfDe - D2kappa1Df2;
    DDkappa1(9:11, 9:11) =                                               D2kappa1Df2;
    
    % Twist terms
    DDkappa1(4, 4)     =   D2kappa1Dthetae2;
    DDkappa1(8, 8)     =   D2kappa1Dthetaf2;
    
    % Curvature-twist coupled terms
    DDkappa1(1:3, 4)   = - D2kappa1DeDthetae;
    DDkappa1(5:7, 4)   =   D2kappa1DeDthetae - D2kappa1DfDthetae;
    DDkappa1(9:11,4)   =                       D2kappa1DfDthetae;
    DDkappa1(4, 1:3)   =   transpose(DDkappa1(1:3, 4));
    DDkappa1(4, 5:7)   =   transpose(DDkappa1(5:7, 4));
    DDkappa1(4, 9:11)  =   transpose(DDkappa1(9:11,4));
    
    % Curvature-twist coupled terms
    DDkappa1(1:3, 8)   = - D2kappa1DeDthetaf;
    DDkappa1(5:7, 8)   =   D2kappa1DeDthetaf - D2kappa1DfDthetaf;
    DDkappa1(9:11, 8)  =                       D2kappa1DfDthetaf;
    DDkappa1(8, 1:3)   =   transpose(DDkappa1(1:3, 8));
    DDkappa1(8, 5:7)   =   transpose(DDkappa1(5:7, 8));
    DDkappa1(8, 9:11)  =   transpose(DDkappa1(9:11,8));
    
    % Curvature terms
    DDkappa2(1:3, 1:3) =   D2kappa2De2;
    DDkappa2(1:3, 5:7) = - D2kappa2De2 + D2kappa2DeDf;
    DDkappa2(1:3, 9:11) =               - D2kappa2DeDf;
    DDkappa2(5:7, 1:3) = - D2kappa2De2                + D2kappa2DfDe;
    DDkappa2(5:7, 5:7) =   D2kappa2De2 - D2kappa2DeDf - D2kappa2DfDe + D2kappa2Df2;
    DDkappa2(5:7, 9:11)=                 D2kappa2DeDf                - D2kappa2Df2;
    DDkappa2(9:11, 1:3)=                              - D2kappa2DfDe;
    DDkappa2(9:11, 5:7)=                                D2kappa2DfDe - D2kappa2Df2;
    DDkappa2(9:11, 9:11)=                                               D2kappa2Df2;
    
    % Twist terms
    DDkappa2(4, 4)     = D2kappa2Dthetae2;
    DDkappa2(8, 8)     = D2kappa2Dthetaf2;
    
    % Curvature-twist coupled terms
    DDkappa2(1:3, 4)   = - D2kappa2DeDthetae;
    DDkappa2(5:7, 4)   =   D2kappa2DeDthetae - D2kappa2DfDthetae;
    DDkappa2(9:11,4)   =                       D2kappa2DfDthetae;
    DDkappa2(4, 1:3)   =   transpose(DDkappa2(1:3, 4));
    DDkappa2(4, 5:7)   =   transpose(DDkappa2(5:7, 4));
    DDkappa2(4, 9:11)  =   transpose(DDkappa2(9:11,4));
    
    % Curvature-twist coupled terms
    DDkappa2(1:3, 8)   = - D2kappa2DeDthetaf;
    DDkappa2(5:7, 8)   =   D2kappa2DeDthetaf - D2kappa2DfDthetaf;
    DDkappa2(9:11,8)   =                       D2kappa2DfDthetaf;
    DDkappa2(8, 1:3)   =   transpose(DDkappa2(1:3, 8));
    DDkappa2(8, 5:7)   =   transpose(DDkappa2(5:7, 8));
    DDkappa2(8,9:11)   =   transpose(DDkappa2(9:11,8));
    
    % Store
    DDkappa1All(c, :, :) = DDkappa1;
    DDkappa2All(c, :, :) = DDkappa2;
end

%%
Jb = zeros(length(x), length(x));
for c=2:ne
    % Figure out bending stiffness
    EIMat = [ EI 0; ...
        0 EI];

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
    
    len = voronoiRefLen(c);
    
    relevantPart = gradKappa(c, cig:cfg, :); %units of 1/length
    relevantPart = squeeze(relevantPart); % converts a 3D into 2D matrix
    df = -1.0 / len * relevantPart * EIMat * transpose(relevantPart);

    relevantPart1 = DDkappa1All(c, cig:cfg, cig:cfg);
    relevantPart2 = DDkappa2All(c, cig:cfg, cig:cfg);
    relevantPart1 = squeeze(relevantPart1);
    relevantPart2 = squeeze(relevantPart2);
    kappaL = kappa(c,:) - kappaBar(c,:);
    temp = -1.0 / len * kappaL * EIMat;
    df = df + temp(1) * relevantPart1 + temp(2) * relevantPart2;
    Jb( ci: cf, ci: cf) = Jb( ci: cf, ci: cf) + df;
end

end
