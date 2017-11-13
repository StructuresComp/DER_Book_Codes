function [Fs, Js] = getFs(x)
%%
global EA nv ne refLen

%% compute tangent
tangentL = zeros(ne, 3); % tangent, computed locally
for c=1:ne
    dx = x(4*c+1:4*c+3) - x( 4*(c-1) + 1: 4*(c-1) + 3);
    tangentL(c,:) = dx / norm(dx);
end

%% Edge length
edgeLen = zeros(ne, 1);
for c=1:ne
    dx = x(4*c+1:4*c+3) - x(4*(c-1)+1:4*(c-1)+3);
    edgeLen(c) = norm(dx);
end
voronoiLen = zeros(nv, 1);
for c=1:nv
    if c==1
        voronoiLen(c) = 0.5 * edgeLen(c);
    elseif c==nv
        voronoiLen(c) = 0.5 * edgeLen(c-1);
    else
        voronoiLen(c) = 0.5 * (edgeLen(c-1) + edgeLen(c));
    end
end

%% Compute force
Fs = x * 0;
for c=1:ne
    epsX = edgeLen(c)/refLen(c) - 1;
    f = EA * tangentL(c,:) * epsX;
    Fs(4*(c-1)+1:4*(c-1)+3) = Fs(4*(c-1)+1:4*(c-1)+3) + f';
    Fs(4*c+1:4*c+3) = Fs(4*c+1:4*c+3) - f';
end

%% Jacobian
Id3 = eye(3);
Js = zeros(length(x), length(x));
for c=1:ne
    ci = 4*(c-1) + 1;
    cf = 4*(c-1) + 3;
    len = edgeLen(c);
    refLength = refLen(c);
    dx = x(4*c+1:4*c+3) - x(4*(c-1) + 1: 4*(c-1) + 3);
    
    M = EA * ( ...
        (1/refLength - 1/len) * Id3 + ...
        1/len * (dx*dx')/ (norm(dx))^2 ...
    ); %Note dx * dx' must be 3x3
    
    Js(ci:cf, ci:cf) = Js(ci:cf, ci:cf) - M;
    Js(ci+4:cf+4, ci+4:cf+4) = Js(ci+4:cf+4, ci+4:cf+4) - M;
    Js(ci+4:cf+4, ci:cf) = Js(ci+4:cf+4, ci:cf) + M;
    Js(ci:cf, ci+4:cf+4) = Js(ci:cf, ci+4:cf+4) + M;
end

end
