function [F, J] = objfun(xUncons)
%%
global m dt x0 uUncons nCons xCons
global ScaleSolver Fb Fs Ft Fg
global d1 refTwist

%%
x0Uncons = x0(nCons+1:end);
x = [xCons; xUncons];

%% Update directors to get reference twist
tangent = computeTangent(x);
[d1Iterate, d2Iterate] = computeTimeParallel(d1, x0, x);
refTwistIterate = getRefTwist(d1Iterate, tangent, refTwist); % Compute reference twist
theta = x(4:4:end);
[m1, m2] = computeMaterialDirectors(d1Iterate, d2Iterate, theta);

%%
mUncons = m(nCons+1:end);

[Fb, Jb] = getFb(x, m1, m2);
[Fg, Jg] = getFg(x);
[Fs, Js] = getFs(x);
[Ft, Jt] = getFt(x, refTwistIterate);

Forces = (Fb + Fs + Ft + Fg);
Forces = Forces(nCons+1:end);
% F = mUncons.*uUncons/dt - Forces;
% xNew = x0Uncons - F * dt^2 ./ mUncons;

F = mUncons .* (xUncons - x0Uncons)/dt^2 - mUncons.*uUncons/dt - Forces;
%     + visc * abs(xUncons - x0Uncons)/dt;

mMat = diag(mUncons);
Jforces = Jb + Js + Jt + Jg;
Jforces = Jforces(nCons+1:end, nCons+1:end);
J = mMat/dt^2 - Jforces;

F = F*ScaleSolver;
J = J*ScaleSolver;
end
