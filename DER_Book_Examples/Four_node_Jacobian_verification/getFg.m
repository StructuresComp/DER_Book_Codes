function [Fg, Jg] = getFg(x)
%%
global m garr

%%
Fg =  m .* garr;
Jg = zeros(length(x), length(x));

end
