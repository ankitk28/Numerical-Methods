function [ p, pp, ppp ] = pleg(x, c)
format long
if(numel(c) == 0)
    p(1:numel(x)) = 1;
    pp = zeros(numel(x),1);
    ppp = zeros(numel(x),1);
elseif(numel(c) == 1)
    p = x;
    pp(1:numel(x)) = 1;
    ppp = zeros(numel(x),1);
else
    [p__, pp__, ppp__] =  pleg(x, c(1:numel(c)-2));
    [p_, pp_, ppp_] =  pleg(x, c(1:numel(c)-1));
    p = x.*p_ - c(numel(c)).*p__;
    pp = p_ + x.*pp_ - c(numel(c)).*pp__;
    ppp = 2*pp_ + x.*ppp_ - c(numel(c)).*ppp__;
end

end