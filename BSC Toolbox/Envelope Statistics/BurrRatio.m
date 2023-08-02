function [Ratio] = BurrRatio(b)
Ratio = (b-2)*pi*gamma(b-1.5)^2/(4*gamma(b-1)^2);
end
