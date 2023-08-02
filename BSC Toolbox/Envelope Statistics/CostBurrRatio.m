function [E] = CostBurrRatio(b,M)
E = abs(BurrRatio(b)-M).^2;

end