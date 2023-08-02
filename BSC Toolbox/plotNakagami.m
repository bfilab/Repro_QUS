function [nak_x,nak_y] = plotNakagami(model_nak,hi)
%PLOTNAKAGAMI Return points to plot the Nakagami pdf for given
%parameters.

% 11/15/2021 (THL): Created

nak_x = hi.BinEdges(1:(end-1))+hi.BinWidth/2;
nak_y = pdf(model_nak, nak_x)*hi.BinWidth;
            
end

