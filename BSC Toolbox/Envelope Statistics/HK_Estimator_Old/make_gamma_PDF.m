function [GA] = make_gamma_PDF(normalized_amp,bin)
GA.bin = bin;
GA.para = fitdist(normalized_amp,'gamma');
GA.a = GA.para.a ;
GA.b = GA.para.b ;
GA.pdf = pdf('Gamma',bin,GA.a,GA.b);
end

