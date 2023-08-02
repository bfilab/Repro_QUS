function [ y ] = make_K_PDF( amp, bin )

[para(1), para(2)] = kfit(amp);
% if para(1) >= 0
    y.pdf = Kfunc(bin, para(1), para(2));
    y.bin = bin;
    y.a = double(para(1));
    y.b = double(real(para(2)));
% else
%     y.pdf = NaN;
%     y.bin = bin;
%     y.a = NaN;
%     y.b = NaN;
end

function [ a,b ] = kfit( x )
% For computing the parameters of K distribution
% input : Envelope
% return: a(shape),b(scale) 

% a=2/(mean(x.^4)/(mean(x.^2))^2-2); %wrong statement 2015.11.3 Masaaki
% Omura
E = mean(x.^4)/(mean(x.^2)^2);
a=2/(E-2);
b=sqrt(4*a/mean(x.^2));

end

function y=Kfunc(x,a,b)
% input : x(Envelope),a(shape),b(scale)
% output: p|x 
% G=gamma(a);
K=besselk(a-1,b*x);
% y=(b*x/2).^a*(2*b/G).*K; wrong statement(20151030 Masaaki OMURA)
y = 2*(x./2).^a*(b^(a+1)/gamma(a)).*K;
if isnan(y(1))
    y(1)=0;
end

end
