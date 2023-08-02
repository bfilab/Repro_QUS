function  [Ge_Gamma] = make_geneGamma_PDF(normalized_amp,bin)
Ge_Gamma.bin = bin;
%Ge_Gamma.para = ggammamomfit(normalized_amp);1
Ge_Gamma.para = ggammafit(normalized_amp);
Ge_Gamma.a = Ge_Gamma.para(1) ;
Ge_Gamma.c = Ge_Gamma.para(2) ;
Ge_Gamma.b = Ge_Gamma.para(3) ;
% Roi.Normalized.bin = bin;
Ge_Gamma.pdf = ggammafunc(bin,Ge_Gamma.a,Ge_Gamma.c,Ge_Gamma.b);
end

function paras = ggammafit( x )
% For computing the parameters of generalized gamma distribution
% input : x(Envelope)
% return: a(shape),c(shape adjust),b(scale) 
global A

A=abs(skewness(log(x)));

a=fsolve(@ggammamleequ,0.1); 
c=sqrt(psi(1,a)/var(log(x)));
b=mean(x)*gamma(a)/gamma(a+1/c);

paras=[ a,c,b ];
end

function f = ggammamleequ( a )
% For computing the shape parameter(a) of generalized gamma
% distribution
global A
f=psi(2,a)/((psi(1,a))^1.5)+A;
end

function [ y ] = ggammafunc( x, a, c, b )
% y = (c*b^(c*a)*x.^(c*a-1)/(gamma(a)).*exp(-b*x).^c); 
y = c.*x.^(a*c-1)./(gamma(a)*b^(a*c)).*exp(-(x./b).^c);
end