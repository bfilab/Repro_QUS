function  X  = XHK( gama,alpha )
% computation of the X statistics for given values of the gamma and alpha parameters
% Equation 3.4 (Destrempes 2014)

term0=(1+2*alpha)/(gama+alpha)-2*gama^(alpha/2+1/2)/((gama+alpha)*gamma(alpha))*besselk(alpha+1,2*sqrt(gama));

% term1=hypergeom ([1,1],[2,2,1-alpha],gama);
% term2=-pi/sin(pi*alpha)*gama^(alpha-1)*hypergeom(alpha,[1+alpha,1+alpha],gama)/(gamma(alpha)*gamma(alpha+1));
% term3=pi/sin(pi*(alpha+1))*gama^alpha*hypergeom(1+alpha,[2+alpha,2+alpha],gama)/(gamma(alpha)*gamma(alpha+2)*(1+alpha));
% term4=-alpha/(alpha-1)*hypergeom ([1,1],[2,2,2-alpha],gama);
% if alpha<1
%     term4=pi/sin(pi*(1-alpha))*alpha/(gamma(alpha)*gamma(2-alpha))*hypergeom ([1,1],[2,2,2-alpha],gama);
% end

% % hypergeom is part of the symbolic math toolbox. I don't have that.
% % Replacing with genHyper that I found on the file exchange:
% % https://www.mathworks.com/matlabcentral/fileexchange/5616-generalized-hypergeometric-function
% term1=genHyper ([1,1],[2,2,1-alpha],gama);
% term2=-pi/sin(pi*alpha)*gama^(alpha-1)*genHyper(alpha,[1+alpha,1+alpha],gama)/(gamma(alpha)*gamma(alpha+1));
% term3=pi/sin(pi*(alpha+1))*gama^alpha*genHyper(1+alpha,[2+alpha,2+alpha],gama)/(gamma(alpha)*gamma(alpha+2)*(1+alpha));
% term4=-alpha/(alpha-1)*genHyper ([1,1],[2,2,2-alpha],gama);
% if alpha<1
%     term4=pi/sin(pi*(1-alpha))*alpha/(gamma(alpha)*gamma(2-alpha))*genHyper ([1,1],[2,2,2-alpha],gama);
% end

term1=simplif_hyper_term1( alpha,gama );
term2=-pi/sin(pi*alpha)*gama^(alpha-1)*simplif_hyper_term2( alpha,gama )/(gamma(alpha)*gamma(alpha+1));
term3=pi/sin(pi*(alpha+1))*gama^alpha*simplif_hyper_term3( alpha,gama )/(gamma(alpha)*gamma(alpha+2)*(1+alpha));
term4=-alpha/(alpha-1)*simplif_hyper_term4( alpha,gama );

if alpha<1
    term4=pi/sin(pi*(1-alpha))*alpha/(gamma(alpha)*gamma(2-alpha))*simplif_hyper_term4( alpha,gama );
end

X = term0+gama/(gama+alpha)*(term1+term2+term3+term4);

end