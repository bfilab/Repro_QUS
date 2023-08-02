function  U  = UHK( gama,alpha )

% computation of the U statistics for given values of the gamma and alpha parameters
% by aproximation of the formula
%  U = - eulergamma - log(gama+alpha) + psi(alpha) ...
%      - gama^alpha * (gamma(-alpha)/(alpha*gamma(alpha))) * hypergeom (alpha,[1+alpha,1+alpha],gama) ...
%      + gama * (gamma(alpha-1)/gamma(alpha)) * hypergeom ([1,1],[2,2,2-alpha],gama) ;

eulergamma=0.5772156649;

term1=-eulergamma-log(gama+alpha)+psi(alpha);

% term2=pi/sin(pi*alpha)*gama^alpha*hypergeom(alpha,[1+alpha,1+alpha],gama)/(alpha*gamma(alpha)*gamma(alpha+1)); 
% term3=gama/(alpha-1)*hypergeom ([1,1],[2,2,2-alpha],gama);
% if alpha<1
%     term3=-pi/sin(pi*(1-alpha))*gama/(gamma(alpha)*gamma(2-alpha))*hypergeom ([1,1],[2,2,2-alpha],gama);
% end

term2=pi/sin(pi*alpha)*gama^alpha*simplif_hyper_term2( alpha,gama )/(alpha*gamma(alpha)*gamma(alpha+1)); 
term3=gama/(alpha-1)*simplif_hyper_term4( alpha,gama );

if alpha<1
    term3=-pi/sin(pi*(1-alpha))*gama/(gamma(alpha)*gamma(2-alpha))*simplif_hyper_term4( alpha,gama );
end

U = term1+term2+term3;


end

