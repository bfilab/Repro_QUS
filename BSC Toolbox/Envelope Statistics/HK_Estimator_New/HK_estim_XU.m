function [ gama, alpha, epsilon, sigma ] = HK_estim_XU( data )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% F. Destrempes UMB 2014
% compute_gama -> algo figure 4
% compute_gama_alpha -> algo figure 5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data=data(:);                                                               % Get envelope amplitudes (A)
I=data.^2;                                                                  % Get envelope intensity (I)

I(I==0)=I(I==0)+eps;                                                        % Remove any 0 values

moy=mean(I);                                                                % Get mean intensity (I-bar)
X=mean(I.*log(I))/mean(I)-mean(log(I));                                     % X-statistics: Section 3
U=mean(log(I))-log(mean(I));                                                % U-statistics: Section 3

if X<=1                                                                     %% initialisation d'alpha en fonction de la structure du speckle 
    alpha_star=80;  % fully developped speckle 
else
    alpha0=1/(X-1); % partially developped speckle                          % Corollary 3.4
    alpha_star=min(alpha0,80);  
end

gama = compute_gama(alpha_star,X);                                          % gamma

if UHK_loop(gama,alpha_star) <= U                                           %% pourquoi? figure?
   [gama,alpha]=compute_gama_alpha(X,U);
else
    % gama=compute_gama(alpha_star,X); % redundant
    alpha=alpha_star;
end

epsilon=sqrt(moy*gama/(gama+alpha));                                        % Eqn. 3.5 (Destrempes 2014)
sigma=sqrt(moy/(2*(gama+alpha)));

end

