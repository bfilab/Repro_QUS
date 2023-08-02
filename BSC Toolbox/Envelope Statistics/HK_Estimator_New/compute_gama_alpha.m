function [gama , alpha] = compute_gama_alpha (X,U)

alpha0=1/(X-1);
alpha1=0;

format long 
if X<=1  %% fully dev
    alpha2=0;
    alpha2=alpha2+1;
    nbpas=0;
    gama=compute_gama(alpha2,X);
    %     while (UHK_loop(gama,alpha2)>U) && (alpha2<59.5)
    while (UHK_loop(gama,alpha2)>U)
%         UHK_loop(gama,alpha2)-U;
        nbpas=nbpas+1;
%         disp (['number of initialization steps in the simultaneous alpha and gamma computation ',num2str(nbpas)])
        alpha2=alpha2+1
        gama=compute_gama(alpha2,X)
    end
else %% partically dev
    alpha2=min(alpha0,80);
%     alpha2=alpha0;
end


tol=1e-2;

nbpas2=0;
while abs(alpha1-alpha2)>tol
    
    alpha=(alpha1+alpha2)/2;
    nbpas2=nbpas2+1;
%     disp (['number of steps in the simultaneous alpha and gamma computation ',num2str(nbpas2)])
    gama=compute_gama(alpha,X);
    if UHK_loop(gama,alpha)<=U
        alpha2=alpha;
    else
        alpha1=alpha;
    end
    
end






end