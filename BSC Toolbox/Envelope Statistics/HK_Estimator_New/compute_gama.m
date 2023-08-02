function gama = compute_gama(alpha,X)
% Compute gamma, given alpha.
% Figure 4 (Destrempes 2014)

gama1=0;
gama2=0;

gama2=gama2+1;

nbpas=0;
while XHK_loop(gama2,alpha)>X
    nbpas=nbpas+1;
    disp (['number of initialization steps in the gamma computation ',num2str(nbpas)])
    gama2=gama2+1;
end


tol=1e-4;
gama=(gama1+gama2)/2;

while abs(gama1-gama2)>tol
    
    gama=(gama1+gama2)/2;
    
    if XHK_loop(gama,alpha)<=X
        gama2=gama;
    else
        gama1=gama;
    end
    
end



end