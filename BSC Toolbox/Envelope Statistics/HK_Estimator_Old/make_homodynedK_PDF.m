function  [HomodynedK] = make_homodynedK_PDF(normalized_amp,bin,eq)
% 'Equation of Hruska' is referenced by Hruska, Oelze et al.(2009)
% 'Equation of Dutt' is referenced by Dutt, Greenleaf et al.(1994)  

global A2 A4 A6
switch eq
    case 'Hruska' 
        [HomodynedK.k,HomodynedK.mu,HomodynedK.err] = estimator_RSK(normalized_amp,0);
        % HomodynedK.sigma = 1/(HomodynedK.k^2+2);        % s^2+2*É–^2=1 Ç∆Å@k=s/É– Ç…ÇÊÇÈ
        % HomodynedK.s = HomodynedK.k*sqrt(HomodynedK.sigma);     % HK_Sigma=É–^2
        % HomodynedK.pdf = homok_func(bin,HomodynedK.mu,HomodynedK.sigma,HomodynedK.s);
        % ver2
        % HomodynedK.sigma = sqrt(mean(normalized_amp)/(2*size(normalized_amp,1)));
        % HomodynedK.s = HomodynedK.k*HomodynedK.sigma;
        HomodynedK.sigma = sqrt(1/(HomodynedK.k^2*HomodynedK.mu+2));
        HomodynedK.s = sqrt(HomodynedK.mu)*HomodynedK.k*HomodynedK.sigma;
    case 'Dutt'
        A2 = mean(normalized_amp.^2);
        A4 = mean(normalized_amp.^4);
        A6 = mean(normalized_amp.^6);        
        F = fsolve(@HKequation,[0; 0.01; 0.2]);
        HomodynedK.s = sqrt(abs(F(1)));
        HomodynedK.mu = 1/F(2);
        HomodynedK.sigma = sqrt(F(3));
        HomodynedK.k = HomodynedK.s/HomodynedK.sigma;
end
HomodynedK.bin = bin;
HomodynedK.pdf = homok_func(bin,HomodynedK.mu,HomodynedK.s,HomodynedK.sigma^2);

end

%k = mu = err = -1111: when min(Datakurtosis) < min(TheoreticalKurtosis)


% 
% maxn=90;
% s=10;
% sigma = 1;
% mu = 1;
% x = 0:0.1:6;
% y2=zeros(maxn+1,max(size(x)));
% 
% term1=sqrt(2*mu/(pi*sigma)).*(x./s).^0.5./gamma(mu);
% % term2= (-1)^n * gamma(n+0.5)./(factorial(n)*gamma(0.5-n)).*(sigma/(s.*x.*mu).^n);
% term3=mu/2/sigma;
% term4=sqrt(2*mu/sigma).*abs(s-x);
% for n=0:maxn
%     y2(n+1,:)=term1.*...
%         ((-1)^n * gamma(n+0.5)./(factorial(n)*gamma(0.5-n)).*(sigma./(s.*x.*mu).^n)).*...
%         term3^((mu+n-1)/4)....
%         .*besselk(n+mu-1/2,term4);
% end
% y2(isnan(y2))=0;
% y2(y2==inf)=0;
% y2(y2==-inf)=0;
% y=sum(y2,1);


