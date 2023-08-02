function y=homok_func(x,mu,s,sigma)
% 2011.11.27

if size(x,1)~=1
    x=x(:);
    if size(x,1)~=1
        x=x';
    end
end

maxn=90;
y2=zeros(maxn+1,max(size(x)));

term1=2*mu*x/sigma/gamma(mu);
term2=mu*s*x/2/sigma;
term3=sqrt(mu*(s^2+x.^2)/2/sigma);
term4=sqrt(2*mu*(s^2+x.^2)/sigma);
for n=0:maxn
    y2(n+1,:)=term1/factorial(n)^2.*term2.^(2*n).*...
        term3.^(mu-1-2*n).*besselk(2*n+1-mu,term4);
end

y2(isnan(y2))=0;
y2(y2==inf)=0;
y2(y2==-inf)=0;
y=sum(y2,1);

end