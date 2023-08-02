

function [ value ] = simplif_hyper_term2( alpha,gama )


value=0;
oldvalue=value;
k=0;
vect=[(1:k),k]+alpha;
value=value+alpha/factorial(k)*1/prod(vect)*gama^k;

% figure,hold on
% plot(k,value,'*k')

while abs(value-oldvalue)>0.001
    
    oldvalue=value;
    k=k+1;
    vect=[(1:k),k]+alpha;
    value=value+alpha/factorial(k)*1/prod(vect)*gama^k;
%     plot(k,value,'*k')
    
end

% hold off

end

