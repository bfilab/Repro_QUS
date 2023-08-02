

function [ value ] = simplif_hyper_term3( alpha,gama )


value=0;
oldvalue=value;
k=0;
vect=[(2:k+1),k+1]+alpha;
value=value+(alpha+1)/factorial(k)*1/prod(vect)*gama^k;

% figure,hold on
% plot(k,value,'*k')

while abs(value-oldvalue)>0.001
    
    oldvalue=value;
    k=k+1;
    vect=[(2:k+1),k+1]+alpha;
    value=value+(alpha+1)/factorial(k)*1/prod(vect)*gama^k;
%     plot(k,value,'*k')
    
end

% hold off

end

