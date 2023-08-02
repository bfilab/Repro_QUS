function y=myfun_KL(X1,X2)


temp = X1./X2;
% 
% temp2=temp(find(~isnan(temp)));  % method 1
% X1=X1(find(~isnan(temp)));
% temp3=temp2(find(temp2~=0));
% X1=X1(find(temp2~=0));
% y=sum(X1.*log(temp3));

k=1;
for i=1:length(temp)
    if temp(i)~=0 && ~isnan(temp(i)) && ~isinf(temp(i))
        tempX(k)=temp(i);
        X(k)=X1(i);
        k=k+1;
    end
end
y=sum(X.*log(tempX));
