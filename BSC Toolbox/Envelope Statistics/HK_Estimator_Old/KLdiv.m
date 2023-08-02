function [KL] = KLdiv(p,p_model)

if size(p,2)~=size(p_model,2)
    p=p';
end

KL=p.*log(p./p_model);
KL(isfinite(KL)~=1)=0;
KL = sum(KL);

end