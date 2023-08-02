function [KL]=calc_KL(p_data,modelpdf,div)
% p_data:histcで計算されるデータの「確率分布」．N-1×1
% modelpdf:モデルの 「確率密度分布」.N×1
% div : histc および modelpdf計算に用いた刻み

% モデルの確率密度分布を台形近似で確率密度に変形する
p_model=(modelpdf(1:end-1)+modelpdf(2:end))/2;% 区間 [x(n), x(n+1)] の y 中間値　点数はN-1になる
p_model=p_model*div;% 台形積分による区間[x(n), x(n+1)] の確率
%disp(sum(p)) % 和が 1 になるかの確認用

zz=p_data./p_model;

KL=p_data.*log(zz);
KL=KL.*(KL~=inf);
KL=KL.*(KL~=-inf);
KL(isnan(KL)) = [];
KL = sum(KL);
%disp(['  KL=' num2str(KL) '[bit]']) % 結果の確認用
end