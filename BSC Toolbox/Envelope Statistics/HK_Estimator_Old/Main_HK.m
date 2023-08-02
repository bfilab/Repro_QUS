
clear all
load('ROI_env')




%%
[val bin] = hist(ROI_env(:),100);
nor_val = val/sum(val*(bin(2)-bin(1)));

%%
clear hk mumu err;
[hk Homodyned_K.mumu err] = estimator_RSK(double(ROI_env(:)),0);          %　二番目の変数は1にすると，等高線のグラフを出力，一般に0にする

% Homodyned_K.omegaHK = 1/(hk^2+2);    % s^2+2*σ^2=1と　k=s/σ による
Homodyned_K.omegaHK = mean(ROI_env(:).^2)/(hk^2+2);
Homodyned_K.s = hk*sqrt(Homodyned_K.omegaHK);
Homodyned_KPDF = homok_func(bin,Homodyned_K.mumu,Homodyned_K.s,Homodyned_K.omegaHK);

%%
PDF = nor_val;
figure; % このままだと ROIの数だけ出てしまいパニックになります．
bar(bin ,PDF,'displayname','echo probabirity density')
hold on; 
% plot(bin,Homodyned_KPDF,'displayname',sprintf('Fitted HK PDF \mu = %f, \alpha = %f',log10(paramsN(1)),log10(paramsN(2))))                
plot(bin,Homodyned_KPDF,'displayname',sprintf('Fitted HK PDF \mu = %f, s = %f  omega = %f',Homodyned_K.mumu,Homodyned_K.s,Homodyned_K.omegaHK))              

%%
% [kk mumu err] = estimator_RSK_v4_LN(ROI_current_Cyl(:), 0);
% copy form Kazuki's Liver analysis code
% clear kkk mumu err;
% [kkk Homodyned_K.mumu(i,j,k) err] = estimator_RSK(double(AmpDataCut(:)),0);          %　二番目の変数は1にすると，等高線のグラフを出力，一般に0にする
% Homodyned_K.omegaHK(i,j,k) = 1/(kkk^2+2);    % s^2+2*σ^2=1と　k=s/σ による
% Homodyned_K.s(i,j,k) = kkk*sqrt(Homodyned_K.omegaHK(i,j,k));
% 
% Homodyned_KPDF = homok_func(PDFbin,Homodyned_K.mumu(i,j,k),Homodyned_K.s(i,j,k),Homodyned_K.omegaHK(i,j,k));
% Homodyned_K.RMSE(i,j,k) = RMSE(Homodyned_KPDF,PDF);
% Homodyned_K.KL(i,j,k) = KL(Homodyned_KPDF,PDF);
% Homodyned_K.KS(i,j,k) = KS(Homodyned_KPDF,PDF);