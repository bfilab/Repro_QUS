              %% HK distribution
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                [val bin] = hist(ROI_env(:),100);
                nor_val = val/sum(val*(bin(2)-bin(1)));

%                 HomodynedK = make_homodynedK_PDF(ROI_env,bin);
% 
%                 HomodynedK.alpha = HomodynedK.mu;
%                 HomodynedK.sigma = sqrt(1/( HomodynedK.k^2 * HomodynedK.alpha +2));  
%                 HomodynedK.epsilon = sqrt(HomodynedK.alpha) *  HomodynedK.k*HomodynedK.sigma;
% 
%                 HomodynedK.pdf = homok_func(bin,HomodynedK.alpha,HomodynedK.sigma^2,HomodynedK.epsilon);
%                           
%                 paramsHK = HomodynedK;
%                 muHK(QUS_Vox_Z,QUS_Vox_X) = log10(paramsHK.mu);
%                 kHK(QUS_Vox_Z,QUS_Vox_X) = paramsHK.k;
%               
                try
                clear kkk mumu err;
                    [kkk Homodyned_K.mumu err] = estimator_RSK(double(ROI_env(:)),0);          %�@��Ԗڂ̕ϐ���1�ɂ���ƁC�������̃O���t���o�́C��ʂ�0�ɂ���
                    Homodyned_K.omegaHK = 1/(kkk^2+2);    % s^2+2*��^2=1�Ɓ@k=s/�� �ɂ��
                    Homodyned_K.s = kkk*sqrt(Homodyned_K.omegaHK);
                    Homodyned_KPDF = homok_func(bin,Homodyned_K.mumu,Homodyned_K.s,Homodyned_K.omegaHK);
                catch               
                    Homodyned_K.RMSE = RMSE(Homodyned_KPDF,PDF);
                    Homodyned_K.KL = KL(Homodyned_KPDF,PDF);
                    Homodyned_K.KS = KS(Homodyned_KPDF,PDF);   
                end
                                
                PDF = nor_val;
                figure; % ���̂܂܂��� ROI�̐������o�Ă��܂��p�j�b�N�ɂȂ�܂��D
                bar(bin ,PDF,'displayname','echo probabirity density')
%                 bar(HomodynedK.bin,nor_val,'displayname','echo probabirity density')
                hold on; plot(bin,Homodyned_KPDF,'displayname',sprintf('Fitted Nakagami PDF \mu = %f, \alpha = %f',log10(paramsN(1)),log10(paramsN(2))))                
                                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%