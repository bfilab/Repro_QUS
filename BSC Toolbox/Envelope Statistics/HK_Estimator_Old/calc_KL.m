function [KL]=calc_KL(p_data,modelpdf,div)
% p_data:histc�Ōv�Z�����f�[�^�́u�m�����z�v�DN-1�~1
% modelpdf:���f���� �u�m�����x���z�v.N�~1
% div : histc ����� modelpdf�v�Z�ɗp��������

% ���f���̊m�����x���z���`�ߎ��Ŋm�����x�ɕό`����
p_model=(modelpdf(1:end-1)+modelpdf(2:end))/2;% ��� [x(n), x(n+1)] �� y ���Ԓl�@�_����N-1�ɂȂ�
p_model=p_model*div;% ��`�ϕ��ɂ����[x(n), x(n+1)] �̊m��
%disp(sum(p)) % �a�� 1 �ɂȂ邩�̊m�F�p

zz=p_data./p_model;

KL=p_data.*log(zz);
KL=KL.*(KL~=inf);
KL=KL.*(KL~=-inf);
KL(isnan(KL)) = [];
KL = sum(KL);
%disp(['  KL=' num2str(KL) '[bit]']) % ���ʂ̊m�F�p
end