% ����Լ��֮��ĳͷ�����
function y=sum_constraint(x)
% variable=[R_wingeff,C_avereff,xr,C_maxy];     % ����
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ������Լ���߽硪�����������������R_wingeff, C_avereff, xr��x_0(Ťת��ǰԵ���ֵ��ľ���)
% nvars = 4;         % Number of variables 
% R_wingeff��[2,4]*10^(-3);
% C_avereff��[0.5,2]*10^(-3);
% xr��[0,2]*10^(-3);
%%%%%%%%%%%%%%%%%%%%%%%%%
% x_0��[0,0.5];  % ��ǰԵ�������پ���
% C_maxylb =0.464385778290230;    % C_maxylb =0.4644;
% C_maxy25 =0.138924474377504;   % C_maxy25 =0.1389; 
% C_maxyub =-0.186536829535222; % C_maxyub =-0.1865
%%%%%%%%%%%%%%%%%%%%%%%%%
% variable=[R_wingeff,C_avereff,xr,C_maxy];     % ����
LB = [2, 0.5, 0, 0];      % Lower bound       % ��λ: mm & Ťת���ǰԵ�������پ���[0,0.5]
UB = [4, 2, 2, 0.5];      % Upper bound      % ��λ: mm & Ťת���ǰԵ�������پ���[0,0.5]
% X_AspR=[R_wingeff,C_avereff,xr];                   % AR=(R_wingeff+xr)/C_avereff; % AR��[2,8]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=length(x);
con_min=LB;
con_max=UB;
zeta=zeros(N,1);
% y=zeros(N,1);
for i=1:N
    if x(i) < con_min(i)
        zeta(i)=abs(con_min(i)-x(i))/(con_max(i)-con_min(i));
    elseif x(i) > con_max(i)  % x=[4.0000,0.5000,1.3364,0.5000];  penaltyfun2=663.6000=2000*sum(zeta);��y=sum(zeta)=0.3318;
        zeta(i)=abs(x(i)-con_max(i))/(con_max(i)-con_min(i));
    else % x(i)>=con_min(i) && x(i)<=con_max(i);  % ���������ʽ
        zeta(i)=0;  % disp('����δ�����߽�Լ��');  
    end
end
 y=sum(zeta);
end