%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function obj_function=morphology_para_opt_obj_function(x)
% morphology_para_opt_obj_function
% morphology parametric optimization objective function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% variable=[R_wingeff,C_avereff,xr,C_maxy];     % ����
% LB = [2, 0.5, 0, 0];      % Lower bound            % ��λ: mm & Ťת���ǰԵ�������پ���[0,0.5]
% UB = [4, 2, 2, 0.5];     % Upper bound            % ��λ: mm & Ťת���ǰԵ�������پ���[0,0.5]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ��һ���֡�����������
% ԭʼ��Ӭ��� R_wing=3.004;  C_aver=0.8854;   
% C_lead_ymax =0.4644;  C_trail_ymin =-0.8374;  C_max_LtoT = 1.3018; @C_maxy=0;
% Ratio_leadmax=C_lead_ymax/C_max_LtoT;
% %Ratio_leadmax=0.4644/1.3018=0.356737; %��Գ���ͳ�����ߵ�Ťת�ᶨ����ǰԵ/���ǰԵ������º�Ե��������=0.3567
% x=[3.004,0.8854,0.3289,0.356737]; % ��ʼֵ. % δ�Ż��Ĺ�Ӭ�����ò����������ͳ�����ߵ�Ťת�ᶨ����ǰԵ/���ǰԵ������º�Ե��������=0.3567
% C_lead_ymax =0.3255;  C_trail_ymin =-0.9764;  C_max_LtoT =1.3018;
% Ratio_leadmax=C_lead_ymax/C_max_LtoT; %Ratio_leadmax=0.3255/1.3018=0.25;%���Ťת��λ�����ǰԵ������º�Ե��������0.25��ʱ
% x=[3.004,0.8854,0.3289,0.25]; % ��ʼֵ. % δ�Ż��Ĺ�Ӭ�����ò������Ťת��λ�����ǰԵ������º�Ե��������0.25��ʱ
% x=[3.004,0.8854,0.3289,0.5]; % ��ʼֵ. % δ�Ż��Ĺ�Ӭ�����ò������Ťת��λ�����ǰԵ������º�Ե��������0.5��ʱ
% x=[100,33.7,0.3289,0.356737];  %��֤���Ŵ��ĳ���ת�����������ǳ���
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ������λ: mm & Ťת���ǰԵ�����پ���[0,0.6509]����ע��������ܲ��á���XXX
% x=[3.004,0.8854,0.3289,0.1389];     % δ�Ż��Ĺ�Ӭ�����ò����
% C_maxy25_nd=(0.464385778290230-0.138924474377504)/1.3018=0.3255/1.3018;  %C_maxy25_nd=0.25;
% x=[3.004,0.8854,0.3289,0.3255]; %��򱻷Ŵ�֮ǰ,Ťת���ڳ�����ǰԵ�㵽��С��Ե��ľ���1/4λ��@0.25*C_max_LtoT, ����һ��ָ��Ľ���
% Ratio_leadmax=C_lead_ymax/C_max_LtoT; %Ratio_leadmax=0.4644/1.3018=0.356737; %��Գ���ͳ�����ߵ�Ťת�ᶨ����ǰԵ
% x=[3.004,0.8854,0.3289,0.4644];     %��򱻷Ŵ�֮ǰ,Ťת���ڳ�����������@0.356737*C_max_LtoT, ����һ��ָ��Ľ���
% x=[3.004,0.8854,0.3289,0];          % ��򱻷Ŵ�֮ǰ,Ťת���ڳ��ǰԵ���ֵ�㴦@0*C_max_LtoT
% x=[3.004,0.8854,0.3289,0.6509]; % ��򱻷Ŵ�֮ǰ,Ťת���ڳ�����ǰԵ�㵽��С��Ե��ľ���һ��λ��@0.5*C_max_LtoT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%@ options = gaoptimset('InitialPopulation',x_start,'TimeLimit',14400);
% x =[3.6470,1.9176,0.6487,0.4258]; 
%@options=gaoptimset('InitialPopulation',x_start,'TimeLimit',14400,'HybridFcn',{@fminsearch,fminsearchOptions});
% x =[4.0000,1.3086,0.9363,0.4644]; 
% @options=gaoptimset('InitialPopulation',x_start,'TimeLimit',14400,'HybridFcn',{@fminsearch,fminsearchOptions});
% x =[4.0000,1.2895,0.9238,0.4644]; %ConstraintFunction=@aspectratio_constraint; 
% @options=gaoptimset('InitialPopulation',x_start,'TimeLimit',14400,'HybridFcn',{@fmincon,fminconOptions});
% x =[2.0000,1.1161,1.0711,-0.1865];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x=[4.0000,0.5000,1.3364,0.5000];  % ������
% x=[4.0000,1.4283,0.9992,0.5000];  % ������
x=[2.2899,1.0561,0.3788,0.0001];  % 20150204�����������Ŷ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P_total=Aero_M_fruitfly2_exp(x);     % �����ʡ���ƽ�����ʺ�Ťת���ʡ���Ťת�Ṧ�ʺ��Ĵ��Ṧ��
% size(P_total)                                      % (1000*2)
P_totalx=P_total(:,1);    % (1000*1)
P_totalz=P_total(:,2);    % (1000*1)
% xlswrite('D:\KXJ\PassiveRot_dynamic_Science_fruitfly\optimal_wing_para_motion\P_total.xlsx',P_total,'sheet1','A1:B2000');
%%%%%%%%%%%%%%%%%%%%%%%%%
N=length(P_totalx);
for i=1:N
    if P_totalx(i)<=0
        P_totalx(i)=0;
    end
end
P_totalx_posi=P_totalx;
%%%%%%%%%%%%%%%%%%%%%%%%%
for j=1:N  % N=length(P_totalx);
    if P_totalz(j)<=0
        P_totalz(j)=0;
    end
end
P_totalz_posi=P_totalz;
%%%%%%%%%%%%%%%%%%%%%%%%%
% P_total_posi=[P_totalx_posi,P_totalz_posi]; 
% xlswrite('D:\KXJ\PassiveRot_dynamic_Science_fruitfly\optimal_wing_para_motion\P_total_positive.xlsx',P_total_posi,'sheet1','A1:B1000');
%%%%%%%%%%%%%%%%%%%%%%%%%
f=188.7; T=1/f;  %����Ƶ�� (Hz)������  % w =1185.6; 
t=linspace(0.0052824335,0.0052824335+3*T,1000);  % t_steady1 
P_totalx_aver=trapz(t,P_totalx_posi)/(3*T);  % ʱ�������ʡ���ƽ�����ʺ�Ťת���ʡ���Ťת�Ṧ�ʺ��Ĵ��Ṧ��
P_totalz_aver=trapz(t,P_totalz_posi)/(3*T);  % ʱ�������ʡ���ƽ�����ʺ�Ťת���ʡ���Ťת�Ṧ�ʺ��Ĵ��Ṧ��
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m_insect=1.8;         % M_fly=1.8mg;  % 2014-Science-MH Dickinson     % P_asterisk =31.8058; % ���ٵ�λ=uW/mg;
% m_insect=0.72;    % 0.72 mg  % 2007-JFM-Wang ZJ                           % P_asterisk =80.6538; % ���ٵ�λ=uW/mg;
% m_insect=0.96;    % 0.96 mg  % 2005-JEB-MH Dickinson                    % P_asterisk =60.4903; % ���ٵ�λ=uW/mg;
% m_insect=1.05;    % 1.05mg; % 1997-JEB-Fritz-Olaf lehmann              % P_asterisk =55.3054; % ���ٵ�λ=uW/mg;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P_asterisk=(P_totalx_aver+P_totalz_aver)/m_insect    % ���ٵ�λ=uW/mg;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% �ڶ����֡���penaltyfun1�������ر�Լ�������ͷ�����
% % x=[3.004,0.8854,0.3289,0.4644];  % δ�Ż��Ĺ�Ӭ�����ò������������ͳ�����ߵ�Ťת�ᶨ����ǰԵ
% x=[3.004,0.8854,0.3289,0.356737]; % ��ʼֵ. % δ�Ż��Ĺ�Ӭ�����ò����������ͳ�����ߵ�Ťת�ᶨ����ǰԵ/���ǰԵ������º�Ե��������=0.3567
% % x=[3.004,0.8854,0.3289,0.1389];  % δ�Ż��Ĺ�Ӭ�����ò����
% x=[3.004,0.8854,0.3289,0.25];     % ��ʼֵ. % δ�Ż��Ĺ�Ӭ�����ò������Ťת��λ�����ǰԵ������º�Ե��������0.25��ʱ
% @ options = gaoptimset('InitialPopulation',x_start,'TimeLimit',14400);
% x =[3.6470,1.9176,0.6487,0.4258]; 
%@options=gaoptimset('InitialPopulation',x_start,'TimeLimit',14400,'HybridFcn',{@fminsearch,fminsearchOptions});
% x =[4.0000,1.3086,0.9363,0.4644]; 
% @options=gaoptimset('InitialPopulation',x_start,'TimeLimit',14400,'HybridFcn',{@fminsearch,fminsearchOptions});
% x =[4.0000,1.2895,0.9238,0.4644]; %ConstraintFunction=@aspectratio_constraint;  
% @options=gaoptimset('InitialPopulation',x_start,'TimeLimit',14400,'HybridFcn',{@fmincon,fminconOptions});
% x =[2.0000,1.1161,1.0711,-0.1865];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x=[4.0000,0.5000,1.3364,0.5000];  % ������
% x=[4.0000,1.4283,0.9992,0.5000];  % L =-2.6131;  % ������
% x =[4,2,0.8188,0.5];                         % L = 1.0009;
x=[2.2899,1.0561,0.3788,0.0001];  % 20150204�����������Ŷ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% F_verticalaver=trapz(t,F_vertical)/(3*T);
F_vertical_aver=Aero_F3_fruitfly_exp(x); %���ú���:x�Ǳ���variable=[R_wingeff,C_avereff,xr,C_maxy]; % F_verticalaver=12.3074uN; 
g=9.81;   % N*kg^-1
m_insect=1.8; 
M_insectweight=m_insect*g;                  % M_weight=10.3��1.27uN@1997-JEB-Fritz-Olaf lehmann
% L��ı�,��ΪŤת���λ��Ϊ�ı�ת����������ϵ��I6z�Ĵ�С;
% ���⻹��������������ϵ��(Z_rnd,M_xrdcoeff)������������ϵ��(I_xzam,I_xxam)�ı仯;
L=2*F_vertical_aver/M_insectweight      % ������: L =1.3940; 
delta=L-1  % delta =-1.0332e-004;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% �����ͷ���������penaltyfun1
r=2000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (L-1)<2.9*10^(-15)&&(L-1)>0
    penaltyfun1=0;
% elseif L<0 
%     penaltyfun1=r*abs((L-1));
elseif (L-1)>=2.9*10^(-15)||(L-1)<=0
    penaltyfun1=r*abs((L-1));
else
    disp('��ֱ�������д�');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % �÷��������е�Ŷ
% if ((L-1)>(2.9*10^(-15)) || (L-1)<0)
%     penaltyfun1=r*(L-1);    % penaltyfun1 =-1.6561e+004; 
% else
%    penaltyfun1=0;
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% penaltyfun1=r*heaviside(1-L);        % heaviside(x) has the value 0 for x < 0;         1 for x > 0;       and 0.5 for x = 0.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% �������֡���penaltyfun2�����Ż��������������½�Լ�������ͷ�����
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x=[4.0000,0.5000,1.3364,0.5000];  % ������
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fcn_con=sum_constraint(x);
% �����ͷ���������penaltyfun2
s=2000;
penaltyfun2=s*fcn_con;  % penaltyfun2 =663.6000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ���Ĳ��֡���penaltyfun3����չ�ұ�Լ������Rossby��Լ�������ͷ�����
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x=[4.0000,0.5000,1.3364,0.5000];  % ������
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% penaltyfun3=aspectratio_constraint(x);
% obj_function=P_asterisk+penaltyfun1+penaltyfun2+penaltyfun3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
obj_function=P_asterisk+penaltyfun1+penaltyfun2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%