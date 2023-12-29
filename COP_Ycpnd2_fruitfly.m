function Y_rcpnd=COP_Ycpnd2_fruitfly(alpha,xr0,C_maxyaxis)   % ���������C_maxyaxis��xr0, ���ú���R��C_aver;
% �޸�ʱ�䡪��2014��12��21��,23:36
% �޸�ʱ�䡪��2015��01��20��,11:16
%ת����ƫ�������е��ƫ��������,��Ťת������ƫ��C_maxy֮��
%% (1) �����ת���������ء�����⾻ѹ�ĵ�������λ��Y_cpnd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% R=16.0148-0.88=15.1348;
% �� r_nd��(0.88/15.1348,13.4427/15.1348)
% clear all; clc;  
R_wingeff=3.004;          %��Ч��򳤶�(mm) 
C_avereff=0.8854;         % mm
% xr=0.3289;                 % x-root offset  \mm
xr=xr0;                           % ����������������������������������������������
xr_nd=xr/R_wingeff;      % x-root offset  ������չ��ƫ�þ���
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ���-ƫ��-����ϵԭ��ľ���
R_proximal=xr;                                                    % xr=3.19;     %RJ Wood��Ƶĳ��\mm
R_distal=R_wingeff+xr;                                        % yr=0.73;    %RJ Wood��Ƶĳ��\mm
x=linspace(R_proximal,R_distal,200);
x_mod_Root=0.636;                                            % mm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% C_maxylb=0.464385778290230;
% C_maxy25=0.138924474377504;  % ��Գ���wing_model_88_yaxis��: ��122��; C_maxy =0.1389; 
% C_maxyub=-0.186536829535222; 
% C_maxy=C_maxylb;
% C_maxy=C_maxy25;
C_maxy=C_maxyaxis;      % ����������������������������������������������
% C_maxy=C_maxyub;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ԭʼ��Ӭ��� R_wing=3.004;  C_aver=0.8854;   
C_lead_ymax=0.4644;   % C_trail_ymin =-0.8374;  
C_max_LtoT= 1.3018;    % @C_maxy=0;
% Ratio_leadmax=C_lead_ymax/C_max_LtoT; %Ratio_leadmax=0.4644/1.3018=0.356737; %��Գ���ͳ�����ߵ�Ťת�ᶨ����ǰԵ
% x_start=[3.004,0.8854,0.3289,0.356737]; % ��ʼֵ. % δ�Ż��Ĺ�Ӭ�����ò������������ͳ�����ߵ�Ťת�ᶨ����ǰԵ
% C_lead_ymax =0.3255;  C_trail_ymin =-0.9764;  C_max_LtoT =1.3018;
% Ratio_leadmax=C_lead_ymax/C_max_LtoT; %Ratio_leadmax=0.3255/1.3018=0.25;%���Ťת��λ�����ǰԵ������º�Ե��������0.25��ʱ
% x_start=[3.004,0.8854,0.3289,0.25];     % ��ʼֵ. % δ�Ż��Ĺ�Ӭ�����ò��������Ťת��λ�����ǰԵ������º�Ե��������0.25��ʱ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C_maxy=C_lead_ymax-C_maxy*C_max_LtoT;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yr_lead=-0.08249*x.^6+0.9167*x.^5-4.04*x.^4+8.872*x.^3-10.06*x.^2+5.674*x-0.413-x_mod_Root-C_maxy;  
yr_trail=-0.0333*x.^6+0.504*x.^5-2.795*x.^4+7.258*x.^3-8.769*x.^2+3.739*x+0.1282-x_mod_Root-C_maxy;
% C_rx=yr_lead-yr_trail;      % ��ȷ�������ٻ�ʵ���ҳ��ֲ�
% ����ǰԵ��Ϻ�����⡪��ʵ��ƽ���ҳ�=���/R_wingeff
% wing_aera=trapz(x,C_rx);             %���: wing_aera =2.6597; % mm^2
% C_aver=wing_aera/R_wingeff    % ������ٻ�ƽ���ҳ�: C_avereff=C_aver =0.8854; % mm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (a) ������ǰԵ�ֲ�����
r_nd=(x-xr)/R_wingeff;
yr_leadnd0=yr_lead/C_avereff;
P_coeff_lead=polyfit(r_nd,yr_leadnd0,6);
% (b)  �����ٺ�Ե�ֲ�����
yr_trailnd0=yr_trail/C_avereff;
% P_coeff_trail=polyfit(r_nd,yr_trailnd0,6);
% (c)  �������ҳ��ֲ�����
Cr_nd=yr_leadnd0-yr_trailnd0;
P_coeff_Cr=polyfit(r_nd,Cr_nd,6);  % ����ʽϵ��  % Cr_nd2=polyval(Coeff,r_nd1);
syms r_nd   % �������ҳ��ֲ�Ϊ6�׶���ʽ����ת������������ָ��
yr_leadnd=vpa(poly2sym(P_coeff_lead,r_nd),6); 
% yr_trailnd=vpa(poly2sym(P_coeff_trail,r_nd),6);
C_nd=vpa(poly2sym(P_coeff_Cr,r_nd),6);  
% % ����(1)������ǰ��Ե�����������ٻ�yr_leadnd����yr_trailnd�������
% �����ǳ�ǰԵ�����������Ťת���λ�ò�ͬ��Ҫ���зֶκ���������?
% yr_leadnd =-68.4639*r_nd^6+208.297*r_nd^5-245.231*r_nd^4+137.467*r_nd^3-36.8594*r_nd^2+4.79235*r_nd-0.156071;
% yr_trailnd =-27.6379*r_nd^6+121.093*r_nd^5-185.804*r_nd^4+125.603*r_nd^3-33.1053*r_nd^2-0.145527*r_nd-0.156013;
% �������ҳ��ֲ�Ϊ6�׶���ʽ
% C_nd =-40.826*r_nd^6+87.204*r_nd^5-59.4267*r_nd^4+11.8645*r_nd^3-3.75408*r_nd^2+4.93788*r_nd-0.0000578215;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% wing_kenimatics=kenimatics_wing_and_AoA();        %���ú���kenimatics_wing_and_AoA
% % size(wing_kenimatics)     %  (400*16)
% alpha=wing_kenimatics(:,4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d_cprnd=0.82*abs(alpha)/pi+0.05;                    %��ת���������ص�����ѹ��λ��;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % alpha1=(0.25-0.05)*pi/0.82*180/pi                % alpha1=43.9024;
% Ratio_axistr =0.2834;   % ƽ������չ��ѹ�Ķ�Ӧ��Ƭ����Ťת�ᵽǰԵ�������پ���
% alpha1=(Ratio_axistr-0.05)*pi/0.82*180/pi;        % alpha1=51.2341;  
% Ratio_axisrot=0.2706;  % ת������չ��ѹ�Ķ�Ӧ��Ƭ����Ťת�ᵽǰԵ�������پ���
% alpha2=(Ratio_axisrot-0.05)*pi/0.82*180/pi;      % alpha2=48.4244;  
% alpha_crit=[+19.5295,-26.2696, 93.28874372]*pi/180;  % ���ι��ǵļ����ٽ�ֵ
% d_cprnd_crit=0.82*abs(alpha_crit)/pi+0.05;        % ��ת���������ص�����ѹ��λ��;
% % output: d_cprnd_crit =[0.1390    0.1697    0.4750];  % ����ѹ�ı仯������������
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% yr_cpnd=yr_nd+yr_leadnd-C_nd*d_cprnd;  
yr_cpnd=yr_leadnd-C_nd*d_cprnd;  
fx2=(r_nd+xr_nd)^2*C_nd;                  % ������������F_nd��ԭʼ��������
fx4=expand(fx2);
F_nd=double(int(fx4,r_nd,0,1));           % Result: F_nd=0.46392
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ����������������(nondimention_aerodynamic_component)����⡪��F_nd
% ע�⡪���öγ����мǲ����޸ģ�ǰ��ֻҪ��֤������ȷ���������ҳ��ֲ����ɡ�
%���µĹ�ʽӦʹ�ú���������ٵ��ҳ��ֲ���ʽC_nd
% R2nd2=double(int(r_nd^2*C_nd,r_nd,0,1)); %��������صĻ�ת�뾶��ƽ��
% R1nd1=double(int(r_nd*C_nd,r_nd,0,1));     %һ������صĻ�ת�뾶
% % S_nd=double(int(C_nd,r_nd,0,1));                %�����ٳ����
% F_nd2=R2nd2+2*xr_nd*R1nd1+xr_nd^2;    %ʹ����������Ҳ��ȷ; ���:F_nd2 =0.5024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Y_rcpnd=int(yr_cpnd*fx4,r_nd,0,1)/F_nd;
% disp(['��ѹ�ĵ�������λ��Y_cpnd(alpha)=' num2str(Y_rcpnd)  ' ���ٵ�λ������mm'])
% Y_rcpnd=abs(double(int(yr_cpnd*fx4,r_nd,0,1))/F_nd);
Y_rcpnd=double(int(yr_cpnd*fx4,r_nd,0,1))/F_nd;

