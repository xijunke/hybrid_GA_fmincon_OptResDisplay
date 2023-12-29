 function [c, ceq]=aspectratio_constraint(x)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % function penaltyfun3=aspectratio_constraint(x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2009-JEB-David Lentink: AR��[1,5], ���������AR��[2.5,3.5]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% variable=[R_wingeff,C_avereff,xr,C_maxy];     % ����
  R_wing=x(1);
  C_aver=x(2); 
  xr0=x(3);
%  %  C_maxyaxis=x(4);
%  %%%%%%%%%%%%%%%%%%%
%  % չ�ұ�Լ�� % AR=(R_wingeff+xr)/C_avereff; % AR��[1,5]
% %   c = [R_wing+xr0- 5*C_aver;
% %           -R_wing-xr0+C_aver];
% %   ceq=[];
% %%%%%%%%%%%%%%%%%%%%%%
% չ�ұ�Լ�� % AR=(R_wingeff+xr)/C_avereff; % AR��[2.5,3.5]
  c = [R_wing+xr0-3.5*C_aver;
          -R_wing-xr0+2.5*C_aver];
  ceq=[];
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % LB = [2, 0.5, 0, 0];      % Lower bound       % ��λ: mm & Ťת���ǰԵ�������پ���[0,0.5]
% UB = [4, 2, 2, 0.5];      % Upper bound      % ��λ: mm & Ťת���ǰԵ�������پ���[0,0.5]
% AR=(R_wingeff+xr)/C_avereff;
%%%%%%%%%%%%%%%%%%%%%%%%%%
% AR1=(4+2)/0.5=12;
% AR2=(4+2)/2=3;
% AR3=(4+0)/0.5=8;
% AR4=(4+0)/2=2;
%%%%%%%%%%%%%%%%%%%%%%%%%%
% AR11=(2+2)/0.5=8;
% AR22=(2+2)/2=2;
% AR33=(2+0)/0.5=4;
% AR44=(2+0)/2=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%
% չ�ұ�Լ�� % AR=(R_wingeff+xr)/C_avereff;  % AR��(1,12);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %  x=[3.004,0.8854,0.3289,0.25]; % AR=3.76429;
% if (R_wing+xr0-3.5*C_aver)>0 || (-R_wing-xr0+2.5*C_aver)>0
%     penaltyfun3=2000;  % �����ͷ���������penaltyfun3  ����% disp('չ�ұ�Լ��ʧЧ,��Rossby��������');
% else
%     penaltyfun3=0;
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
