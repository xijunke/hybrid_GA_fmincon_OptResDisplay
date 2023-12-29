% optimal_wing_para_motion
%% ��һ���֡���������Լ���߽�
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (1) ��򳤶Ⱥ�չ�ұȡ���R_wingeff,C_avereff��AR
% R_wingeff=3.004;      % ��Ч��򳤶�(mm)  
% C_avereff=0.8854;     % mm
% (2) ���Կ�����֪չ�ұ�, ������������ŵĳ���ò
% AR=R_wingeff/C_avereff;     %aspect ratio: AR=R^2/A_w;  �������Ϊ: AR=3.40158;  % Science������: AR=3.1519;
% C_avereff=R_wingeff/AR;     %mean chord length: C_aver=A_w/R=R^2/AR/R=R/AR; % C_aver=0.884mm;
% A_w=R_wingeff^2/AR;        %Area of wing: �������Ϊ: A_w=2.66mm^2;   %RJ Wood��Ƶĳ��: A_w=2.84345 mm^2 
% (3) �����ƫ����롪��xr
% xr=0.3289;                     % x-root offset  \mm
% xr_nd=xr/R_wingeff;      % x-root offset  ������չ��ƫ�þ���
% (4) Ťת���λ�á���C_maxy
% C_max_LtoT=C_lead_ymax-C_trail_ymin;       % C_max_LtoT =1.3018;
% x_0lb=0*C_max_LtoT;                                    % x_0lb=0;
% x_025=0.25*C_max_LtoT;                              % x_025=0.3255;
% x_0ub=0.5*C_max_LtoT;                                % x_0ub=0.6509;
% C_lead_ymax=max(f_x_lead2); % ���: C_lead_ymax=0.4644; k_leadmax=644;%ԭʼ��Ӭ���@���ǰԵ�㵽���������ߵľ���;@C_maxy=0;
% % C_maxylb =0.464385778290230;
% C_maxy25 =0.138924474377504;  % ��Գ���wing_model_88_yaxis��: ��122��; C_maxy =0.1389; 
% C_maxyub =-0.186536829535222; 
% C_maxy=C_maxylb;
% C_maxy=C_maxy25;
% C_maxy=C_maxyub;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% variable=[R_wingeff,C_avereff,xr,C_maxy];        % ����
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (5) ������Լ���߽硪�����������������R_wingeff, C_avereff, xr��x_0(Ťת��ǰԵ���ֵ��ľ���)
% nvars = 4;         % Number of variables 
% R_wingeff��[2,4]*10^(-3);
% C_avereff��[0.5,2]*10^(-3);
% xr��[0,2]*10^(-3);
%%%%%%%%%%%%%%%%%%%%%%%%%
% x_0��[0,0.5];  % Ťת���ǰԵ�������پ��롪��% ԭʼ��Ӭ��� R_wing=3.004;  C_aver=0.8854; 
% C_maxylb =0.464385778290230;                   % C_maxylb=0.4644;
% C_maxy25 =0.138924474377504;                  % C_maxy25=0.1389; 
% C_maxyub =-0.186536829535222;                % C_maxyub=-0.1865
%%%%%%%%%%%%%%%%%%%%%%%%%
% variable=[R_wingeff,C_avereff,xr,C_maxy];     % ����
% LB = [2, 0.5, 0, 0];      % Lower bound            % ��λ: mm & Ťת���ǰԵ�������پ���[0,0.5]
% UB = [4, 2, 2, 0.5];     % Upper bound            % ��λ: mm & Ťת���ǰԵ�������پ���[0,0.5]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% �ڶ����֡���Ŀ�꺯��
tic
ObjectiveFunction=@morphology_para_opt_obj_function;  % fitness and constraint functions��������Ŀ�꺯��
nvars=4;   % Number of variables 
LB = [2, 0.5, 0, 0];      % Lower bound       % ��λ: mm & Ťת���ǰԵ�������پ���[0,0.5]
UB = [4, 2, 2, 0.5];      % Upper bound      % ��λ: mm & Ťת���ǰԵ�������پ���[0,0.5]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% �������֡���չ�ұȶԱ���������Լ��
% X_AspR=[R_wingeff,C_avereff,xr];               % AR=(R_wingeff+xr)/C_avereff; % AR��[1,5]
ConstraintFunction=@aspectratio_constraint;  % nonlinear inequality constraints�������÷����Բ���ʽԼ������XXX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ���Ĳ��֡�������Ŵ��㷨(Hybrid Function)-GA+fminsearch
% ԭʼ��Ӭ��� R_wing=3.004;  C_aver=0.8854;   
% C_lead_ymax =0.4644;  C_trail_ymin =-0.8374;  C_max_LtoT = 1.3018; @C_maxy=0;
% Ratio_leadmax=C_lead_ymax/C_max_LtoT; %Ratio_leadmax=0.4644/1.3018=0.356737; %��Գ���ͳ�����ߵ�Ťת�ᶨ����ǰԵ
% x_start=[3.004,0.8854,0.3289,0.356737]; % ��ʼֵ. % δ�Ż��Ĺ�Ӭ�����ò������������ͳ�����ߵ�Ťת�ᶨ����ǰԵ
% C_lead_ymax =0.3255;  C_trail_ymin =-0.9764;  C_max_LtoT =1.3018;
% Ratio_leadmax=C_lead_ymax/C_max_LtoT; %Ratio_leadmax=0.3255/1.3018=0.25;%���Ťת��λ�����ǰԵ������º�Ե��������0.25��ʱ
x_start=[3.004,0.8854,0.3289,0.25];     % ��ʼֵ. % δ�Ż��Ĺ�Ӭ�����ò��������Ťת��λ�����ǰԵ������º�Ե��������0.25��ʱ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (1) �����Ŵ��㷨��С���Ż�
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% options = gaoptimset('InitialPopulation',x_start);  % (1)�������С������Ĵμ���(1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% options = gaoptimset('PlotFcns',{@gaplotbestf,@gaplotmaxconstr},'Display','iter');  % (2)��������
% options = gaoptimset(options,'InitialPopulation',x_start,'TimeLimit',14400);  % (2)�������С������Ĵμ���(2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (2) ����Ŵ��㷨(Hybrid Function)-GA+fminsearch
% % �Զ���fminsearch��ת��fmincon�����������С���������á������Ĵμ���(3)
% options = gaoptimset('InitialPopulation',x_start,'TimeLimit',10800,'PlotFcns',{@gaplotbestf,@gaplotmaxconstr});
% % fminsearchOptions = optimset('algorithm','Nelder-Mead simplex direct search');
% % options=gaoptimset(options,'HybridFcn',{@fminsearch,x_start, fminsearchOptions}); 
% fminsearchOptions=optimset('Display','iter'); % fminsearch������Լ�������������С��
% options=gaoptimset(options,'HybridFcn',{@fminsearch,fminsearchOptions}); % fminsearch������Լ�������������С��
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % ����fminsearch�ֲ���������������
% options=gaoptimset('PopulationSize',200,'InitialPopulation',[],'TimeLimit',14400,'PlotFcns',{@gaplotbestf,@gaplotmaxconstr});
% fminsearchoptions = optimset('Display','iter','TolFun',1e-10);% fminsearch������Լ�������������С��
% options =gaoptimset(options,'Hybridfcn',{@fminsearch,fminsearchoptions});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % ����fminunc�ֲ���������������
% options = gaoptimset('InitialPopulation',x_start,'TimeLimit',14400,'PlotFcns', {@gaplotbestf,@gaplotstopping});
% fminuncOptions=optimset('Display','iter', 'LargeScale','off'); % fminunc������Լ�������������С��
% options=gaoptimset(options,'HybridFcn',{@fminunc, fminuncOptions});  % fminunc������Լ�������������С��
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % ����fmincon�ֲ������������С���������á������Ĵμ���(4)
options = gaoptimset('InitialPopulation',x_start,'TimeLimit',10800,'PlotFcns',{@gaplotbestf,@gaplotmaxconstr});
fminconOptions=optimset('Display','iter','TolCon',1e-10); % fmincon������Լ�������Զ������С�������㷨��ѡ��
% fminconOptions=optimset('Display','iter','Algorithm','active-set'); % fmincon������Լ�������Զ������С�������㷨ѡ��
options=gaoptimset(options,'HybridFcn',{@fmincon,fminconOptions}); % fmincon������Լ�������Զ������С��
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Next we run the GA solver.
[x,fval]=ga(ObjectiveFunction,nvars,[],[],[],[],LB,UB, ConstraintFunction,options)
% [x,fval]=ga(ObjectiveFunction,nvars,[],[],[],[],LB,UB,[],options)
% [x,fval]=ga(ObjectiveFunction,nvars,[],[],[],[],[],[],[],options)
toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%