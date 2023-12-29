function F_verticalaver=Aero_F3_fruitfly_exp(x) % x�Ǳ�������% variable=[R_wingeff,C_avereff,xr,C_maxy];
%% ������: Solution of the aerodynamic force
% �޸�ʱ�䡪��2014��12��20��,11:54����ת����ƫ�������е��ƫ��������,��Ťת������ƫ��C_maxy֮��
% clear all;clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ���ú�����ò��������������ϵ���ĺ���
% wing_para_output=zeros(1,16);
% wing_para_output=[F_ndTrans,Coeff_liftdragF_N,M_xaercoeff,I1y,Z_rnd, M_xrdcoeff,...
%     F_ndRot,F_yrotcoeff,M_xRotcoeff, I2y,...
%     I_xzam,I_xxam, I5z,I6z,I7y,M_zrdcoeff];
R_wing=x(1);
C_aver=x(2); 
xr0=x(3);
C_maxyaxis=x(4);
wing_para=wing_shape_fruitfly_sixteen_good(R_wing,C_aver,xr0,C_maxyaxis);   %���ú���wing_shape_fruitfly;  % size(wing_para)
% F_ndTrans=wing_para(1,1);             % F_ndTrans=0.46392;  ������, ���ٻ���λΪmm^4
Coeff_liftdragF_N=wing_para(1,2);      % Coeff_liftdragF_N=0.00682;  %��λ��mg*mm
% M_xaercoeff=wing_para(1,3);          % M_xaercoeff=0.006038;   %��λ��: mg.mm^2
% I1z=wing_para(1,4);                         % I1y=0.016158   % ��λ�� mg.mm^2
% Z_rnd=wing_para(1,5);                     % Z_rnd=0.16265;  ������, ���ٻ���λΪmm
% M_xrdcoeff=wing_para(1,6);           % M_xrdcoeff=0.0001839; % ��λ��mg.mm^2
% F_ndRot=wing_para(1,7);                % F_ndRot=0.74851;  ������, ���ٻ���λΪmm^4
F_yrotcoeff=wing_para(1,8);               % F_yrotcoeff =0.0032433;  % ��λ�� mg.mm
% M_xRotcoeff=wing_para(1,9);         % M_xRotcoeff=0.002871;   % ��λ�� mg.mm^2
% I2z=wing_para(1,10);                      % I2y=0.006943;        % ��λ�� mg.mm^2
% I_xzam=wing_para(1,11);                % I_xzam=0.001424  % ��λ�� mg.mm^2
% I_xxam=wing_para(1,12);                % I_xxam=0.000338  % ��λ�� mg.mm^2
I5y=wing_para(1,13);                          % I5z=0.0050926   % ��λ�� mg.mm
I6y=wing_para(1,14);                          % I6z=0.00077164  % ��λ�� mg.mm
% I7z=wing_para(1,15);                      % I7y=0.0109056;        % ��λ�� mg.mm^2
% M_zrdcoeff=wing_para(1,16);         % M_zrdcoeff=0.0011169; % ��λ�� mg.mm % ת�������������ز������Ƴ�ƽ���µ�������
% C_max_LtoT=wing_para(1,17);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I5y=wing_para(1,13);                       % I5z=0.0050945   % ��λ�� mg.mm        % ������I3yӦ�ø�ΪI5y
% I6y=wing_para(1,14);                       % I6z=0.0011         % ��λ�� mg.mm         % ������I4yӦ�ø�ΪI6y
% I7z=wing_para(1,15);                       % I7y=0.0109;        % ��λ�� mg.mm^2    % ������I5zӦ�ø�ΪI7z
% I1z=wing_para(1,4);                         % I1y=0.0162;        % ��λ�� mg.mm^2    % ������I7zӦ�ø�ΪI1z
% I2z=wing_para(1,10);                       % I2y=0.0069;        % ��λ�� mg.mm^2    % ������I6zӦ�ø�ΪI2z 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ���ó���˶�ѧ-���˶����ɺͼ��ι���(AOA)������
% wing_m_output=[t',phi',psi',alpha',dphi',dpsi',ddphi',ddpsi',C_L',C_D',C_N1',C_T'];
% wing_kenimatics=kenimatics_wing_and_AoA_fruitfly_sim();  %���ú���kenimatics_wing_and_AoA;  % size(wing_kenimatics)  % (1000,11)
wing_kenimatics=kenimatics_wing_and_AoA_fruitfly_exp();
t=wing_kenimatics(:,1);                % ��λ��ms
% phi=wing_kenimatics(:,2);        % �Ĵ�ǡ�����λ��rad
psi=wing_kenimatics(:,3);            % �Ĵ�ǡ�����λ��rad
alpha=wing_kenimatics(:,4);      % alpha=atan2(omega_z,-omega_y);  % ���ι��ǡ���������   %������������и�
dphi=wing_kenimatics(:,5);          % ��λ��rad/s
dpsi=wing_kenimatics(:,6);          % ��λ��rad/s
ddphi=wing_kenimatics(:,7);       % ��λ��rad/s^2
ddpsi=wing_kenimatics(:,8);     % ��λ��rad/s^2
C_L=wing_kenimatics(:,9);          
C_D=wing_kenimatics(:,10);     
C_N1=wing_kenimatics(:,11);   
C_T=wing_kenimatics(:,12);   
C_N=C_N1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% �������ϵ�µĽ����ʺͽǼ����ʡ����������������������Գ�2DOF�˶�
% �������ϵ�µĽ��ٶ�
omega_x=dpsi;                     % չ��
omega_y=dphi.*sin(psi);       % ����(��ʼ����)
omega_z=dphi.*cos(psi);      % ������
omega_h=dphi;       % �����Ľ��ٶ�% omega_h=-sqrt(omega_y^2+omega_z^2);  % ���ַ���˳ʱ��
% �������ϵ�µĽǼ��ٶȡ����������������ļ���
domega_x=ddpsi;
domega_y=ddphi.*sin(psi)+dphi.*dpsi.*cos(psi);
domega_z=ddphi.*cos(psi)-dphi.*dpsi.*sin(psi); 
%% �������Ǻ͵��������ٶȵļ���
v_y_nonr=omega_z;   % v_y=r*dphi*cos(psi)
v_z_nonr=-omega_y;    % v_z=-r*dphi*sin(psi)
% alpha2=atan2(v_y_nonr,v_z_nonr);   % ��ȷ����ע�������ĵ�alpha=atan2(omega_z,-omega_y)*180/pi; ��ͬ
% % ����alpha=atan2(cot(psi))=atan2(cot(pi/2-alpha))=atan2(tan(alpha)); %����atan2������������ֵ��������alpha>pi/2ʱ
V_nonr=sqrt(v_y_nonr.^2+v_z_nonr.^2); % ���������ٶ�V_nonr=omega_h=dphi;   % ��λ�� rad/s
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% xr=0.3289;  % x-root offset  \mm      % R_wingeff=wing_para(1,1);  % R_wingeff =3.0040;   % ��λ�� mm
% R_ref=(xr+R_wingeff)*10^(-3);  
% R_ref_non=3.3293*10^(-3);
R_ref_non=1;   % r2_nd=0.5801;   
v_y=R_ref_non*omega_z;
v_z=-R_ref_non*omega_y;
V_ref=sqrt(v_y.^2+v_z.^2);
f=188.7;  T=1/f;          % Hz
V_ref_aver=trapz(t,V_ref)/(3*T);    % ���ο��ٶ�: V_ref_aver=867.0901rad/s;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ��һ������������������ϵ�£�ƽ�������������������������������õ㣺ƬԪ�ҳ��е�
%% �������Ϊ��4��
% C_avereff;  R_wingeff;    F_nd;                �����������������������Գ���ò������;% F_nd=0.50236;
% omega_h;  C_L(alpha);  C_D(alpha);       �����������������������Գ�2DOF�˶�����������, ����
%%  Coeff_liftdragF_N   %��λ  mg.mm=*10^(-9)kg.m
% g=9.821*10^(-6);   % �����������ٶ�:g=9.821N/kg=9.821*10^(-6)=9.821e-006   N/mg  ����g�Ĺ��ʵ�λ��m*s^-2��N*kg^-1
g=9.821;   % �����������ٶ�:g=9.821N/kg=9.821*10^(-6)=9.821e-006   N/mg  ����g�Ĺ��ʵ�λ��m*s^-2��N*kg^-1
lift_inst=V_nonr.^2.*C_L*Coeff_liftdragF_N*10^(-3)/g;      % ��λ��rad^2*s^-2*kg*m / (m*s^-2)=mg
drag_inst=V_nonr.^2.*C_D*Coeff_liftdragF_N*10^(-3)/g;  % ��λ��rad^2*s^-2*kg*m / (m*s^-2)=mg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%��˲ʱ����������ʹ�����λ��ֺ���trapz������ֵ���֣���������, ���ƽ��������
lift_aver=trapz(t,lift_inst)/(3*T)                            % lift_aver =1.0913 
drag_aver=trapz(t,drag_inst)/(3*T)                      % drag_aver =0.9887
drag_instabs_aver=trapz(t,abs(drag_inst))/(3*T);        
lift2drag_aver=lift_aver/drag_instabs_aver         % ������������ֵ: lift2drag_aver =1.1037������������ϵ���ı�ֵ��û��ʵ�����壿
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(10)
F_LD=plot(t/T,lift_inst,'r-',t/T,drag_inst,'b:','LineWidth',2);   
xlabel('\itNormalized time')
ylabel('\itAerodynamic force (mg)')                    %��λ���㵽mg��Ŷ
title('Aerodynamic lift and drag force \itvs. \itt \rm for flapping wing')
% legend('\itinstantaneous lift (\itF_L)','\itinstantaneous drag (\it|F_D|)')
legend('\itinstantaneous lift (\itF_L)','\itinstantaneous drag (\itF_D)')
set(F_LD,'LineWidth',2)
grid on
axis([0.9,3.0,-inf,inf])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ������ϵ�£�ƽ������������(1)���������������õ㣺ƬԪ�ҳ��е�
% C_N=C_N1=cos(abs(alpha)).*C_L+sin(abs(alpha)).*C_D; %��������ϵ���ϳɡ���2010-JFM-RJ Wood
% Coeff_liftdragF_N=6.8201e-012������λ��:kg/m^3*mm*mm^3= 10^(-12) kg.m
% F_ytran=-sign(alpha).*V_nonr.^2.*C_N*Coeff_liftdragF_N/g;  %��λ��rad^2*s^-2*kg*m / (m*s^-2)=mg������������ʾ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(11)
% plot(t/T,F_ytran,'k-','LineWidth',2)
% xlabel('\itNormalized time')
% ylabel('\itF_{y,tran} (mg)')
% legend('F_{y,tran}')
% title('ƽ������������(����)��ʱ��ı仯����') 
% grid on
% axis([0.9,4.05,-inf,inf])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ��λ��rad^2*s^-2*kg*m / (m*s^-2)=mg
F_ytran=-sign(alpha).*V_nonr.^2.*C_N*Coeff_liftdragF_N*10^(-3);   % ��λ��rad^2*s^-2*mg*mm=10^(-9)N=10^(-9)*10^6uN
F_ytran_aver=trapz(t,abs(F_ytran))/(3*T)   % F_ytran_aver =14.5887uN;
F_ztran=-sign(alpha).*V_nonr.^2.*C_T*Coeff_liftdragF_N*10^(-3);  
% F_ytran=V_nonr.^2.*C_N*Coeff_liftdragF_N*10^(-3); 
% F_ztran=V_nonr.^2.*C_T*Coeff_liftdragF_N*10^(-3);  
figure(12)
plot(t/T,F_ytran,'k-',t/T,F_ztran,'r-','LineWidth',2)
xlabel('\itNormalized time')
ylabel('\itF_{y,tran} & F_{z,tran} (uN)')
legend('F_{y,tran}','F_{z,tran}')
title('ƽ������������(���������)��ʱ��ı仯����') 
grid on
axis([0.9,4.05,-inf,inf])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% �ڶ������������� ������ϵ�£� ��ת�����������������������õ㣺ƬԪ�ҳ��е�
%% �������Ϊ��3��
% C_avereff;   R_wingeff;   F_ndRot;  �����������������������Գ���ò������
% omega_h;  omega_x;                    �����������������������Գ�2DOF�˶�����������, �������
% F_zrot=(1/2)*Rou*R_wingeff^2*C_avereff^2*C_R*omega_x.*omega_h*F_ndRot*(10^(-12)*10^3);%��λ��mN����δ���Ƿ���
% F_y_rot=C_R*sign(alpha2).*abs(omega_x).*abs(omega_h)*F_y_rotcoeff*10^6;    % ��λ��uN
C_R=1.55;
% F_yrot=-C_R*sign(alpha).*omega_x.*V_nonr*F_yrotcoeff*10^(-3); % 
F_yrot=C_R*omega_x.*V_nonr*F_yrotcoeff*10^(-3);      % ��λ��rad^2*s^-2*kg*m=10^(-9)N=10^(-9)*10^6uN
F_yrot_aver=trapz(t,abs(F_yrot))/(3*T)                             % F_yrot_aver =2.8944uN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(13)                  % ����������ЩС  
hold on
plot(t/T,F_ytran,'k-',t/T,F_yrot,'r-','LineWidth',2)
xlabel('\itNormalized time')
ylabel('\itF_{y,tran} &  F_{y,rot} (uN)')
legend('\itF_{y,tran}','\itF_{y,rot}')
title('ƽ������������(����)����ת����������(����)��ʱ��ı仯����')   
grid on
axis([0.9,4.05,-inf,inf])
set(gca,'XTick',(0.9:0.1:4.05))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ����������������������ϵ�£��������������������������õ㣺ƬԪ�ҳ��е�
%%% �Ĵ�ͱ���Ťת��Эͬ���ã������������������ο�SK Agrawal�������Ƶ����̣��Ա��������������
%% �������Ϊ��10��
% I5y;   I6y;     % ��λ�� mg.mm=10^(-9) kg.m �����������������������Գ���ò������
% domega_z; omega_x; omega_y;                      �����������������������Գ�2DOF�˶���Ǽ�����, ������
% dphi; dalpha; ddphi; ddalpha; alpha;               �����������������������Գ�2DOF�˶�������, �Ǽ�����, ����
%%  (1) 2010_JFM_Aeromechanics-��������������ʽ����% I5y��I6y��λ�� mg.mm=*10^(-9) kg.m
% F_yadd=sign(phi).*(I5y*(domega_z-omega_x.*omega_y)+I6y*domega_x)*(10^(-9)*10^6);   %��λ��uN
% F_yadd1=(I5y*(domega_z-omega_x.*omega_y)+I6y*domega_x)*10^(-3); % ԭ���Ƴ��Ĺ�ʽ������Ϊ����ЩŶ
% ע�⹫ʽǰ(-)���������ݹ�Ӭ�������ϵ�޸��˽��ٶȺͽǼ��ٶ�
% F_yadd1=-(I5y*(domega_z+omega_x.*omega_y)-I6y*domega_x/4)*10^(-3); %��λ��rad^2*s^-2*kg*m=N=10^6uN % ��2001-JEBƥ��
% F_yadd1=-(I5y*(domega_z+omega_x.*omega_y)+I6y*domega_x)*10^(-3); % �޸���ԭ���Ƴ��Ĺ�ʽ��������,�Ҳ���/4  % ������һ����
F_yadd1=(-I5y*(domega_z+omega_x.*omega_y)+I6y*domega_x)*10^(-3);    % �޸���ԭ���Ƴ��Ĺ�ʽ��������,�Ҳ���/4  % �����ڶ�����
F_yadd1_aver=trapz(t,abs(F_yadd1))/(3*T)    % F_yadd1_aver=6.0929uN(old);  F_yadd1_aver=4.5472;����ƽ��ֵ
% (2)������2001-JEB-Sanjay SP-��������������ʽ
% F_yadd2=(I5y*(ddphi.*sin(alpha)+dphi.*dalpha.*cos(alpha))-(I6y/4)*ddalpha/4)*10^(-3);% ԭ���Ƴ��Ĺ�ʽ��������
% �����alpha1Ϊȫ�����ι��ǣ������и���dalpha1��ddalpha1
% F_yadd2=-(I5y*(ddphi.*sin(alpha1)+dphi.*dalpha1.*cos(alpha1))-(I6y/4)*ddalpha1/4)*10^(-3); 
% F_yadd2=-(I5y*(ddphi.*sin(alpha1)+dphi.*dalpha1.*cos(alpha1))+(I6y/4)*ddalpha1/4)*10^(-3); 
% F_yadd2=-(I5y*(ddphi.*sin(alpha1)+dphi.*dpsi.*cos(alpha1))-(I6y/4)*ddpsi/4)*10^(-3); %�����dpsi��ddpsi����仯��Ҫȷ��
% F_yadd2=-(I5y*(ddphi.*sin(alpha1)+dphi.*dpsi.*cos(alpha1))+(I6y/4)*ddpsi/4)*10^(-3); %�����dpsi��ddpsi����仯��Ҫȷ��
% �����alpha2�����и�, dalpha��ddalpha��Ҫ��
% F_yadd2=-(I5y*(ddphi.*sin(alpha)+dphi.*dalpha.*cos(alpha))-(I6y/4)*ddalpha/4)*10^(-3); 
i=length(alpha);   % size(alpha)=(2000*1)
dalpha=[diff(alpha);alpha(i,1)];    
j=length(dalpha);
ddalpha=[diff(dalpha);dalpha(j,1)];  % size(ddalpha)
F_yadd2=-(I5y*(ddphi.*sin(abs(alpha))+dphi.*dalpha.*cos(abs(alpha)))-(I6y/4)*ddalpha/4)*10^(-3);% ������һ����
%%%%%%%%%%%%%%%%%%%%%%%%%
figure(14)  % ��λ�� mg.mm*rad^2.*s^-2=10^(-9) kg.m*s^-2=10^(-9)N=10^(6)uN;   %��λ��uN
F_yadd1_tran=-(I5y*(domega_z+omega_x.*omega_y))*10^(-3);  % ����������������ƽ�������������
F_yadd1_rot=(I6y*domega_x)*10^(-3);                                         % ����������������ת�������������
plot(t/T,F_yadd1_tran,'r-.',t/T,F_yadd1_rot,'g:',t/T,F_yadd1,'k-','LineWidth',2) 
ylabel('��������������F_{norm}(t) (uN)'); 
legend('F_{y,add1,tran}(t)','F_{y,add1,rot}(t)','F_{y,add1}(t)');  
title('��������������������F_{y,add1}(t))��ʱ��ı仯')  
grid on
axis([0.9,4.05,-inf,inf])
set(gca,'XTick',(0.9:0.1:4.05))
%%%%%%%%%%%%%%%%%%%%%%%%%
figure(15)
hold on
plot(t/T,F_yadd1,'r-',t/T,F_yadd2,'b-','LineWidth',2)    % ע�������domega_x=ddpsi��ddalpha�ķ��Ŵ��ڲ�ͬŶ
xlabel('\itNormalized time')
ylabel('\itF_{y,add} (uN)')
legend('F_{y,add1}','F_{y,add2}')
title('������������(����)��ʱ��ı仯����')  
grid on
axis([0.9,4.05,-inf,inf])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ����������������
% % �������ϵ�µĽ��ٶ�
% omega_x=dpsi;                     % չ��
% omega_y=dphi.*sin(psi);       % ����(��ʼ����)
% omega_z=dphi.*cos(psi);      % ������
% omega_h=dphi;       % �����Ľ��ٶ�% omega_h=-sqrt(omega_y^2+omega_z^2);  % ���ַ���˳ʱ��
% % �������ϵ�µĽǼ��ٶȡ����������������ļ���
% domega_x=ddpsi;
% % domega_y=ddphi.*sin(psi)+dphi.*dpsi.*cos(psi);
% domega_z=ddphi.*cos(psi)-dphi.*dpsi.*sin(psi); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m_wing=2.4*10^-9;                    % kg
x_com=1.920243385*10^-3;      %  m                           % ���ĵ�չ������
% z_com=-0.149785466+0.636=0.486215*10^-3;      % ��Ťת����������
z_com=0.149785466*10^-3;       %  m                          % ��Ťת����������
F_inert_y=-m_wing*(domega_z*x_com-domega_x*z_com+omega_y.*(omega_x*x_com+omega_z*z_com))*10^6;  % ������������������
F_inert_z=-m_wing*(-domega_y*x_com-omega_y.^2*z_com+omega_x.*(omega_z*x_com-omega_x*z_com))*10^6; % ������������������
% F_inert_y=-m_wing*(ddphi.*cos(psi)*x_com-(ddpsi-dphi.^2.*sin(2*psi)/2)*z_com)*10^6;       % ������������������
% F_inert_z=-m_wing*(-ddphi.*sin(psi)*x_com-(dpsi.^2+dphi.^2.*(sin(psi)).^2)*z_com)*10^6;  % ������������������
figure(16) 
plot(t/T,F_inert_y,'r-',t/T,F_inert_z,'g-','LineWidth',2) 
ylabel('��������������������F_{inert,y}(t) & F_{inert,z}(t) (uN)'); 
legend('F_{inert,y}(t)','F_{inert,z}(t)');  
title('��������������������F_{inert,y}(t) & F_{inert,z}(t)��ʱ��ı仯')  
grid on
axis([0.9,4.05,-inf,inf])
set(gca,'XTick',(0.9:0.1:4.05))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ������������ ������ϵ�£���(�ϳ�)������������ƽ������������+��ת����������+�����������������õ㣺ƬԪ�ҳ��е�
F_N=F_ytran+F_yadd1+F_yrot+F_inert_y;    % ��λ��rad^2*s^-2*kg*m=N=10^6uN
figure(17)
plot(t/T,F_N,'k-','LineWidth',2)
xlabel('\itNormalized time')
ylabel('\itF_N (uN)')
legend('F_N')
title('ƽ������������+��ת����������+������������(����)��ʱ��ı仯����')
grid on
axis([0.9,4.05,-inf,inf])
hold on
L=length(t);
plot([0,t(L)/T],[0,0],'k-','LineWidth',2);     %��x-axis
%% ������ϵ�£�ƬԪ��������ֽ⵽ƬԪ��ֱ�����ˮƽ�������õ㣺ƬԪ�ҳ��е�
F_vertical=-sign(alpha).*F_N.*cos(alpha);                % ��ֱ����vertical, %��λ��uN  %  
F_horizontal=F_N.*sin(abs(alpha));                          % ˮƽ����horizontal,%��λ��uN
% % F_horizontal=-sign(alpha).*F_N.*sin(alpha);     % ˮƽ����horizontal,%��λ��uN
% F_vertical=abs(F_N.*cos(alpha));                          % ��ֱ����vertical, %��λ��uN
% F_horizontal=abs(sign(alpha).*F_N.*sin(alpha));  % ˮƽ����horizontal,%��λ��uN
F_verticalaver=trapz(t,F_vertical)/(3*T)               % F_verticalaver=12.3074uN;   
g_uN=9.821;   % �����������ٶ�:g=9.821N/kg=9.821*10^6/10^6=9.821 uN/mg  ����g�Ĺ��ʵ�λ��m*s^-2��N*kg^-1
F_vaver=trapz(t,F_vertical)/(3*T)/g_uN               % F_vaver =1.2532; %��λ��mg      1.2532(��������)
F_haver=trapz(t,F_horizontal)/(3*T)/g_uN          % F_haver =-0.1109; %��λ��mg     -0.1109(��������)
F_haverabs=trapz(t,abs(F_horizontal))/(3*T);    
F_v2haver=F_vaver/F_haverabs                % ������������ֵ: F_v2haver =0.1150;  0.1045(��������)������������ϵ���ı�ֵ��û��ʵ�����壿
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(18)
hold on
plot(t/T,F_vertical,'r-',t/T,F_horizontal,'b:','LineWidth',2)
xlabel('\itNormalized time')
ylabel('\itF_{vertical} & F_{horizontal} (uN)')
legend('F_{vertical}','F_{horizontal}')
title('��������ֽ⵽��ֱ�����ˮƽ����ķ�������ʱ��ı仯����')   
grid on
axis([0.9,4.05,-inf,inf])
set(gca,'XTick',(0.9:0.1:4.05))
hold on
L=length(t);
plot([0,t(L)/T],[0,0],'k-','LineWidth',2);     %��x-axis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ������ϵ�£� �ɷ���������ϵ�����: ��ֱ�����ˮƽ����������ϵ��
% Coeff_liftdragF_N=0.0068201������λ��:mg*mm^-3*mm*mm^3= 10^(-9) kg*m
% Rou=1.225*10^(-3);           %��λ��Kg/m^3=10^6/(10^3)^3=10^(-3)mg/mm^3
% Coeff_liftdragF_N=(1/2)*Rou*C_avereff*R_wingeff^3*F_nd*10^(-12);  % kg*m
C_N_total=F_N*10^(-6)./(V_ref_aver.^2*Coeff_liftdragF_N*10^(-9));          % ��λ: N/(rad^2*s^-2*kg*m)=һ������
% C_v=F_vertical*10^(-6)./(V_ref_aver.^2*Coeff_liftdragF_N*10^(-9)); 
% C_h=F_horizontal*10^(-6)./(V_ref_aver.^2*Coeff_liftdragF_N*10^(-9));    % *10^3); 
C_v=-sign(alpha).*C_N_total.*cos(alpha); 
C_h=C_N_total.*sin(abs(alpha)); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%��˲ʱ������ϵ��ʹ�����λ��ֺ���trapz������ֵ���֣���������, ���ƽ����ֱ�����ˮƽ����������ϵ��
C_vaver=trapz(t,C_v)/(3*T);                 % C_vaver =2.3852;       2.4003(��������)   
C_haver=trapz(t,C_h)/(3*T);                 % C_haver =-0.2655;   -0.2125(��������)  
C_habsaver=trapz(t,abs(C_h))/(3*T);     
C_v2haver=C_vaver/C_habsaver;           % ��ֱ����ϵ����ˮƽ����ϵ���ı�ֵ: C_v2haver = 1.1298;  1.026(��������)
figure(19)
plot(t/T,C_N_total,'g-','LineWidth',2)
xlabel('\itNormalized time')
ylabel('\itC_{N,total}')
legend('C_{N,total}')
title('�����������ϵ����ʱ��ı仯����')
grid on
axis([0.9,4.05,-inf,inf])

figure(20)
plot(t/T,C_v,'r-',t/T,C_h,'b-','LineWidth',2)
xlabel('\itNormalized time')
ylabel('\itC_v & C_h')
legend('C_v','C_h')
title('��ֱ�����ˮƽ����������ϵ����ʱ��ı仯����')
grid on
axis([0.9,4.05,-inf,inf])
% close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

