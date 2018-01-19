%
% Description
% ===========
%
% ���ļ����ڲ��Լ���ֹ����е���������ϵ�����������ѧģ��������Gunhee Jang
% �����ס�Determination of the dynamic coefficients of the coupled journal
% and thrust bearings by the perturbation method����������û�п��Ǵ��ȡ�����
% �п�����ֹ����е��������ã������������õ���Eskild Storteig���˵����ף���
% Dynamic characteristics of hydrodynamically lubricated fixed-pad thrust
% bearings��.
%
% ��ⷽ��ʱ�����˶�ά����Ԫ����ʹ���Ľڵ�Ȳ�����Ԫ����������������̬ƽ��λ��
% ����⣬�������������̫Reynolds���̵��ӳ���ReyThrustFunc���ĳ��������һ��
% ��������λ�õ�ֹ����е���Reynolds���̡�
%
% ��⾲ƽ��λ������Newton-Raphson������������ʽ�ο�������ѧ����Ƽ��ʿѧλ����
% ����������ж���ѧ��ģ�����������о���һ���й��ھ������ƽ��λ�õ������̡�
%
% - 3��31�գ����������ڻ�������˳�������������ѧģ���Ƶ�����������������⡣
% �������е��ķ������a����Ϊa��atd�����������ֱ����ڼ���Kzz��Czz��
% 
% ע�⣬�ڻ���ʱ�����Ǳ߽�Ļػ���
%
% - 4��1�գ������NR������Ϻ���ʾ���յ�ѹ���ֲ������Ӽ���ʱ����ԡ�
% - 4��9�գ�������ת�ǸնȺ�������㡣
%
%
% Author
% ======
%
% Yaoyu HU <huyaoyu@sjtu.edu.cn>
%
% Date
% ====
%
% ������2012��3��30��
% �޸ģ�2012��3��31��
% �޸ģ�2012��4��1��
% �޸ģ�2012��4��9��




% ========================= �������ռ� ==================================

close ALL
clear
clc

disp('Reynolds Equation Solver, for thrust bearing coefficients, with FEM');

ThrustInput;

% Newton-Raphson ����

qk  = HP;
qk1 = 0; % ��һ��

% ������Ϣǰ׺
iter_count = 1;
iter_prefix = '';

% �����������ɢ����������Ԫ����
TH_DIM = 100;     % Բ�ܷ������ɢ��ע��������ָ�ĵ�Ԫ��
RA_DIM = 100;      % �������ɢ��ע��������ָ�ĵ�Ԫ��

% ===== 2012.3.24 ======
% ʹ�ô���ػ����������ɺ���
% ������Ϣ
% [ns,es] = RecField2DIso([0,0],[2*pi,0.1],[TH_DIM,RA_DIM]);
% [ns,es] = RecField2DIsoWrap([0,0],[2*pi,0.1],[TH_DIM,RA_DIM],[1,0]);
% ===== 2012.3.31 =======
% ע����ʼλ��
[ns,es] = RecField2DIsoWrap(...
    [PAD_DIM(1,1),PAD_DIM(1,3)],...
    [PAD_DIM(1,2)-PAD_DIM(1,1),PAD_DIM(1,4)-PAD_DIM(1,3)],...
    [TH_DIM,RA_DIM],[0,0]);

% ======== 2012.3.26 =========
% Jacobian ����ĶԽ�Ԫ�أ���ʱ����
% ����Ŀǰ�̶�������
Dt = (PAD_DIM(1,2)-PAD_DIM(1,1))/TH_DIM; % ��Ԫ�ĵ�һ���곤��
Dr = (PAD_DIM(1,4)-PAD_DIM(1,3))/RA_DIM; % ��Ԫ�ĵڶ����곤��

while(abs(qk-qk1) > NEWTON_RAPHSON_NORM)
    % ���һ������ն�
    
    % ��Ϣǰ׺
    iter_prefix = ['NR ',num2str(iter_count),': '];
    
    [p,pz,A,DIA_IN,DIA_OUT,idx_boundary_in,idx_boundary_out,fx,fz,Kzz] = ReyThrustStiffFunc(...
    PAD_DIM,TH_DIM,RA_DIM,ns,es,Dt,Dr,qk,AS,VISCO,VIS_EN,TURB_SWITCH,RHO,PB,ALPHA,iter_prefix,0);

    qa  = -1 * (fz - W) / (-1 * Kzz);
    qk1 = qk;
    qk  = qk + qa;
    
    iter_count = iter_count + 1;
end % abs(qk-qk1) < NEWTON_RAPHSON_NORM

% ========================= ��ʾѹ��������� ==========================

% �ع���������γ�һ����������Ϊ��TH_DIM+1���У���RA_DIM+1����
pr = reshape(p,TH_DIM+1,RA_DIM+1); % ѹ�����ع�
th_idx_re = ((1:1:(TH_DIM+1))-1) .* Dt ./ (2*pi) .* 360; % re for reconstruction
ra_idx_re = ((1:1:(RA_DIM+1))-1) .* Dr + PAD_DIM(1,3);     % re for reconstruction

disp('Display the final pressure field in 3D.');

surf(ra_idx_re,th_idx_re,pr,'LineStyle','none')

% ========================= ��װ�����Ҷ� ==============================

disp('Solve the damping coeffiecents ...');

% ת������Ԫ��Ԫ���ֲ���
Dtn = 1/Dt;
Drn = 1/Dr;
det_J = Dt*Dr;

% ���ڱ�ʾ����ı���
% a    = 1;
% atx  = -1; % x�ڽǵ�alphaֵ
% aty  = 1; % y�ڽǵ�alphaֵ
a    = 1;
atx  = -1; % x�ڽǵ�alphaֵ
aty  = 1; % y�ڽǵ�alphaֵ
atd  = 1;
atxd = -1;
atyd = 1;
ang_off = 0; % ���������ϵ�ĽǶ�ƫ��

tRe = RHO*AS/(VISCO*VIS_EN); % ��������ʱ�ֲ���ŵ��������ʱ����

% �����Ҷ���ʱ������������
rhs_idx = zeros(1,4); % �����Ҷ�����

% ����������ʱ�����Ҷ�
tRHSzt  =  atd * det_J * [1,1,1,1]';
tRHStxt = atxd * det_J * [1,1,1,1]';
tRHStyt = atyd * det_J * [1,1,1,1]';

% RHS
RHSzt = zeros((TH_DIM+1)*(RA_DIM+1),1);
% ======= 2012.4.9����ʼ�� ==========
th_x = 0; % x �����ת��
th_y = 0; % y �����ת��
cz0 = 0; % center of mass �� z ����
rz  = 0; % �ڶ����ĵ�z����
tphi1 = (cz0 - rz)*sin(th_x);
tphi2 = -1*cos(th_x);
tpsi1 = (cz0 - rz)*sin(th_y);
tpsi2 = cos(th_y);
RHStx = RHSzt; % x �����ת�ǣ�theta x
RHSty = RHSzt; % y �����ת�ǣ�theta y
% tRHSas = -1*0.5 * AS * [-1,1,1,-1]' .*...
%     [Dtn,Dtn,Dtn,Dtn]' * det_J; % 'as' for angluar speed
tRHSas = 0.5 * AS * [-1,1,1,-1]' .*...
    [Dtn,Dtn,Dtn,Dtn]' * det_J; % 'as' for angluar speed
RHStxt = RHSzt;
RHStyt = RHSzt;

% ���򣨵�һ���꣩�;��򣨵ڶ����꣩�ķ���RHS��ʱ����
% tRHSAr = (VISCO*VIS_EN) * Drn^2 * det_J * [1/3,-1/3,1/3,-1/3]';
% tRHSAt = (VISCO*VIS_EN) * Dtn^2 * det_J * [1/3,-1/3,1/3,-1/3]';
% tRHSB  = (VISCO*VIS_EN) * Dtn^2 * det_J * [-1,1,1,-1]'; % t����
% tRHSC  = (VISCO*VIS_EN) * Drn^2 * det_J * [-1,-1,1,1]'; % r����
tRHSAr = -1*(VISCO*VIS_EN) * Drn^2 * det_J * [1/3,-1/3,1/3,-1/3]';
tRHSAt = -1*(VISCO*VIS_EN) * Dtn^2 * det_J * [1/3,-1/3,1/3,-1/3]';
tRHSB  = -1*(VISCO*VIS_EN) * Dtn^2 * det_J * [-1,1,1,-1]'; % t����
tRHSC  = -1*(VISCO*VIS_EN) * Drn^2 * det_J * [-1,-1,1,1]'; % r����


% ֹ���������ϵ�������õĲ���
MS  = -0.25;
NS  = 0.066;
tGx = 2^(1+MS) / (NS * (2 + MS));
tGz = 2^(1+MS) / NS;
if(IS_DEBUG)
    check_reynolds = zeros(1,TH_DIM*RA_DIM); % ��¼�ֲ���ŵ�����ڵ���
    check_rac      = zeros(1,TH_DIM*RA_DIM); % ��¼��Ԫ����λ�����ڵ���
    check_h        = zeros(1,TH_DIM*RA_DIM); % ��¼ҺĤ������ڵ���
end

% ======= 2012.4.9�������� ==========


% ֹ���߿����ʼ
ths = PAD_DIM(1,1);
thend = PAD_DIM(1,2);

% ��ʾ��Ϣ
disp('Assemble RHSs of perturbation equations...');

for I = 1:1:TH_DIM*RA_DIM
    % ��װÿһ����Ԫ
    
    % ȡ���ĸ��ڵ�
    for J = 1:1:ELE_NODES_NUM
        tn(J,:) = [es(I,J),ns(es(I,J),1),ns(es(I,J),2)];
    end % J
    
    % ���㵥Ԫ����������, theta center��axis center
    thc = (tn(1,2) + tn(2,2))/2 + ang_off;
    rac = (tn(1,3) + tn(4,3))/2;
    
    % ����ҺĤ���
    h = qk1 + ALPHA*rac*sin(thend - thc);
    
    % ����phi��psi
    phi_t = tphi1 + rac * tphi2 * sin(thc);
    psi_t = tpsi1 + rac * tpsi2 * cos(thc);
    
    % ========== 2012.4.9����ʼ�� ========
    % ֹ����е�����
    if(TURB_SWITCH == 0)
        Gx = 1 / (4 * rac);
        Gz = rac / 4;
    elseif(TURB_SWITCH == 1)
        lRe = tRe * h * rac;
        if(IS_DEBUG)
            check_reynolds(1,I) = lRe;
            check_rac(1,I) = rac;
            check_h(1,I) = h;
        end
        gRe = lRe^(1+MS);
        if(lRe > 2000)
%             Gx = 1 / ((12+0.0136*lRe^0.9) * rac / 3);
            Gx = tGx / gRe / rac * 3;
%             Gz = rac / ((12+0.0043*lRe^0.96) / 3);
            Gz = tGz / gRe * rac * 3;
        else
            Gx = 1 / (4 * rac);
            Gz = rac / 4; 
        end % lRe > 2000
    end
    
    % ��ʱ�ڵ�ѹ��
    p0t = p(tn(:,1),1);
    
    p0A = sum(0.25 * [1,-1,1,-1]' .* p0t);
    p0B = sum(0.25 * [-1,1,1,-1]' .* p0t);
    p0C = sum(0.25 * [-1,-1,1,1]' .* p0t);
    
    ttRHS = (tRHSAt .* p0A .* h^2 + tRHSB .* p0B .* h^2) * Gx;
    trRHS = (tRHSAr .* p0A .* h^2 + tRHSC .* p0C .* h^2) * Gz;

  
    % ���㷽���Ҷ�
    rhs_idx = [tn(1,1),tn(2,1),tn(3,1),tn(4,1)];
    RHSReadyToAdd  = trRHS + ttRHS + tRHSas .* rac;
    RHSReadyToAddX = atx * RHSReadyToAdd;
    RHSReadyToAddY = aty * RHSReadyToAdd;
    RHStx(rhs_idx,1) = RHStx(rhs_idx,1) + RHSReadyToAddX * phi_t;
    RHSty(rhs_idx,1) = RHSty(rhs_idx,1) + RHSReadyToAddY * psi_t;
    RHStxt(rhs_idx,1) = RHStxt(rhs_idx,1) + tRHStxt * rac * phi_t;
    RHStyt(rhs_idx,1) = RHStyt(rhs_idx,1) + tRHStyt * rac * psi_t;
    % ========== 2012.4.9�������� ========
      
    % ���㷽���Ҷ�
    RHSzt(rhs_idx,1) = RHSzt(rhs_idx,1) + tRHSzt * rac;
end

% �������Ҷˣ��ر�ע��߽磡������δ�޸�
RHSzt( idx_boundary_in) = PB*DIA_IN;
RHSzt(idx_boundary_out) = PB*DIA_OUT;

% ========== 2012.4.9����ʼ�� ========
RHStx( idx_boundary_in) = PB*DIA_IN;
RHStx(idx_boundary_out) = PB*DIA_OUT;

RHSty( idx_boundary_in) = PB*DIA_IN;
RHSty(idx_boundary_out) = PB*DIA_OUT;

RHStxt( idx_boundary_in) = PB*DIA_IN;
RHStxt(idx_boundary_out) = PB*DIA_OUT;

RHStyt( idx_boundary_in) = PB*DIA_IN;
RHStyt(idx_boundary_out) = PB*DIA_OUT;
% ========== 2012.4.9�������� ========

% ��ʾ��Ϣ
disp('RHSs assembled.');

% =========================== Solve ==================================

% z������ٶ��㶯ѹ����
pzt = A\RHSzt;
disp('pzt solved.');

res_vec = A*pzt - RHSzt;
res_norm = norm(res_vec);

disp(['The norm of residual vector of pzt is ',num2str(res_norm),'.']);

% theta x ������㶯ѹ����
ptx = A\RHStx;
disp('ptx solved.');

res_vec = A*ptx - RHStx;
res_norm = norm(res_vec);

disp(['The norm of residual vector of ptx is ',num2str(res_norm),'.']);
% theta y ������㶯ѹ����
pty = A\RHSty;
disp('pty solved.');

res_vec = A*pty - RHSty;
res_norm = norm(res_vec);

disp(['The norm of residual vector of pty is ',num2str(res_norm),'.']);

% theta x ������ٶ��㶯ѹ����
ptxt = A\RHStxt;
disp('pty solved.');

res_vec = A*ptxt - RHStxt;
res_norm = norm(res_vec);

disp(['The norm of residual vector of ptxt is ',num2str(res_norm),'.']);

% theta y ������ٶ��㶯ѹ����
ptyt = A\RHStyt;
disp('pty solved.');

res_vec = A*ptyt - RHStyt;
res_norm = norm(res_vec);

disp(['The norm of residual vector of ptyt is ',num2str(res_norm),'.']);

% ============================= ��ʾ��� ==============================

% �ع���������γ�һ����������Ϊ��TH_DIM+1���У���RA_DIM+1����
% z ��������
pztr = reshape(pzt,TH_DIM+1,RA_DIM+1); % �ٶ��㶯ѹ�����ع�
pzr = reshape(pz,TH_DIM+1,RA_DIM+1); % �ٶ��㶯ѹ�����ع�

disp('Display the results in 3D.');

figure

surf(ra_idx_re,th_idx_re,pztr,'LineStyle','none')

view(0,90);
colorbar;
title('pressure field of velocity perturbation along z');

% theta x ����ն�
ptxr = reshape(ptx,TH_DIM+1,RA_DIM+1); % �ٶ��㶯ѹ�����ع�

disp('Display the displacement perturbation of theta x in 3D.');

figure

surf(ra_idx_re,th_idx_re,ptxr,'LineStyle','none')
view(0,90);
colorbar;
title('pressure field of a-displacement perturbation along theta x');

% theta y ����ն�
ptyr = reshape(pty,TH_DIM+1,RA_DIM+1); % �ٶ��㶯ѹ�����ع�

disp('Display the displacement perturbation of theta y in 3D.');

figure

surf(ra_idx_re,th_idx_re,ptyr,'LineStyle','none')
view(0,90);
colorbar;
title('pressure field of a-displacement perturbation along theta y');

% theta x ����ն�
ptxtr = reshape(ptxt,TH_DIM+1,RA_DIM+1); % �ٶ��㶯ѹ�����ع�

disp('Display the velocity perturbation of theta x in 3D.');

figure

surf(ra_idx_re,th_idx_re,ptxtr,'LineStyle','none')
view(0,90);
colorbar;
title('pressure field of a-velocity perturbation along theta x');

% theta y ����ն�
ptytr = reshape(ptyt,TH_DIM+1,RA_DIM+1); % �ٶ��㶯ѹ�����ع�

disp('Display the velocity perturbation of theta y in 3D.');

figure

surf(ra_idx_re,th_idx_re,ptytr,'LineStyle','none')
view(0,90);
colorbar;
title('pressure field of a-velocity perturbation along theta y');

% ========================= ���㶯ѹ��������� ===========================
% % ���λ���
% Czz = 0; % ���غ�ƽ�е�����
% Czz_temp = 0; % ��ʱ����
% fth_begin = -1*Dt/2; % ��ʱ������force theta
% fth = 0; % ��ʱ�������Ƕ�λ�ã���λΪrad
% 
% for I = 1:1:RA_DIM
%     fth = fth_begin;
%     for J = 1:1:TH_DIM
%         fth = fth + Dt;
%         
% %         if(J ~= TH_DIM)
%         if(1)
%             % δ���ػ��߽�
%             Czz_temp = pztr(  J,I) + pztr(  J,I+1) + ...
%                        pztr(J+1,I) + pztr(J+1,I+1);
%         else
%             % ����ػ��߽�
%             Czz_temp = pztr(J,I) + pztr(J,I+1) + ...
%                        pztr(1,I) + pztr(1,I+1);
%         end
%         
%         Czz = Czz + (-1)*Czz_temp;
% 
%     end % J
% end % I
% 
% % ���Ի���ϵ��
% Czz = Czz * (Dt * Dr)/4;

alpha = [0,0,a,atx,aty]';
K = DouIntegration(alpha,pzr,ptxr,ptyr,Dt,Dr,PAD_DIM(1,3));
C = DouIntegration(alpha,pztr,ptxtr,ptytr,Dt,Dr,PAD_DIM(1,3));

% ��ʾ�ն���Ϣ
disp(['Czz = ',num2str(C(3,3)/DAMPING_BASE),'x10^6 Ns/m']);

if(IS_DEBUG)
    disp('Debug mode is on.');
end
