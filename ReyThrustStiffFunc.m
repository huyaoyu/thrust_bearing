function [p,pz,A,DIA_IN,DIA_OUT,idx_boundary_in,idx_boundary_out,fx,fz,Kzz] = ReyThrustStiffFunc(...
PAD_DIM,TH_DIM,RA_DIM,ns,es,Dt,Dr,HP,AS,VISCO,VIS_EN,TURB_SWITCH,RHO,PB,ALPHA,DISP_PREFIX,ANG_OFF)

% ���Կ���
IS_DEBUG = 1;

% ���㳣��
ELE_NODES_NUM = 4;        % ÿ����Ԫ�Ľڵ�������Ľڵ�

STIFFNESS_BASE = 1e9;     % ���ڽ���õĸն�ֵ��10^9��ʽ��ʾ
% DAMPING_BASE   = 1e6;     % ���ڽ���õ�����ֵ��10^6��ʽ��ʾ

% ���һ����̬ѹ����
STABLE_FIELD_DISP_PREFIX = [DISP_PREFIX,'Stable pressure field: '];
[p,A,DIA_IN,DIA_OUT,idx_boundary_in,idx_boundary_out,fx,fz] = ReyThrustFunc(...
PAD_DIM,TH_DIM,RA_DIM,ns,es,Dt,Dr,HP,AS,VISCO,VIS_EN,TURB_SWITCH,RHO,PB,ALPHA,0,STABLE_FIELD_DISP_PREFIX,ANG_OFF);

% ========================= ��װ�����Ҷ� ==============================

% ת������Ԫ��Ԫ���ֲ���
Dtn = 1/Dt;
Drn = 1/Dr;
det_J = Dt*Dr;

% ���ڱ�ʾ����ı���
a   = -1;
atd = 1;

tRe = RHO*AS/(VISCO*VIS_EN); % ��������ʱ�ֲ���ŵ��������ʱ����

% �����Ҷ���ʱ������������
rhs_idx = zeros(1,4); % �����Ҷ�����

% ���򣨵�һ���꣩�;��򣨵ڶ����꣩�ķ���RHS��ʱ����
tRHSAr = -1*a/(VISCO*VIS_EN) * Drn^2 * det_J * [1/3,-1/3,1/3,-1/3]';
tRHSAt = -1*a/(VISCO*VIS_EN) * Dtn^2 * det_J * [1/3,-1/3,1/3,-1/3]';
tRHSB  = -1*a/(VISCO*VIS_EN) * Dtn^2 * det_J * [-1,1,1,-1]'; % t����
tRHSC  = -1*a/(VISCO*VIS_EN) * Drn^2 * det_J * [-1,-1,1,1]'; % r����

% ����������ʱ�����Ҷ�
% tRHSzt = atd * det_J * [1,1,1,1]';

% RHS
RHSz  = zeros((TH_DIM+1)*(RA_DIM+1),1);
% RHSzt = RHSz;

% ֹ���߿����ʼ
ths = PAD_DIM(1,1);
thend = PAD_DIM(1,2);

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

% ��ʾ��Ϣ
disp([DISP_PREFIX,'Assemble RHSs of perturbation equations...']);

for I = 1:1:TH_DIM*RA_DIM
    % ��װÿһ����Ԫ
    
    % ȡ���ĸ��ڵ�
    for J = 1:1:ELE_NODES_NUM
        tn(J,:) = [es(I,J),ns(es(I,J),1),ns(es(I,J),2)];
    end % J
    
    % ���㵥Ԫ����������, theta center��axis center
    thc = (tn(1,2) + tn(2,2))/2 + ANG_OFF;
    rac = (tn(1,3) + tn(4,3))/2;
    
    % ����ҺĤ���
    h = HP + ALPHA*rac*sin(thend - thc);
    
    % ����Բ�ܷ����ϵ��Gx��Gz
    % ������е�����
%     if(TURB_SWITCH == 0)
%         Gx = 1 / (4 * rac);
%         Gz = rac / 4;
%     elseif(TURB_SWITCH == 1)
%         lRe = tRe * h * rac;
%         if(lRe > 2000)
%             Gx = 1 / ((12+0.0136*lRe^0.9) * rac / 3);
%             Gz = rac / ((12+0.0043*lRe^0.96) / 3);
%         else
%             Gx = 1 / (4 * rac);
%             Gz = rac / 4; 
%         end % lRe > 2000
%     end
    % ===== 2012.3.31 =========
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
    RHSz(rhs_idx,1) = RHSz(rhs_idx,1) + trRHS + ttRHS;
%     RHSzt(rhs_idx,1) = RHSzt(rhs_idx,1) + tRHSzt * rac;
end

% �������Ҷˣ��ر�ע��߽磡������δ�޸�
RHSz( idx_boundary_in) = PB*DIA_IN;
RHSz(idx_boundary_out) = PB*DIA_OUT;

% RHSzt( idx_boundary_in) = PB*DIA_IN;
% RHSzt(idx_boundary_out) = PB*DIA_OUT;

% ��ʾ��Ϣ
disp([DISP_PREFIX,'RHSs assembled.']);

% =========================== Solve ==================================
% % z������ٶ��㶯ѹ����
% pzt = A\RHSzt;
% disp('pzt solved.');

% z�����λ���㶯ѹ����
pz  = A\RHSz;
disp([DISP_PREFIX,'pz solved.']);

% �鿴���
res_vec = A*pz - RHSz;
res_norm = norm(res_vec);

disp([DISP_PREFIX,'The norm of residual vector of pz is ',num2str(res_norm),'.']);

% res_vec = A*pzt - RHSzt;
% res_norm = norm(res_vec);
% 
% disp(['The norm of residual vector of pzt is ',num2str(res_norm),'.']);

if(IS_DEBUG)
    % ���ֲ���ŵ����ֵ
    disp([DISP_PREFIX,'The maximun Reynolds number is ',...
        num2str(max(check_reynolds)),'.']);
    
    disp([DISP_PREFIX,'The maximun r number is ',...
        num2str(max(check_rac)),'.']);
    
    disp([DISP_PREFIX,'The maximun h number is ',...
        num2str(max(check_h)),'.']);

    % ������ʱ����
    clear check_reynolds
    clear check_rac
    clear check_h
end % IS_DEBUG

% ============================= ��ʾ��� ==============================

% �ع���������γ�һ����������Ϊ��TH_DIM+1���У���RA_DIM+1����
pzr  = reshape(pz,TH_DIM+1,RA_DIM+1);  % λ���㶯ѹ�����ع�
% pztr = reshape(pzt,TH_DIM+1,RA_DIM+1); % �ٶ��㶯ѹ�����ع�
% th_idx_re = ((1:1:(TH_DIM+1))-1) .* Dt ./ (2*pi) .* 360; % re for reconstruction
% ra_idx_re = ((1:1:(RA_DIM+1))-1) .* Dr + PAD_DIM(1,3);     % re for reconstruction
% 
% disp([DISP_PREFIX,'Display the results in 3D.']);
% 
% figure
% 
% surf(ra_idx_re,th_idx_re,pzr,'LineStyle','none')

% figure
% 
% surf(ra_idx_re,th_idx_re,pztr,'LineStyle','none')

% ========================= ���㶯ѹ��������� ===========================
% ���λ���
Kzz = 0; % ���غ�ƽ�еĸն�
% Czz = 0; % ���غ�ƽ�е�����
Kzz_temp = 0; % ��ʱ����
% Czz_temp = 0; % ��ʱ����
% fth_begin = -1*Dt/2; % ��ʱ������force theta
% fth = 0; % ��ʱ�������Ƕ�λ�ã���λΪrad
fa = PAD_DIM(1,3);

for I = 1:1:RA_DIM
%     fth = fth_begin;
    for J = 1:1:TH_DIM
%         fth = fth + Dt;
        
%         if(J ~= TH_DIM)
        if(1)
            % δ���ػ��߽�
            Kzz_temp = pzr(  J,I)*fa + pzr(  J,I+1)*(fa+Dr) + ...
                       pzr(J+1,I)*fa + pzr(J+1,I+1)*(fa+Dr);
%             Czz_temp = pztr(  J,I) + pztr(  J,I+1) + ...
%                        pztr(J+1,I) + pztr(J+1,I+1);
        else
            % ����ػ��߽�
            Kzz_temp = pzr(J,I) + pzr(J,I+1) + ...
                       pzr(1,I) + pzr(1,I+1);
%             Czz_temp = pztr(J,I) + pztr(J,I+1) + ...
%                        pztr(1,I) + pztr(1,I+1);
        end
        

       Kzz = Kzz + (-1)*Kzz_temp;
%        Czz = Czz + (-1)*Czz_temp;

    end % J
    fa = fa + Dr;
end % I

% ���Ի���ϵ��
Kzz = Kzz * (Dt * Dr)/4;
% Czz = Czz * (Dt * Dr)/4;

% ��ʾ�ն���Ϣ
disp([DISP_PREFIX,'Kzz = ',num2str(Kzz/STIFFNESS_BASE),'x10^9 N/m']);
% disp(['Czz = ',num2str(Czz/DAMPING_BASE),'x10^6 Ns/m']);

if(IS_DEBUG)
    disp([DISP_PREFIX,'Debug mode is on.']);
end