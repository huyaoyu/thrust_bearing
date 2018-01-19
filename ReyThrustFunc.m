function [p,A,DIA_IN,DIA_OUT,idx_boundary_in,idx_boundary_out,fx,fz] =...
    ReyThrustFunc(PAD_DIM,TH_DIM,RA_DIM,ns,es,Dt,Da,HP,AS,VISCO,VIS_EN,TURB_SWITCH,RHO,PB,ALPHA,coeff_check,DISP_PREFIX,ANG_OFF)


% ��ʾ�����װ���ȵ���ֵ
MAX_NONDISP_DIM = 100*100;

% ϡ�������ÿ�еķ���Ԫ�ظ���
ROW_NONZERO_NUM = 16;

ELE_NODES_NUM = 4; % ÿ����Ԫ�Ľڵ������������4

% ����ѹ���߽�ʱ���Խ����ϵĴ���
BIG_COEFF = 1e16;

% =========================== ����������ɢ ================================

% ��ʾ��Ϣ
disp([DISP_PREFIX,'Pre-allocate memory and generate the grid ...']);

% ϡ����������������,������
row_idx = zeros(1,TH_DIM*RA_DIM*ELE_NODES_NUM^2);

% ϡ������������������������
col_idx = zeros(1,TH_DIM*RA_DIM*ELE_NODES_NUM^2);

% ϡ�����ķ���ֵ������������
elems = zeros(1,TH_DIM*RA_DIM*ELE_NODES_NUM^2);
if(coeff_check == 1)
    elems_h0 = elems;
else
    elems_h0 = 0;
end

% ===== 2012.3.24 ======
% �����˵�һ���꣬����Ľڵ���
% �����Ҷˣ�RHS��������
% RHS = zeros((TH_DIM+1)*(RA_DIM+1),1);
% ===== 2012.3.28 ======
% �����˵�һ����
RHS = zeros((TH_DIM+1)*(RA_DIM+1),1);

% ���ɵĽ�
p = 0;

disp([DISP_PREFIX,'Memory allocated and grid generated.']);

% ========================= ��װϵ������ =================================

% ��ʾ��Ϣ
disp([DISP_PREFIX,'Begin assemble the system...']);

% ��ʱ�ڵ����,���д洢����һ���ǽڵ���
tn = zeros(4,3);

% ѭ����ʱ����
t_th = -1  / (VISCO*VIS_EN);  % �������ʱ����
t_ra = -1  / (VISCO*VIS_EN); % �������ʱ����
h = 0; % ҺĤ���
Ke = zeros(ELE_NODES_NUM,ELE_NODES_NUM); % ��Ԫ�ĸնȾ���
OT  =  1/3; % one thirds
MOT = -1/3; % negtive one thirds, ��ΪNOT������ؼ����ظ�������MOT
OS  =  1/6; % one sixths
NOS = -1/6; % negtive one sixths
idx_offset = 0; % ȫ�����������ĵ�ǰ��Ԫƫ��
idx_off_base = ELE_NODES_NUM^2; % ȫ�����������ĵ�Ԫƫ�ƻ�
tRe = RHO*AS/(VISCO*VIS_EN); % ��������ʱ�ֲ���ŵ��������ʱ����

% % ======== 2012.3.26 =========
% % Jacobian ����ĶԽ�Ԫ�أ���ʱ����
% % ����Ŀǰ�̶�������
% Dt = (PAD_DIM(1,2)-PAD_DIM(1,1))/TH_DIM; % ��Ԫ�ĵ�һ���곤��
% Da = (PAD_DIM(1,4)-PAD_DIM(1,3))/RA_DIM; % ��Ԫ�ĵڶ����곤��
Dtn = 1/Dt;
Dan = 1/Da;
det_J = Dt*Da;
% ���¸�������ֵ
t_th = t_th*Dtn^2*det_J;
t_ra = t_ra*Dan^2*det_J;

% ��ʱ����ϵ������
MatrixTh = [ 
     OT,MOT,NOS, OS;
    MOT, OT, OS,NOS;
    NOS, OS, OT,MOT;
     OS,NOS,MOT, OT
];

% ��ʱ����ϵ������
MatrixRa = [
     OT, OS,NOS,MOT;
     OS, OT,MOT,NOS;
    NOS,MOT, OT, OS;
    MOT,NOS, OS, OT
];

% �����Ҷ���ʱ������������
tRHS = -1*0.5 * AS * [-1,1,1,-1]';
rhs_idx = zeros(1,4); % �����Ҷ�����
% ======== 2012.3.26 =========
% �����Ҷ˿���Jacobian���󣬵�һ���귽��
tRHS = tRHS .* [Dtn,Dtn,Dtn,Dtn]' * det_J; 

% ֹ���߿����ʼ
ths = PAD_DIM(1,1);
thend = PAD_DIM(1,2);

% ֹ���������ϵ�������õĲ���
MS  = -0.25;
NS  = 0.066;
tGx = 2^(1+MS) / (NS * (2 + MS));
tGz = 2^(1+MS) / NS;

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
%     h = HP + ALPHA*rac*sin(thc - ths);
    h = HP + ALPHA*rac*sin(thend - thc);
    
    h3 = h^3;
    % ����Բ�ܷ����ϵ��Gx��Gz
    if(TURB_SWITCH == 0)
        Gx = t_th * h3 / 12 / rac;
        Gz = t_ra * h3 / 12 * rac;
    elseif(TURB_SWITCH == 1)
        lRe = tRe * h * rac;
        gRe = lRe^(1+MS);
        if(lRe > 2000)
            % ====== 2012.3.31 =======
            % ����Ϊֹ����е�����ģ��
%             Gx = t_th * h3 / (12+0.0136*lRe^0.9) / rac;
%             Gz = t_ra * h3 / (12+0.0043*lRe^0.96) * rac;
            Gx = t_th * h3 * tGx / gRe / rac;
            Gz = t_ra * h3 * tGz / gRe * rac;
        else
            Gx = t_th * h3 / 12 / rac;
            Gz = t_ra * h3 / 12 * rac;
        end % lRe > 2000
    end   

    % ���������Ԫ�Ĺ���նȾ���
    Ke = Gx * MatrixTh + Gz * MatrixRa;
    % ������նȾ���������
    Ke_row = reshape(Ke,1,idx_off_base);
    
    % ��ֵ����նȾ�������
    row_idx(1,(idx_offset+1):1:(idx_offset+idx_off_base)) = [
        tn(1,1),tn(1,1),tn(1,1),tn(1,1),...
        tn(2,1),tn(2,1),tn(2,1),tn(2,1),...
        tn(3,1),tn(3,1),tn(3,1),tn(3,1),...
        tn(4,1),tn(4,1),tn(4,1),tn(4,1)
    ];

    col_idx(1,(idx_offset+1):1:(idx_offset+idx_off_base)) = [
        tn(1,1),tn(2,1),tn(3,1),tn(4,1),...
        tn(1,1),tn(2,1),tn(3,1),tn(4,1),...
        tn(1,1),tn(2,1),tn(3,1),tn(4,1),...
        tn(1,1),tn(2,1),tn(3,1),tn(4,1)
    ];
    
    elems(1,(idx_offset+1):1:(idx_offset+idx_off_base)) = Ke_row;
    
    if(coeff_check == 1)
        h0 = HP;
    
        h3 = h0^3;
        % ����Բ�ܷ����ϵ��Gx��Gz
        if(TURB_SWITCH == 0)
            Gx = t_th * h3 / 12 / rac;
            Gz = t_ra * h3 / 12 * rac;
        elseif(TURB_SWITCH == 1)
            lRe = tRe * h * rac;
            if(lRe > 2000)
                % ====== 2012.3.31 =======
                % ����Ϊֹ����е�����ģ��
%                 Gx = t_th * h3 / (12+0.0136*lRe^0.9) / rac;
%                 Gz = t_ra * h3 / (12+0.0043*lRe^0.96) * rac;
                Gx = tGx / gRe / rac;
                Gz = tGz / gRe * rac;
            else
                Gx = t_th * h3 / 12 / rac;
                Gz = t_ra * h3 / 12 * rac;
            end % lRe > 2000
        end   

        % ���������Ԫ�Ĺ���նȾ���
        Ke = Gx * MatrixTh + Gz * MatrixRa;
        % ������նȾ���������
        Ke_row = reshape(Ke,1,idx_off_base);

        elems_h0(1,(idx_offset+1):1:(idx_offset+idx_off_base)) = Ke_row;
    end

    idx_offset = idx_offset + idx_off_base;
    
    % ���㷽���Ҷ�
    rhs_idx = [tn(1,1),tn(2,1),tn(3,1),tn(4,1)];
    RHS(rhs_idx,1) = RHS(rhs_idx,1) + h .* tRHS .* rac;
end

% ����ϡ��ϵ������
A = sparse(row_idx,col_idx,elems,...
    (TH_DIM+1)*(RA_DIM+1),(TH_DIM+1)*(RA_DIM+1)...
    );
if(coeff_check == 1)
    A0 = sparse(row_idx,col_idx,elems_h0,...
        (TH_DIM+1)*(RA_DIM+1),(TH_DIM+1)*(RA_DIM+1)...
        );
end

% ��ʾ��Ϣ
disp([DISP_PREFIX,'Handle the boundaries...']);

% ===== 2012.3.24 ======
% ��С����Ľڵ���
% % ����߽�������������ѹ���߽�
% % ��������Ϊ������ֵʱ��������ѹ��ֵ����֪��
% % �ҵ�ѹ���߽�������
% idx_boundary_in  = 1:1:(TH_DIM+1);
% idx_boundary_out = ((TH_DIM+1)*RA_DIM+1):1:((TH_DIM+1)*(RA_DIM+1));
% % ���ϵ������ĶԽ���
% DIA_A = diag(A);
%
% ����߽�������������ѹ���߽�
% ��������Ϊ������ֵʱ��������ѹ��ֵ����֪��
% �ҵ�ѹ���߽�������
idx_boundary_in  = [1:1:(TH_DIM+1),1:(TH_DIM+1):((TH_DIM+1)*RA_DIM+1)];
idx_boundary_out = [((TH_DIM+1)*RA_DIM+1):1:((TH_DIM+1)*(RA_DIM+1)),...
    ( (TH_DIM+1):(TH_DIM+1):((TH_DIM+1)*(RA_DIM+1)) )];
% ���ϵ������ĶԽ���
DIA_A = diag(A);

if(coeff_check == 1)
    DIA_A0 = diag(A0);
end

% ��õ��������ڱ߽�ͳ��ڱ߽�ĶԽ���ϵ��
DIA_IN  = DIA_A( idx_boundary_in) + BIG_COEFF;
DIA_OUT = DIA_A(idx_boundary_out) + BIG_COEFF;

% �������Ҷˣ���ڱ߽�
RHS( idx_boundary_in) = PB*DIA_IN;
RHS(idx_boundary_out) = PB*DIA_OUT;

if(coeff_check == 1)
    DIA_IN  = DIA_A0( idx_boundary_in) + BIG_COEFF;
    DIA_OUT = DIA_A0(idx_boundary_out) + BIG_COEFF;  
end

% �ͷ�ϡ�������ڴ�
clear A
clear A0

% �ӳ���������������
final_row_idx = [row_idx,idx_boundary_in,idx_boundary_out];
final_col_idx = [col_idx,idx_boundary_in,idx_boundary_out];
final_elems   = [elems,ones(1,TH_DIM+1+RA_DIM+1)*BIG_COEFF,ones(1,TH_DIM+1+RA_DIM+1)*BIG_COEFF];

if(coeff_check == 1)
    final_elems_h0 = [elems_h0,ones(1,TH_DIM+1+RA_DIM+1)*BIG_COEFF,ones(1,TH_DIM+1+RA_DIM+1)*BIG_COEFF];
else
    final_elems_h0 = 0;
end

% ��������ϡ��ϵ������
A = sparse(final_row_idx,final_col_idx,final_elems,...
    (TH_DIM+1)*(RA_DIM+1),...
    (TH_DIM+1)*(RA_DIM+1));

if(coeff_check == 1)
A0 = sparse(final_row_idx,final_col_idx,final_elems_h0,...
    (TH_DIM+1)*(RA_DIM+1),...
    (TH_DIM+1)*(RA_DIM+1));
end

% �ͷŲ�����ʱ�ڴ�
clear row_idx
clear col_idx
clear elems
clear elems_h0
clear final_row_idx
clear final_col_idx
clear final_elems
clear final_elems_h0
% clear idx_boundary_in
% clear idx_boundary_out
clear DIA_A
% clear DIA_IN
% clear DIA_OUT

% ��ʾ��Ϣ
disp([DISP_PREFIX,'System matrix assembled.']);

% �鿴ϡ�����ķ���Ԫ��
% spy(A)

% ============================== ��� ==================================

% ��ʾ��Ϣ
disp([DISP_PREFIX,'Solve the system...']);

% ֱ�����
p = A\RHS;

% ��ʾ��Ϣ
disp([DISP_PREFIX,'System solved.'])

% �鿴���
res_vec = A*p - RHS;
res_norm = norm(res_vec);

disp([DISP_PREFIX,'The norm of residual vector is ',num2str(res_norm),'.']);

% ============================= ��ʾ��� ==============================

% �ع���������γ�һ����������Ϊ��TH_DIM+1���У���RA_DIM+1����
pr = reshape(p,TH_DIM+1,RA_DIM+1);
% th_idx_re = ((1:1:(TH_DIM+1))-1) .* Dt ./ (2*pi) .* 360; % re for reconstruction
% ra_idx_re = ((1:1:(RA_DIM+1))-1) .* Da + PAD_DIM(1,3);     % re for reconstruction
% disp('Display the pressure in 3D.');
% surf(ra_idx_re,th_idx_re,pr,'LineStyle','none')

% ========================== ��ѹ��������� ============================
% ���λ���
fx = 0; % ���غɴ�ֱ�ĺ���
fz = 0; % ���غ�ƽ�еĺ���
ft = 0; % ��ʱ����
% fth_begin = -1*Dt/2; % ��ʱ������force theta
% fth = 0;
fa = PAD_DIM(1,3);

for I = 1:1:RA_DIM
%     fth = fth_begin;
    for J = 1:1:TH_DIM
%         fth = fth + Dt;
        
%         if(J ~= TH_DIM)
        if(1)
            % δ���ػ��߽�
            ft = pr(  J,I)*fa + pr(  J,I+1)*(fa+Da) + ...
                 pr(J+1,I)*fa + pr(J+1,I+1)*(fa+Da);
        else
            % ����ػ��߽�
            ft = pr(J,I) + pr(J,I+1) + ...
                 pr(1,I) + pr(1,I+1);
        end
        
        if(ft > 0)
%             fx = fx +      sin(fth)*ft;
%             fz = fz + (-1)*cos(fth)*ft;
            % ==== 2012.3.29 =======
            % ����ֹ����е�������
            fx = fx + sin(ALPHA)*ft;
            fz = fz + cos(ALPHA)*ft;
        end
%         fx = fx +      sin(fth)*ft;
%         fz = fz + (-1)*cos(fth)*ft;
    end % J
    fa = fa + Da;
end % I


fx = fx * (Dt * Da)/4;
fz = fz * (Dt * Da)/4;

% =================== check coeff_check ==============================

if(coeff_check == 1)
    clear A
    A = A0;
end