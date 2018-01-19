function [p,A,DIA_IN,DIA_OUT,idx_boundary_in,idx_boundary_out,fx,fz] =...
    ReyThrustFunc(PAD_DIM,TH_DIM,RA_DIM,ns,es,Dt,Da,HP,AS,VISCO,VIS_EN,TURB_SWITCH,RHO,PB,ALPHA,coeff_check,DISP_PREFIX,ANG_OFF)


% 显示输出组装进度的限值
MAX_NONDISP_DIM = 100*100;

% 稀疏矩阵内每行的非零元素个数
ROW_NONZERO_NUM = 16;

ELE_NODES_NUM = 4; % 每个单元的节点个数，这里是4

% 处理压力边界时主对角线上的大数
BIG_COEFF = 1e16;

% =========================== 计算区域离散 ================================

% 显示信息
disp([DISP_PREFIX,'Pre-allocate memory and generate the grid ...']);

% 稀疏矩阵的行索引向量,行向量
row_idx = zeros(1,TH_DIM*RA_DIM*ELE_NODES_NUM^2);

% 稀疏矩阵的列索引向量，行向量
col_idx = zeros(1,TH_DIM*RA_DIM*ELE_NODES_NUM^2);

% 稀疏矩阵的非零值向量，行向量
elems = zeros(1,TH_DIM*RA_DIM*ELE_NODES_NUM^2);
if(coeff_check == 1)
    elems_h0 = elems;
else
    elems_h0 = 0;
end

% ===== 2012.3.24 ======
% 减少了第一坐标，周向的节点数
% 方程右端，RHS，列向量
% RHS = zeros((TH_DIM+1)*(RA_DIM+1),1);
% ===== 2012.3.28 ======
% 增加了第一坐标
RHS = zeros((TH_DIM+1)*(RA_DIM+1),1);

% 方成的解
p = 0;

disp([DISP_PREFIX,'Memory allocated and grid generated.']);

% ========================= 组装系数矩阵 =================================

% 显示信息
disp([DISP_PREFIX,'Begin assemble the system...']);

% 临时节点矩阵,按行存储，第一列是节点编号
tn = zeros(4,3);

% 循环临时变量
t_th = -1  / (VISCO*VIS_EN);  % 周向的临时变量
t_ra = -1  / (VISCO*VIS_EN); % 轴向的临时变量
h = 0; % 液膜厚度
Ke = zeros(ELE_NODES_NUM,ELE_NODES_NUM); % 单元的刚度矩阵
OT  =  1/3; % one thirds
MOT = -1/3; % negtive one thirds, 因为NOT可能与关键字重复，故用MOT
OS  =  1/6; % one sixths
NOS = -1/6; % negtive one sixths
idx_offset = 0; % 全局索引向量的当前单元偏移
idx_off_base = ELE_NODES_NUM^2; % 全局索引向量的单元偏移基
tRe = RHO*AS/(VISCO*VIS_EN); % 湍流计算时局部雷诺数计算临时变量

% % ======== 2012.3.26 =========
% % Jacobian 矩阵的对角元素，临时变量
% % 对于目前固定的网格
% Dt = (PAD_DIM(1,2)-PAD_DIM(1,1))/TH_DIM; % 单元的第一坐标长度
% Da = (PAD_DIM(1,4)-PAD_DIM(1,3))/RA_DIM; % 单元的第二坐标长度
Dtn = 1/Dt;
Dan = 1/Da;
det_J = Dt*Da;
% 重新给变量赋值
t_th = t_th*Dtn^2*det_J;
t_ra = t_ra*Dan^2*det_J;

% 临时周向系数矩阵
MatrixTh = [ 
     OT,MOT,NOS, OS;
    MOT, OT, OS,NOS;
    NOS, OS, OT,MOT;
     OS,NOS,MOT, OT
];

% 临时轴向系数矩阵
MatrixRa = [
     OT, OS,NOS,MOT;
     OS, OT,MOT,NOS;
    NOS,MOT, OT, OS;
    MOT,NOS, OS, OT
];

% 方程右端临时向量，列向量
tRHS = -1*0.5 * AS * [-1,1,1,-1]';
rhs_idx = zeros(1,4); % 方程右端索引
% ======== 2012.3.26 =========
% 方程右端考虑Jacobian矩阵，第一坐标方向
tRHS = tRHS .* [Dtn,Dtn,Dtn,Dtn]' * det_J; 

% 止推瓦块的起始
ths = PAD_DIM(1,1);
thend = PAD_DIM(1,2);

% 止推轴承湍流系数计算用的参数
MS  = -0.25;
NS  = 0.066;
tGx = 2^(1+MS) / (NS * (2 + MS));
tGz = 2^(1+MS) / NS;

for I = 1:1:TH_DIM*RA_DIM
    % 组装每一个单元
    
    % 取出四个节点
    for J = 1:1:ELE_NODES_NUM
        tn(J,:) = [es(I,J),ns(es(I,J),1),ns(es(I,J),2)];
    end % J
    
    % 计算单元的中心坐标, theta center和axis center
    thc = (tn(1,2) + tn(2,2))/2 + ANG_OFF;
    rac = (tn(1,3) + tn(4,3))/2;
    
    % 计算液膜厚度
%     h = HP + ALPHA*rac*sin(thc - ths);
    h = HP + ALPHA*rac*sin(thend - thc);
    
    h3 = h^3;
    % 计算圆周方向的系数Gx和Gz
    if(TURB_SWITCH == 0)
        Gx = t_th * h3 / 12 / rac;
        Gz = t_ra * h3 / 12 * rac;
    elseif(TURB_SWITCH == 1)
        lRe = tRe * h * rac;
        gRe = lRe^(1+MS);
        if(lRe > 2000)
            % ====== 2012.3.31 =======
            % 修正为止推轴承的湍流模型
%             Gx = t_th * h3 / (12+0.0136*lRe^0.9) / rac;
%             Gz = t_ra * h3 / (12+0.0043*lRe^0.96) * rac;
            Gx = t_th * h3 * tGx / gRe / rac;
            Gz = t_ra * h3 * tGz / gRe * rac;
        else
            Gx = t_th * h3 / 12 / rac;
            Gz = t_ra * h3 / 12 * rac;
        end % lRe > 2000
    end   

    % 计算这个单元的广义刚度矩阵
    Ke = Gx * MatrixTh + Gz * MatrixRa;
    % 将广义刚度矩阵变成向量
    Ke_row = reshape(Ke,1,idx_off_base);
    
    % 将值插入刚度矩阵向量
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
        % 计算圆周方向的系数Gx和Gz
        if(TURB_SWITCH == 0)
            Gx = t_th * h3 / 12 / rac;
            Gz = t_ra * h3 / 12 * rac;
        elseif(TURB_SWITCH == 1)
            lRe = tRe * h * rac;
            if(lRe > 2000)
                % ====== 2012.3.31 =======
                % 修正为止推轴承的湍流模型
%                 Gx = t_th * h3 / (12+0.0136*lRe^0.9) / rac;
%                 Gz = t_ra * h3 / (12+0.0043*lRe^0.96) * rac;
                Gx = tGx / gRe / rac;
                Gz = tGz / gRe * rac;
            else
                Gx = t_th * h3 / 12 / rac;
                Gz = t_ra * h3 / 12 * rac;
            end % lRe > 2000
        end   

        % 计算这个单元的广义刚度矩阵
        Ke = Gx * MatrixTh + Gz * MatrixRa;
        % 将广义刚度矩阵变成向量
        Ke_row = reshape(Ke,1,idx_off_base);

        elems_h0(1,(idx_offset+1):1:(idx_offset+idx_off_base)) = Ke_row;
    end

    idx_offset = idx_offset + idx_off_base;
    
    % 计算方程右端
    rhs_idx = [tn(1,1),tn(2,1),tn(3,1),tn(4,1)];
    RHS(rhs_idx,1) = RHS(rhs_idx,1) + h .* tRHS .* rac;
end

% 生成稀疏系数矩阵
A = sparse(row_idx,col_idx,elems,...
    (TH_DIM+1)*(RA_DIM+1),(TH_DIM+1)*(RA_DIM+1)...
    );
if(coeff_check == 1)
    A0 = sparse(row_idx,col_idx,elems_h0,...
        (TH_DIM+1)*(RA_DIM+1),(TH_DIM+1)*(RA_DIM+1)...
        );
end

% 显示信息
disp([DISP_PREFIX,'Handle the boundaries...']);

% ===== 2012.3.24 ======
% 缩小周向的节点数
% % 处理边界条件，进出口压力边界
% % 轴向坐标为零和最大值时，待求解的压力值是已知的
% % 找到压力边界点的索引
% idx_boundary_in  = 1:1:(TH_DIM+1);
% idx_boundary_out = ((TH_DIM+1)*RA_DIM+1):1:((TH_DIM+1)*(RA_DIM+1));
% % 获得系数矩阵的对角线
% DIA_A = diag(A);
%
% 处理边界条件，进出口压力边界
% 轴向坐标为零和最大值时，待求解的压力值是已知的
% 找到压力边界点的索引
idx_boundary_in  = [1:1:(TH_DIM+1),1:(TH_DIM+1):((TH_DIM+1)*RA_DIM+1)];
idx_boundary_out = [((TH_DIM+1)*RA_DIM+1):1:((TH_DIM+1)*(RA_DIM+1)),...
    ( (TH_DIM+1):(TH_DIM+1):((TH_DIM+1)*(RA_DIM+1)) )];
% 获得系数矩阵的对角线
DIA_A = diag(A);

if(coeff_check == 1)
    DIA_A0 = diag(A0);
end

% 获得调整后的入口边界和出口边界的对角线系数
DIA_IN  = DIA_A( idx_boundary_in) + BIG_COEFF;
DIA_OUT = DIA_A(idx_boundary_out) + BIG_COEFF;

% 处理方程右端，入口边界
RHS( idx_boundary_in) = PB*DIA_IN;
RHS(idx_boundary_out) = PB*DIA_OUT;

if(coeff_check == 1)
    DIA_IN  = DIA_A0( idx_boundary_in) + BIG_COEFF;
    DIA_OUT = DIA_A0(idx_boundary_out) + BIG_COEFF;  
end

% 释放稀疏矩阵的内存
clear A
clear A0

% 接长行索引和列索引
final_row_idx = [row_idx,idx_boundary_in,idx_boundary_out];
final_col_idx = [col_idx,idx_boundary_in,idx_boundary_out];
final_elems   = [elems,ones(1,TH_DIM+1+RA_DIM+1)*BIG_COEFF,ones(1,TH_DIM+1+RA_DIM+1)*BIG_COEFF];

if(coeff_check == 1)
    final_elems_h0 = [elems_h0,ones(1,TH_DIM+1+RA_DIM+1)*BIG_COEFF,ones(1,TH_DIM+1+RA_DIM+1)*BIG_COEFF];
else
    final_elems_h0 = 0;
end

% 生成最终稀疏系数矩阵
A = sparse(final_row_idx,final_col_idx,final_elems,...
    (TH_DIM+1)*(RA_DIM+1),...
    (TH_DIM+1)*(RA_DIM+1));

if(coeff_check == 1)
A0 = sparse(final_row_idx,final_col_idx,final_elems_h0,...
    (TH_DIM+1)*(RA_DIM+1),...
    (TH_DIM+1)*(RA_DIM+1));
end

% 释放部分临时内存
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

% 显示信息
disp([DISP_PREFIX,'System matrix assembled.']);

% 查看稀疏矩阵的非零元素
% spy(A)

% ============================== 求解 ==================================

% 显示信息
disp([DISP_PREFIX,'Solve the system...']);

% 直接求解
p = A\RHS;

% 显示信息
disp([DISP_PREFIX,'System solved.'])

% 查看结果
res_vec = A*p - RHS;
res_norm = norm(res_vec);

disp([DISP_PREFIX,'The norm of residual vector is ',num2str(res_norm),'.']);

% ============================= 显示结果 ==============================

% 重构结果向量形成一个矩阵矩阵的为（TH_DIM+1）行，（RA_DIM+1）列
pr = reshape(p,TH_DIM+1,RA_DIM+1);
% th_idx_re = ((1:1:(TH_DIM+1))-1) .* Dt ./ (2*pi) .* 360; % re for reconstruction
% ra_idx_re = ((1:1:(RA_DIM+1))-1) .* Da + PAD_DIM(1,3);     % re for reconstruction
% disp('Display the pressure in 3D.');
% surf(ra_idx_re,th_idx_re,pr,'LineStyle','none')

% ========================== 对压力场求积分 ============================
% 梯形积分
fx = 0; % 与载荷垂直的合力
fz = 0; % 与载荷平行的合力
ft = 0; % 临时变量
% fth_begin = -1*Dt/2; % 临时变量，force theta
% fth = 0;
fa = PAD_DIM(1,3);

for I = 1:1:RA_DIM
%     fth = fth_begin;
    for J = 1:1:TH_DIM
%         fth = fth + Dt;
        
%         if(J ~= TH_DIM)
        if(1)
            % 未到回环边界
            ft = pr(  J,I)*fa + pr(  J,I+1)*(fa+Da) + ...
                 pr(J+1,I)*fa + pr(J+1,I+1)*(fa+Da);
        else
            % 到达回环边界
            ft = pr(J,I) + pr(J,I+1) + ...
                 pr(1,I) + pr(1,I+1);
        end
        
        if(ft > 0)
%             fx = fx +      sin(fth)*ft;
%             fz = fz + (-1)*cos(fth)*ft;
            % ==== 2012.3.29 =======
            % 修正止推轴承的力方向
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