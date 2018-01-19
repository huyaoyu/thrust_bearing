%
% Description
% ===========
%
% 本文件用于测试计算止推轴承的轴向动特性系数。计算的数学模型来自于Gunhee Jang
% 的文献“Determination of the dynamic coefficients of the coupled journal
% and thrust bearings by the perturbation method”。计算中没有考虑传热。计算
% 中考虑了止推轴承的湍流作用，湍流修正采用的是Eskild Storteig等人的文献，“
% Dynamic characteristics of hydrodynamically lubricated fixed-pad thrust
% bearings”.
%
% 求解方程时采用了二维有限元法，使用四节点等参数单元。本计算程序包括稳态平衡位置
% 的求解，将借助与求解稳太Reynolds方程的子程序ReyThrustFunc，改程序负责求解一次
% 给定轴向位置的止推轴承单瓦Reynolds方程。
%
% 求解静平衡位置利用Newton-Raphson迭代，迭代格式参考复旦大学王丽萍博士学位论文
% 《可倾瓦轴承动力学建模及动力特性研究》一文中关于径向轴承平衡位置的求解过程。
%
% - 3月31日，修正了由于积分坐标顺序搞错而带来的数学模型推导错误。增加了阻尼求解。
% 将文献中到的方向参数a修正为a和atd两个参数，分别用于计算Kzz和Czz。
% 
% 注意，在积分时不考虑边界的回环。
%
% - 4月1日，在求解NR迭代完毕后显示最终的压力分布。增加计算时间测试。
% - 4月9日，增加了转角刚度和阻尼计算。
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
% 创建：2012年3月30日
% 修改：2012年3月31日
% 修改：2012年4月1日
% 修改：2012年4月9日




% ========================= 清理工作空间 ==================================

close ALL
clear
clc

disp('Reynolds Equation Solver, for thrust bearing coefficients, with FEM');

ThrustInput;

% Newton-Raphson 迭代

qk  = HP;
qk1 = 0; % 上一个

% 迭代信息前缀
iter_count = 1;
iter_prefix = '';

% 计算区域的离散化参数，单元个数
TH_DIM = 100;     % 圆周方向的离散，注意这里是指的单元数
RA_DIM = 100;      % 轴向的离散，注意这里是指的单元数

% ===== 2012.3.24 ======
% 使用处理回环的网格生成函数
% 网格信息
% [ns,es] = RecField2DIso([0,0],[2*pi,0.1],[TH_DIM,RA_DIM]);
% [ns,es] = RecField2DIsoWrap([0,0],[2*pi,0.1],[TH_DIM,RA_DIM],[1,0]);
% ===== 2012.3.31 =======
% 注意起始位置
[ns,es] = RecField2DIsoWrap(...
    [PAD_DIM(1,1),PAD_DIM(1,3)],...
    [PAD_DIM(1,2)-PAD_DIM(1,1),PAD_DIM(1,4)-PAD_DIM(1,3)],...
    [TH_DIM,RA_DIM],[0,0]);

% ======== 2012.3.26 =========
% Jacobian 矩阵的对角元素，临时变量
% 对于目前固定的网格
Dt = (PAD_DIM(1,2)-PAD_DIM(1,1))/TH_DIM; % 单元的第一坐标长度
Dr = (PAD_DIM(1,4)-PAD_DIM(1,3))/RA_DIM; % 单元的第二坐标长度

while(abs(qk-qk1) > NEWTON_RAPHSON_NORM)
    % 求解一次轴向刚度
    
    % 信息前缀
    iter_prefix = ['NR ',num2str(iter_count),': '];
    
    [p,pz,A,DIA_IN,DIA_OUT,idx_boundary_in,idx_boundary_out,fx,fz,Kzz] = ReyThrustStiffFunc(...
    PAD_DIM,TH_DIM,RA_DIM,ns,es,Dt,Dr,qk,AS,VISCO,VIS_EN,TURB_SWITCH,RHO,PB,ALPHA,iter_prefix,0);

    qa  = -1 * (fz - W) / (-1 * Kzz);
    qk1 = qk;
    qk  = qk + qa;
    
    iter_count = iter_count + 1;
end % abs(qk-qk1) < NEWTON_RAPHSON_NORM

% ========================= 显示压力场求解结果 ==========================

% 重构结果向量形成一个矩阵矩阵的为（TH_DIM+1）行，（RA_DIM+1）列
pr = reshape(p,TH_DIM+1,RA_DIM+1); % 压力场重构
th_idx_re = ((1:1:(TH_DIM+1))-1) .* Dt ./ (2*pi) .* 360; % re for reconstruction
ra_idx_re = ((1:1:(RA_DIM+1))-1) .* Dr + PAD_DIM(1,3);     % re for reconstruction

disp('Display the final pressure field in 3D.');

surf(ra_idx_re,th_idx_re,pr,'LineStyle','none')

% ========================= 组装方程右端 ==============================

disp('Solve the damping coeffiecents ...');

% 转换有限元单元积分参数
Dtn = 1/Dt;
Drn = 1/Dr;
det_J = Dt*Dr;

% 用于表示方向的变量
% a    = 1;
% atx  = -1; % x摆角的alpha值
% aty  = 1; % y摆角的alpha值
a    = 1;
atx  = -1; % x摆角的alpha值
aty  = 1; % y摆角的alpha值
atd  = 1;
atxd = -1;
atyd = 1;
ang_off = 0; % 横截面坐标系的角度偏移

tRe = RHO*AS/(VISCO*VIS_EN); % 湍流计算时局部雷诺数计算临时变量

% 方程右端临时向量，列向量
rhs_idx = zeros(1,4); % 方程右端索引

% 阻尼计算的临时方程右端
tRHSzt  =  atd * det_J * [1,1,1,1]';
tRHStxt = atxd * det_J * [1,1,1,1]';
tRHStyt = atyd * det_J * [1,1,1,1]';

% RHS
RHSzt = zeros((TH_DIM+1)*(RA_DIM+1),1);
% ======= 2012.4.9（开始） ==========
th_x = 0; % x 方向的转角
th_y = 0; % y 方向的转角
cz0 = 0; % center of mass 的 z 坐标
rz  = 0; % 摆动中心的z坐标
tphi1 = (cz0 - rz)*sin(th_x);
tphi2 = -1*cos(th_x);
tpsi1 = (cz0 - rz)*sin(th_y);
tpsi2 = cos(th_y);
RHStx = RHSzt; % x 方向的转角，theta x
RHSty = RHSzt; % y 方向的转角，theta y
% tRHSas = -1*0.5 * AS * [-1,1,1,-1]' .*...
%     [Dtn,Dtn,Dtn,Dtn]' * det_J; % 'as' for angluar speed
tRHSas = 0.5 * AS * [-1,1,1,-1]' .*...
    [Dtn,Dtn,Dtn,Dtn]' * det_J; % 'as' for angluar speed
RHStxt = RHSzt;
RHStyt = RHSzt;

% 周向（第一坐标）和径向（第二坐标）的方程RHS临时变量
% tRHSAr = (VISCO*VIS_EN) * Drn^2 * det_J * [1/3,-1/3,1/3,-1/3]';
% tRHSAt = (VISCO*VIS_EN) * Dtn^2 * det_J * [1/3,-1/3,1/3,-1/3]';
% tRHSB  = (VISCO*VIS_EN) * Dtn^2 * det_J * [-1,1,1,-1]'; % t方向
% tRHSC  = (VISCO*VIS_EN) * Drn^2 * det_J * [-1,-1,1,1]'; % r方向
tRHSAr = -1*(VISCO*VIS_EN) * Drn^2 * det_J * [1/3,-1/3,1/3,-1/3]';
tRHSAt = -1*(VISCO*VIS_EN) * Dtn^2 * det_J * [1/3,-1/3,1/3,-1/3]';
tRHSB  = -1*(VISCO*VIS_EN) * Dtn^2 * det_J * [-1,1,1,-1]'; % t方向
tRHSC  = -1*(VISCO*VIS_EN) * Drn^2 * det_J * [-1,-1,1,1]'; % r方向


% 止推轴承湍流系数计算用的参数
MS  = -0.25;
NS  = 0.066;
tGx = 2^(1+MS) / (NS * (2 + MS));
tGz = 2^(1+MS) / NS;
if(IS_DEBUG)
    check_reynolds = zeros(1,TH_DIM*RA_DIM); % 记录局部雷诺数用于调试
    check_rac      = zeros(1,TH_DIM*RA_DIM); % 记录单元中心位置用于调试
    check_h        = zeros(1,TH_DIM*RA_DIM); % 记录液膜厚度用于调试
end

% ======= 2012.4.9（结束） ==========


% 止推瓦块的起始
ths = PAD_DIM(1,1);
thend = PAD_DIM(1,2);

% 显示信息
disp('Assemble RHSs of perturbation equations...');

for I = 1:1:TH_DIM*RA_DIM
    % 组装每一个单元
    
    % 取出四个节点
    for J = 1:1:ELE_NODES_NUM
        tn(J,:) = [es(I,J),ns(es(I,J),1),ns(es(I,J),2)];
    end % J
    
    % 计算单元的中心坐标, theta center和axis center
    thc = (tn(1,2) + tn(2,2))/2 + ang_off;
    rac = (tn(1,3) + tn(4,3))/2;
    
    % 计算液膜厚度
    h = qk1 + ALPHA*rac*sin(thend - thc);
    
    % 计算phi和psi
    phi_t = tphi1 + rac * tphi2 * sin(thc);
    psi_t = tpsi1 + rac * tpsi2 * cos(thc);
    
    % ========== 2012.4.9（开始） ========
    % 止推轴承的描述
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
    
    % 临时节点压力
    p0t = p(tn(:,1),1);
    
    p0A = sum(0.25 * [1,-1,1,-1]' .* p0t);
    p0B = sum(0.25 * [-1,1,1,-1]' .* p0t);
    p0C = sum(0.25 * [-1,-1,1,1]' .* p0t);
    
    ttRHS = (tRHSAt .* p0A .* h^2 + tRHSB .* p0B .* h^2) * Gx;
    trRHS = (tRHSAr .* p0A .* h^2 + tRHSC .* p0C .* h^2) * Gz;

  
    % 计算方程右端
    rhs_idx = [tn(1,1),tn(2,1),tn(3,1),tn(4,1)];
    RHSReadyToAdd  = trRHS + ttRHS + tRHSas .* rac;
    RHSReadyToAddX = atx * RHSReadyToAdd;
    RHSReadyToAddY = aty * RHSReadyToAdd;
    RHStx(rhs_idx,1) = RHStx(rhs_idx,1) + RHSReadyToAddX * phi_t;
    RHSty(rhs_idx,1) = RHSty(rhs_idx,1) + RHSReadyToAddY * psi_t;
    RHStxt(rhs_idx,1) = RHStxt(rhs_idx,1) + tRHStxt * rac * phi_t;
    RHStyt(rhs_idx,1) = RHStyt(rhs_idx,1) + tRHStyt * rac * psi_t;
    % ========== 2012.4.9（结束） ========
      
    % 计算方程右端
    RHSzt(rhs_idx,1) = RHSzt(rhs_idx,1) + tRHSzt * rac;
end

% 处理方程右端，特别注意边界！！！还未修改
RHSzt( idx_boundary_in) = PB*DIA_IN;
RHSzt(idx_boundary_out) = PB*DIA_OUT;

% ========== 2012.4.9（开始） ========
RHStx( idx_boundary_in) = PB*DIA_IN;
RHStx(idx_boundary_out) = PB*DIA_OUT;

RHSty( idx_boundary_in) = PB*DIA_IN;
RHSty(idx_boundary_out) = PB*DIA_OUT;

RHStxt( idx_boundary_in) = PB*DIA_IN;
RHStxt(idx_boundary_out) = PB*DIA_OUT;

RHStyt( idx_boundary_in) = PB*DIA_IN;
RHStyt(idx_boundary_out) = PB*DIA_OUT;
% ========== 2012.4.9（结束） ========

% 显示信息
disp('RHSs assembled.');

% =========================== Solve ==================================

% z方向的速度摄动压力场
pzt = A\RHSzt;
disp('pzt solved.');

res_vec = A*pzt - RHSzt;
res_norm = norm(res_vec);

disp(['The norm of residual vector of pzt is ',num2str(res_norm),'.']);

% theta x 方向的摄动压力场
ptx = A\RHStx;
disp('ptx solved.');

res_vec = A*ptx - RHStx;
res_norm = norm(res_vec);

disp(['The norm of residual vector of ptx is ',num2str(res_norm),'.']);
% theta y 方向的摄动压力场
pty = A\RHSty;
disp('pty solved.');

res_vec = A*pty - RHSty;
res_norm = norm(res_vec);

disp(['The norm of residual vector of pty is ',num2str(res_norm),'.']);

% theta x 方向的速度摄动压力场
ptxt = A\RHStxt;
disp('pty solved.');

res_vec = A*ptxt - RHStxt;
res_norm = norm(res_vec);

disp(['The norm of residual vector of ptxt is ',num2str(res_norm),'.']);

% theta y 方向的速度摄动压力场
ptyt = A\RHStyt;
disp('pty solved.');

res_vec = A*ptyt - RHStyt;
res_norm = norm(res_vec);

disp(['The norm of residual vector of ptyt is ',num2str(res_norm),'.']);

% ============================= 显示结果 ==============================

% 重构结果向量形成一个矩阵矩阵的为（TH_DIM+1）行，（RA_DIM+1）列
% z 方向阻尼
pztr = reshape(pzt,TH_DIM+1,RA_DIM+1); % 速度摄动压力场重构
pzr = reshape(pz,TH_DIM+1,RA_DIM+1); % 速度摄动压力场重构

disp('Display the results in 3D.');

figure

surf(ra_idx_re,th_idx_re,pztr,'LineStyle','none')

view(0,90);
colorbar;
title('pressure field of velocity perturbation along z');

% theta x 方向刚度
ptxr = reshape(ptx,TH_DIM+1,RA_DIM+1); % 速度摄动压力场重构

disp('Display the displacement perturbation of theta x in 3D.');

figure

surf(ra_idx_re,th_idx_re,ptxr,'LineStyle','none')
view(0,90);
colorbar;
title('pressure field of a-displacement perturbation along theta x');

% theta y 方向刚度
ptyr = reshape(pty,TH_DIM+1,RA_DIM+1); % 速度摄动压力场重构

disp('Display the displacement perturbation of theta y in 3D.');

figure

surf(ra_idx_re,th_idx_re,ptyr,'LineStyle','none')
view(0,90);
colorbar;
title('pressure field of a-displacement perturbation along theta y');

% theta x 方向刚度
ptxtr = reshape(ptxt,TH_DIM+1,RA_DIM+1); % 速度摄动压力场重构

disp('Display the velocity perturbation of theta x in 3D.');

figure

surf(ra_idx_re,th_idx_re,ptxtr,'LineStyle','none')
view(0,90);
colorbar;
title('pressure field of a-velocity perturbation along theta x');

% theta y 方向刚度
ptytr = reshape(ptyt,TH_DIM+1,RA_DIM+1); % 速度摄动压力场重构

disp('Display the velocity perturbation of theta y in 3D.');

figure

surf(ra_idx_re,th_idx_re,ptytr,'LineStyle','none')
view(0,90);
colorbar;
title('pressure field of a-velocity perturbation along theta y');

% ========================= 对摄动压力场求积分 ===========================
% % 梯形积分
% Czz = 0; % 与载荷平行的阻尼
% Czz_temp = 0; % 临时变量
% fth_begin = -1*Dt/2; % 临时变量，force theta
% fth = 0; % 临时变量，角度位置，单位为rad
% 
% for I = 1:1:RA_DIM
%     fth = fth_begin;
%     for J = 1:1:TH_DIM
%         fth = fth + Dt;
%         
% %         if(J ~= TH_DIM)
%         if(1)
%             % 未到回环边界
%             Czz_temp = pztr(  J,I) + pztr(  J,I+1) + ...
%                        pztr(J+1,I) + pztr(J+1,I+1);
%         else
%             % 到达回环边界
%             Czz_temp = pztr(J,I) + pztr(J,I+1) + ...
%                        pztr(1,I) + pztr(1,I+1);
%         end
%         
%         Czz = Czz + (-1)*Czz_temp;
% 
%     end % J
% end % I
% 
% % 乘以积分系数
% Czz = Czz * (Dt * Dr)/4;

alpha = [0,0,a,atx,aty]';
K = DouIntegration(alpha,pzr,ptxr,ptyr,Dt,Dr,PAD_DIM(1,3));
C = DouIntegration(alpha,pztr,ptxtr,ptytr,Dt,Dr,PAD_DIM(1,3));

% 显示刚度信息
disp(['Czz = ',num2str(C(3,3)/DAMPING_BASE),'x10^6 Ns/m']);

if(IS_DEBUG)
    disp('Debug mode is on.');
end
