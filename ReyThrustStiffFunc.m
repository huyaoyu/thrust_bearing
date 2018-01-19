function [p,pz,A,DIA_IN,DIA_OUT,idx_boundary_in,idx_boundary_out,fx,fz,Kzz] = ReyThrustStiffFunc(...
PAD_DIM,TH_DIM,RA_DIM,ns,es,Dt,Dr,HP,AS,VISCO,VIS_EN,TURB_SWITCH,RHO,PB,ALPHA,DISP_PREFIX,ANG_OFF)

% 调试开关
IS_DEBUG = 1;

% 计算常量
ELE_NODES_NUM = 4;        % 每个单元的节点个数，四节点

STIFFNESS_BASE = 1e9;     % 用于将求得的刚度值以10^9形式显示
% DAMPING_BASE   = 1e6;     % 用于将求得的阻尼值以10^6形式显示

% 求解一次稳态压力场
STABLE_FIELD_DISP_PREFIX = [DISP_PREFIX,'Stable pressure field: '];
[p,A,DIA_IN,DIA_OUT,idx_boundary_in,idx_boundary_out,fx,fz] = ReyThrustFunc(...
PAD_DIM,TH_DIM,RA_DIM,ns,es,Dt,Dr,HP,AS,VISCO,VIS_EN,TURB_SWITCH,RHO,PB,ALPHA,0,STABLE_FIELD_DISP_PREFIX,ANG_OFF);

% ========================= 组装方程右端 ==============================

% 转换有限元单元积分参数
Dtn = 1/Dt;
Drn = 1/Dr;
det_J = Dt*Dr;

% 用于表示方向的变量
a   = -1;
atd = 1;

tRe = RHO*AS/(VISCO*VIS_EN); % 湍流计算时局部雷诺数计算临时变量

% 方程右端临时向量，列向量
rhs_idx = zeros(1,4); % 方程右端索引

% 周向（第一坐标）和径向（第二坐标）的方程RHS临时变量
tRHSAr = -1*a/(VISCO*VIS_EN) * Drn^2 * det_J * [1/3,-1/3,1/3,-1/3]';
tRHSAt = -1*a/(VISCO*VIS_EN) * Dtn^2 * det_J * [1/3,-1/3,1/3,-1/3]';
tRHSB  = -1*a/(VISCO*VIS_EN) * Dtn^2 * det_J * [-1,1,1,-1]'; % t方向
tRHSC  = -1*a/(VISCO*VIS_EN) * Drn^2 * det_J * [-1,-1,1,1]'; % r方向

% 阻尼计算的临时方程右端
% tRHSzt = atd * det_J * [1,1,1,1]';

% RHS
RHSz  = zeros((TH_DIM+1)*(RA_DIM+1),1);
% RHSzt = RHSz;

% 止推瓦块的起始
ths = PAD_DIM(1,1);
thend = PAD_DIM(1,2);

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

% 显示信息
disp([DISP_PREFIX,'Assemble RHSs of perturbation equations...']);

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
    h = HP + ALPHA*rac*sin(thend - thc);
    
    % 计算圆周方向的系数Gx和Gz
    % 径向轴承的描述
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
    RHSz(rhs_idx,1) = RHSz(rhs_idx,1) + trRHS + ttRHS;
%     RHSzt(rhs_idx,1) = RHSzt(rhs_idx,1) + tRHSzt * rac;
end

% 处理方程右端，特别注意边界！！！还未修改
RHSz( idx_boundary_in) = PB*DIA_IN;
RHSz(idx_boundary_out) = PB*DIA_OUT;

% RHSzt( idx_boundary_in) = PB*DIA_IN;
% RHSzt(idx_boundary_out) = PB*DIA_OUT;

% 显示信息
disp([DISP_PREFIX,'RHSs assembled.']);

% =========================== Solve ==================================
% % z方向的速度摄动压力场
% pzt = A\RHSzt;
% disp('pzt solved.');

% z方向的位移摄动压力场
pz  = A\RHSz;
disp([DISP_PREFIX,'pz solved.']);

% 查看结果
res_vec = A*pz - RHSz;
res_norm = norm(res_vec);

disp([DISP_PREFIX,'The norm of residual vector of pz is ',num2str(res_norm),'.']);

% res_vec = A*pzt - RHSzt;
% res_norm = norm(res_vec);
% 
% disp(['The norm of residual vector of pzt is ',num2str(res_norm),'.']);

if(IS_DEBUG)
    % 检查局部雷诺数的值
    disp([DISP_PREFIX,'The maximun Reynolds number is ',...
        num2str(max(check_reynolds)),'.']);
    
    disp([DISP_PREFIX,'The maximun r number is ',...
        num2str(max(check_rac)),'.']);
    
    disp([DISP_PREFIX,'The maximun h number is ',...
        num2str(max(check_h)),'.']);

    % 清理临时变量
    clear check_reynolds
    clear check_rac
    clear check_h
end % IS_DEBUG

% ============================= 显示结果 ==============================

% 重构结果向量形成一个矩阵矩阵的为（TH_DIM+1）行，（RA_DIM+1）列
pzr  = reshape(pz,TH_DIM+1,RA_DIM+1);  % 位移摄动压力场重构
% pztr = reshape(pzt,TH_DIM+1,RA_DIM+1); % 速度摄动压力场重构
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

% ========================= 对摄动压力场求积分 ===========================
% 梯形积分
Kzz = 0; % 与载荷平行的刚度
% Czz = 0; % 与载荷平行的阻尼
Kzz_temp = 0; % 临时变量
% Czz_temp = 0; % 临时变量
% fth_begin = -1*Dt/2; % 临时变量，force theta
% fth = 0; % 临时变量，角度位置，单位为rad
fa = PAD_DIM(1,3);

for I = 1:1:RA_DIM
%     fth = fth_begin;
    for J = 1:1:TH_DIM
%         fth = fth + Dt;
        
%         if(J ~= TH_DIM)
        if(1)
            % 未到回环边界
            Kzz_temp = pzr(  J,I)*fa + pzr(  J,I+1)*(fa+Dr) + ...
                       pzr(J+1,I)*fa + pzr(J+1,I+1)*(fa+Dr);
%             Czz_temp = pztr(  J,I) + pztr(  J,I+1) + ...
%                        pztr(J+1,I) + pztr(J+1,I+1);
        else
            % 到达回环边界
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

% 乘以积分系数
Kzz = Kzz * (Dt * Dr)/4;
% Czz = Czz * (Dt * Dr)/4;

% 显示刚度信息
disp([DISP_PREFIX,'Kzz = ',num2str(Kzz/STIFFNESS_BASE),'x10^9 N/m']);
% disp(['Czz = ',num2str(Czz/DAMPING_BASE),'x10^6 Ns/m']);

if(IS_DEBUG)
    disp([DISP_PREFIX,'Debug mode is on.']);
end