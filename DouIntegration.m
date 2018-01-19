function M = DouIntegration(alpha,pzr,ptxr,ptyr,Dt,Dr,fa0)
% 专门对止推轴承的摄动压力场进行二重积分的函数。
% alpha - 计算时使用的alpha值向量，列向量，5行1列
% pzr - z方向的计算结果，重构过的矩阵
% ptxr - theta x方向的计算结果，重构过的矩阵
% ptyr - theta y方向的计算结果，重构过的矩阵
% Dt - 角度增量
% Dr - 半径增量
% fa0 - 半径的开始值
% M - 积分得到的矩阵

% 本函数的积分格式采用梯形积分，并使用了Myunggyu Kim等人的文献（2010）"
% Stability analysis of a disk-spindle system supported bu coupled journal
% and thrust bearings considering five degrees of freedom"

ROW = 5;
COL = 5;

% ==================== 维度检查 =========================

[row,col] = size(alpha);
if(row ~= ROW || col ~= 1)
    disp('The dimension of alpha is wrong!');
    M = 0;
    return
end

[row1,col1] = size(pzr);
[row2,col2] = size(ptxr);
[row3,col3] = size(ptyr);

if(row1 ~= row2 || row1 ~= row3 || row2 ~= row3 || ...
        col1 ~= col2 || col1 ~= col3 || col2 ~= col3)
    disp('The dimensions of pressures are wrong!');
    M = 0;
    return
end

% ======================= 变量声明 ========================

TH_DIM = row1 - 1;
RA_DIM = col1 - 1;

M = zeros(ROW,COL);

av = alpha .* [1,1,-1,1,-1]';
pv = zeros(1,COL);

ta = av;

fth_begin = -1*Dt/2;
fr_begin  = -1*Dr/2;
fth = 0;
fr  = 0;
fa  = fa0;

for I = 1:1:RA_DIM
    fth = fth_begin;
    fr  = fr + Dr;
    for J = 1:1:TH_DIM
        fth = fth + Dt;
        
        % 准备alpha向量
        ta = av .* [1,1,1,fr*sin(fth),fr*cos(fth)]';
        
        % 准备p向量
%         pv(1,3) = pzr(  J,I)*fa + pzr(  J,I+1)*(fa+Dr) + ...
%                   pzr(J+1,I)*fa + pzr(J+1,I+1)*(fa+Dr);
%         pv(1,4) = ptxr(  J,I)*fa + ptxr(  J,I+1)*(fa+Dr) + ...
%                   ptxr(J+1,I)*fa + ptxr(J+1,I+1)*(fa+Dr);
%         pv(1,5) = ptyr(  J,I)*fa + ptyr(  J,I+1)*(fa+Dr) + ...
%                   ptyr(J+1,I)*fa + ptyr(J+1,I+1)*(fa+Dr);
        pv(1,3) = pzr(  J,I)*fa*fr + pzr(  J,I+1)*(fa+Dr)*(fr+Dr) + ...
                  pzr(J+1,I)*fa*fr + pzr(J+1,I+1)*(fa+Dr)*(fr+Dr);
        pv(1,4) = ptxr(  J,I)*fa*fr + ptxr(  J,I+1)*(fa+Dr)*(fr+Dr) + ...
                  ptxr(J+1,I)*fa*fr + ptxr(J+1,I+1)*(fa+Dr)*(fr+Dr);
        pv(1,5) = ptyr(  J,I)*fa*fr + ptyr(  J,I+1)*(fa+Dr)*(fr+Dr) + ...
                  ptyr(J+1,I)*fa*fr + ptyr(J+1,I+1)*(fa+Dr)*(fr+Dr);
              
        % 向量相乘
        M = M + ta * pv;
        
    end % J
    fa = fa + Dr;
end % I

M = M * (Dt*Dr)/4;

