function [ Nodes,Elems ] = RecField2DIsoWrap(Origin,Span,Count,WrapCheck)
% 二维平面矩形网格，使用双线性四节点等参数网格。Wrap的含义是将某些坐标的第一
% 和最后一个坐标相等，形成环状的节点坐标分布
% Origin，一个二维行向量，记录原点位置
% Span，一个二维行向量，记录平面上两个坐标方向的长度
% Count，一个二维向量，记录每个坐标方向上的单元个数
% WrapCheck，一个二维向量，几轮需要回环的坐标方向
% Nodes，一个矩阵，每行记录一个点的两个坐标，世界坐标
% Elems， 一个矩阵，每行记录一个单元的四个节点的索引，后手法则
% 矩阵Nodes是以第一坐标为优先的（第一坐标主导），也就是第一坐标最先更改
%
% 2012.3.24 创建
%
% 2012.3.31 修正
% 修正了Origin的使用

% ========================= 函数级常量 ================================

FUN_DISP_PREFIX = 'FUNC RecField2DIso: ';
NODES_PER_ELEM  = 4; % 每个单元的节点个数
COORD_PER_NODE  = 2; % 每个节点的独立坐标个数

% ========================== 初始检查 =================================

% 检查 Origin
[row,col] = size(Origin);
if(row ~= 1 || col ~= COORD_PER_NODE)
    Nodes = 0;
    Elems = 0;
    disp([FUN_DISP_PREFIX,'The dimension of Origin is wrong!']);
    return
end

% 检查 Span
[row,col] = size(Span);
if(row ~= 1 || col ~= COORD_PER_NODE)
    Nodes = 0;
    Elems = 0;
    disp([FUN_DISP_PREFIX,'The dimension of Span is wrong!']);
    return
end

% 检查 Count
[row,col] = size(Count);
if(row ~= 1 || col ~= COORD_PER_NODE)
    Nodes = 0;
    Elems = 0;
    disp([FUN_DISP_PREFIX,'The dimension of Count is wrong!']);
    return
end

% 检查 WrapCheck
[row,col] = size(WrapCheck);
if(row ~= 1 || col ~= COORD_PER_NODE)
    Nodes = 0;
    Elems = 0;
    disp([FUN_DISP_PREFIX,'The dimension of WrapCheck is wrong!']);
    return
end

if(sum(WrapCheck == 0)+sum(WrapCheck == 1) ~= COORD_PER_NODE)
    Nodes = 0;
    Elems = 0;
    disp([FUN_DISP_PREFIX,'The elements of WrapCheck is wrong!']);
    return
end

% ============================== 变量声明 ===============================

% 查看回环,不回环的坐标方向
c = (WrapCheck == 0);

NUM_NODES = (Count(1,1)+c(1,1))*(Count(1,2)+c(1,2)); % 总点数
NUM_ELEMS = Count(1,1)*Count(1,2);         % 总单元数

% 预分配返回值内存
Nodes = zeros(NUM_NODES,COORD_PER_NODE);
Elems = zeros(NUM_ELEMS,NODES_PER_ELEM);

% =========================== 生成节点信息 ===============================

d = Span ./ Count; % 坐标增量，行向量
pos = Origin; % 内部循环的当前坐标位置


for J = 1:1:(Count(1,2) + c(1,2))
    % 第二坐标方向
    
    % 清理当前第一坐标位置
    pos(1,1) = Origin(1,1);
    
    for I = 1:1:(Count(1,1) + c(1,1))
        % 第一坐标方向
        Nodes(I+(J-1)*(Count(1,1)+c(1,1)),:) = pos;
        pos(1,1) = pos(1,1) + d(1,1);     
    end % I
    
    pos(1,2) = pos(1,2) + d(1,2);
end % J
    
% =========================== 生成单元信息 ===============================

N = [0,0,0,0];
for J = 1:1:Count(1,2)
    for I = 1:1:Count(1,1)
        N = [
            I+(J-1)*(Count(1,1)+c(1,1))  ,I+(J-1)*(Count(1,1)+c(1,1))+1,...
            I+    J*(Count(1,1)+c(1,1))+1,I+    J*(Count(1,1)+c(1,1))
            ];
        if(I == Count(1,1) || J == Count(1,2))
            % 到达边界
            if(c(1,1) == 1 && c(1,2) == 1)
                % 不回环
                Elems(I+(J-1)*Count(1,1),:) = N;
            else
                if(c(1,1) == 0)
                    % 第一坐标回环
                    N(2) = N(2) - Count(1,1);
                    N(3) = N(3) - Count(1,1);
                    Elems(I+(J-1)*Count(1,1),:) = N;
                end
                
                if(c(1,2) == 0 && J == Count(1,2))
                    % 第二坐标回环
                    N(3) = N(3) - J*(Count(1,1));
                    N(4) = N(4) - J*(Count(1,1));
                    Elems(I+(J-1)*Count(1,1),:) = N;
                end
            end
        else
            Elems(I+(J-1)*Count(1,1),:) = N;
        end % I == Count(1,1) || J == Count(1,2)
    end % I
end % J



