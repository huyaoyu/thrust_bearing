function [ Nodes,Elems ] = RecField2DIsoWrap(Origin,Span,Count,WrapCheck)
% ��άƽ���������ʹ��˫�����Ľڵ�Ȳ�������Wrap�ĺ����ǽ�ĳЩ����ĵ�һ
% �����һ��������ȣ��γɻ�״�Ľڵ�����ֲ�
% Origin��һ����ά����������¼ԭ��λ��
% Span��һ����ά����������¼ƽ�����������귽��ĳ���
% Count��һ����ά��������¼ÿ�����귽���ϵĵ�Ԫ����
% WrapCheck��һ����ά������������Ҫ�ػ������귽��
% Nodes��һ������ÿ�м�¼һ������������꣬��������
% Elems�� һ������ÿ�м�¼һ����Ԫ���ĸ��ڵ�����������ַ���
% ����Nodes���Ե�һ����Ϊ���ȵģ���һ������������Ҳ���ǵ�һ�������ȸ���
%
% 2012.3.24 ����
%
% 2012.3.31 ����
% ������Origin��ʹ��

% ========================= ���������� ================================

FUN_DISP_PREFIX = 'FUNC RecField2DIso: ';
NODES_PER_ELEM  = 4; % ÿ����Ԫ�Ľڵ����
COORD_PER_NODE  = 2; % ÿ���ڵ�Ķ����������

% ========================== ��ʼ��� =================================

% ��� Origin
[row,col] = size(Origin);
if(row ~= 1 || col ~= COORD_PER_NODE)
    Nodes = 0;
    Elems = 0;
    disp([FUN_DISP_PREFIX,'The dimension of Origin is wrong!']);
    return
end

% ��� Span
[row,col] = size(Span);
if(row ~= 1 || col ~= COORD_PER_NODE)
    Nodes = 0;
    Elems = 0;
    disp([FUN_DISP_PREFIX,'The dimension of Span is wrong!']);
    return
end

% ��� Count
[row,col] = size(Count);
if(row ~= 1 || col ~= COORD_PER_NODE)
    Nodes = 0;
    Elems = 0;
    disp([FUN_DISP_PREFIX,'The dimension of Count is wrong!']);
    return
end

% ��� WrapCheck
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

% ============================== �������� ===============================

% �鿴�ػ�,���ػ������귽��
c = (WrapCheck == 0);

NUM_NODES = (Count(1,1)+c(1,1))*(Count(1,2)+c(1,2)); % �ܵ���
NUM_ELEMS = Count(1,1)*Count(1,2);         % �ܵ�Ԫ��

% Ԥ���䷵��ֵ�ڴ�
Nodes = zeros(NUM_NODES,COORD_PER_NODE);
Elems = zeros(NUM_ELEMS,NODES_PER_ELEM);

% =========================== ���ɽڵ���Ϣ ===============================

d = Span ./ Count; % ����������������
pos = Origin; % �ڲ�ѭ���ĵ�ǰ����λ��


for J = 1:1:(Count(1,2) + c(1,2))
    % �ڶ����귽��
    
    % ����ǰ��һ����λ��
    pos(1,1) = Origin(1,1);
    
    for I = 1:1:(Count(1,1) + c(1,1))
        % ��һ���귽��
        Nodes(I+(J-1)*(Count(1,1)+c(1,1)),:) = pos;
        pos(1,1) = pos(1,1) + d(1,1);     
    end % I
    
    pos(1,2) = pos(1,2) + d(1,2);
end % J
    
% =========================== ���ɵ�Ԫ��Ϣ ===============================

N = [0,0,0,0];
for J = 1:1:Count(1,2)
    for I = 1:1:Count(1,1)
        N = [
            I+(J-1)*(Count(1,1)+c(1,1))  ,I+(J-1)*(Count(1,1)+c(1,1))+1,...
            I+    J*(Count(1,1)+c(1,1))+1,I+    J*(Count(1,1)+c(1,1))
            ];
        if(I == Count(1,1) || J == Count(1,2))
            % ����߽�
            if(c(1,1) == 1 && c(1,2) == 1)
                % ���ػ�
                Elems(I+(J-1)*Count(1,1),:) = N;
            else
                if(c(1,1) == 0)
                    % ��һ����ػ�
                    N(2) = N(2) - Count(1,1);
                    N(3) = N(3) - Count(1,1);
                    Elems(I+(J-1)*Count(1,1),:) = N;
                end
                
                if(c(1,2) == 0 && J == Count(1,2))
                    % �ڶ�����ػ�
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



