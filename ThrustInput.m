%
% Description
% ===========
%
%


% =========================== 计算常数 ====================================

% 调试开关
IS_DEBUG = 1;

% 止推瓦块的结构尺寸
PAD_SPAN = 2*pi/16; % 每个瓦块的张角
% PAD_SPAN = 2*pi/2; % 每个瓦块的张角
PAD_DIM = [
    % 起始角度 结束角度          内半径 外半径
    0.0,       PAD_SPAN+0.0,    0.1,   0.2;
    pi/2,      pi/2+PAD_SPAN,   0.1,   0.2;
    pi,        pi+PAD_SPAN,     0.1,   0.2;
    3/2*pi,    3/2*pi+PAD_SPAN, 0.1,   0.2
    ];

% 轴承参数
HP          = 0.01e-3;    % 轴承轴向基础间隙，单位m
AS_RPM      = 2100;       % 轴转速，单位rpm
AS          = AS_RPM/60*2*pi; % 轴转速，单位rad/s
VISCO       = 470e-6;     % 润滑剂动力粘度
VIS_EN      = 1.0;        % 润滑剂粘度放大系数
TURB_SWITCH = 1;          % 湍流计算开关
RHO         = 1.0e3;      % 润滑剂密度
PB          = 0.0;        % 边界压力，Pa
ALPHA       = 1e-4;       % rad, 瓦斜面偏角
W           = 1e4;        % 稳态轴向力,N

% 计算常量
ELE_NODES_NUM = 4;        % 每个单元的节点个数，四节点
NEWTON_RAPHSON_NORM = 1e-8;

STIFFNESS_BASE = 1e9;     % 用于将求得的刚度值以10^9形式显示
DAMPING_BASE   = 1e6;     % 用于将求得的阻尼值以10^6形式显示

ANG_OFF = pi; % xy坐标轴角度偏移