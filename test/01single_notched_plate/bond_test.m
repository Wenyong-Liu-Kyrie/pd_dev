clc
clear
addpath('/home/sunhh/pd_dev/src/')

elements = readmatrix('elements');
nodes = readmatrix('nodes') / 1e3; %单位mm
mesh1 = mesh(elements, nodes);
elements_pd = elements(1:62740, :);
elements_fem = elements(62741:84680, :);
nodeid_pd = unique(elements_pd);
nodeid_fem = unique(elements_fem);
nodeid_share = intersect(nodeid_fem, nodeid_pd);
nodeid_fem_pure = setdiff(nodeid_fem, nodeid_share);
thick = 1;
fem1 = fem(elements_fem, 'q', 2, 2);
fem1.E = 30e3; %单位MPa
fem1.nu = 1/3;
fem1.thick = thick;
dx = 2 * 1e-3; %质点间距，根据画的网格调整/单位mm
delta = 3.015 * dx;
pd1 = pd(delta, dx, thick, nodeid_pd, mesh1.nodes, @inter_crack);
pd1.c = 6*fem1.E/pi/thick/delta^3/(1-fem1.nu);
G0 = 2700e-3; %断裂能密度/单位J/mm2
sc = sqrt(4*pi*G0/9/fem1.E/delta);
clear elements elements_pd nodeid_pd thick elements_fem delta dx G0

% 初始化固定和加载边界条件
fixed_dofs = [];
err = 1e-7;

% 计算上下边界的位置
bottom_boundary = nodes(:, 2) + max(nodes(:, 2));
top_boundary = nodes(:, 2) - max(nodes(:, 2));

% 确定固定和加载边界条件
fixed_bottom = abs(bottom_boundary) < err;
fixed_top = abs(top_boundary) < err;

% 添加固定边界条件
fixed_dofs = [fixed_dofs; find(fixed_bottom)*2 - 1; find(fixed_bottom)*2];
fixed_dofs = [fixed_dofs; find(fixed_top)*2];

% 添加加载边界条件
load_dofs = find(fixed_top)*2 - 1;

% 设置加载值
load_value = repmat(3e-1, length(load_dofs), 1);

solve1 = solve(fixed_dofs, load_dofs, load_value, mesh1.dofs);
clear nodes err fixed_dofs load_dofs load_value

K_fem = fem1.stiffness_assemble(mesh1.nodes, mesh1.dofs);
K_pd = pd1.stiffness_assemble(mesh1.nodes, mesh1.dofs);
K_fem(2*nodeid_share-1, :) = 0;
K_fem(2*nodeid_share, :) = 0;
K_pd(2*nodeid_fem_pure-1, :) = 0;
K_pd(2*nodeid_fem_pure, :) = 0;
solve1.K =K_fem + K_pd;
clear K_fem K_pd

nincr = 1;
miter = 10;
fileID = fopen('timeseries_data_u.txt', 'w');
fclose(fileID);  % 立即关闭文件
fileID = fopen('timeseries_data_d.txt', 'w');
fclose(fileID);  % 立即关闭文件
for incr = 1:nincr
    incr
    tic
    for iter = 1:miter
        du = solve1.displacement_load(nincr);
        pd1.bond_stretch_gen(mesh1.nodes, solve1.u + du);
        fail_index = pd1.bond_fail_judge(sc, 10000);
        if numel(fail_index) > 0
            solve1.K = pd1.stiffness_update(solve1.K, fail_index, mesh1.nodes, mesh1.dofs);
        else
            solve1.u = solve1.u + du;
            iter
            break
        end
        if iter == miter
            iter
            solve1.u = solve1.u + du;
        end
    end
    d = pd1.damage_gen(size(mesh1.nodes, 1));
    %做一些存储工作

    fileID = fopen('timeseries_data_u.txt', 'a');
    dlmwrite('timeseries_data_u.txt', solve1.u', '-append', 'delimiter', '\t', 'precision', '%.6f');
    fclose(fileID);

    fileID = fopen('timeseries_data_d.txt', 'a');
    dlmwrite('timeseries_data_d.txt', d', '-append', 'delimiter', '\t', 'precision', '%.6f');
    fclose(fileID);
    time = toc
end
