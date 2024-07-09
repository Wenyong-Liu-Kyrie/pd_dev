classdef pd < handle
% 均匀离散网格
    properties
        node_index   % 存储pd点在所有网格节点中的序号/键数量x1
        bond_table   % 存储所有键的节点号/键数量x2
        bond_length  % 键的长度/键数量x1
        bond_stretch % 键的伸长率/键数量x1
        bond_status  % 键连接状态/键数量x1
        bond_number  % 键数量/标量
        delta        % 近场域半径/标量
        thick        % 平面单元的厚度/标量
        r            % pd点半径/标量
        area         % pd点所占据面积/标量
        volume       % pd点所占据体积/标量
        c            % 键刚度/标量
    end

    methods
        function obj = pd(delta, dx, thick, node_index_pd, nodes, inter_crack)
        %{
          生成pd类的一个对象
          输入：
          delta ----------- 近场域半径/标量
          dx -------------- 均匀离散网格质点间距/标量
          thick ----------- 平面单元的厚度
          node_index_pd --- 部分pd点(不包含过渡区网格中的pd点)在在所有网格节点中的序号/列向量
          nodes ----------- 所有网格节点的坐标/节点数x2
          obj.bonds_gen --- 生成键表/函数
          inter_crack ----- 预制裂纹，若两点穿过裂纹，返回真/函数
        %}
            obj.delta = delta;
            obj.r = dx / 2;
            obj.thick = thick;
            obj.area = dx * dx;
            obj.volume = obj.area * thick;
            [obj.bond_table, obj.bond_number] = obj.bonds_gen(node_index_pd, nodes, delta, inter_crack);
            obj.node_index = unique(obj.bond_table);
        end
        function [bond_table, bond_number] = bonds_gen(~, node_index_pd, nodes, delta, inter_crack)
        %{
          生成键表
          输入：
          node_index_pd --- 部分pd点(不包含过渡区网格中的pd点)在在所有网格节点中的序号/点数量x1
          nodes ----------- 所有网格节点的坐标/节点数x2
          delta ----------- 近场域半径/标量
          inter_crack ----- 预制裂纹，若两点穿过裂纹，返回真/函数
          输出：
          bond_table ------ 存储所有键的节点号/键数量x2
          bond_number ----- 键数量/标量
        %}
             [idx, ~] = rangesearch(nodes, nodes(node_index_pd, :), delta);
             bond_table = zeros(30*length(node_index_pd), 2);
             bond_number = 0;
             if nargin == 5
                 % 若有inter_crack函数输入
                 for i = 1:length(idx)
                     for j = 2:length(idx{i, 1}) %idx{i, 1}的第一个值总是i自己
                         if ~inter_crack(nodes(node_index_pd(i), :), nodes(idx{i, 1}(j), :))
                             bond_number = bond_number + 1;
                             bond_table(bond_number, 1) = node_index_pd(i);
                             bond_table(bond_number, 2) = idx{i, 1}(j);
                         end
                     end
                 end
             else
                 %若无需inter_crack函数输入
                 for i = 1:length(idx)
                     for j = 2:length(idx{i, 1}) %idx{i, 1}的第一个值总是i自己
                         bond_number = bond_number + 1;
                         bond_table(bond_number, 1) = node_index_pd(i);
                         bond_table(bond_number, 2) = idx{i, 1}(j);
                     end
                 end
             end
             bond_table = bond_table(1:bond_number, :);
             bond_table=sort(bond_table,2);
             bond_table=unique(bond_table,'rows');
        end
        function K = stiffness_assemble(obj, nodes, dofs)
        %{
          组装pd整体刚度矩阵
          输入：
          dofs --------------------- 所有网格(包括有限元部分)的总自由度数目/标量
          nodes -------------------- 所有网格节点的坐标/节点数x2
          obj.bond_number ---------- 键数量/标量
          obj.bond_stiffness_gen --- 生成键刚度矩阵/函数
          输出：
          K ------------------------ pd整体刚度矩阵/dofs x dofs
        %}
            Is = zeros(16*obj.bond_number, 1);
            Js = zeros(16*obj.bond_number, 1);
            Vs = zeros(16*obj.bond_number, 1);
            for i = 1:obj.bond_number
                [I, J, V] = obj.bond_stiffness_gen(i, nodes);
                Is(16*(i-1)+1:16*i) = I;
                Js(16*(i-1)+1:16*i) = J;
                Vs(16*(i-1)+1:16*i) = V;
            end
            K = sparse(Is, Js, Vs, dofs, dofs);
        end
        function [I, J, V] = bond_stiffness_gen(obj, bond_index, nodes)% TODO
        %{
          生成键刚度矩阵
          输入：
          bond_index --- 键序号/标量
          nodes -------- 所有网格节点的坐标/节点数x2
          obj,delta ---- 近场域半径/标量
          obj.r -------- pd点半径/标量
          obj.volume --- pd点所占据体积/标量
          obj.c -------- 键刚度/标量
          输出：
          I ------------ 键刚度矩阵在整刚中位置的横坐标/列向量
          J ------------ 键刚度矩阵在整刚中位置的纵坐标/列向量
          V ------------ 键刚度矩阵的值/列向量
        %}
        end

        function bond_stretch_gen(obj, nodes, u)% TODO
        %{
          计算键伸长率
          输入：
          nodes -------------- 所有网格节点的坐标/节点数x2
          u ------------------ 所有节点自由度的位移/总自由度数目x1
          obj.bond_table ----- 存储所有键的节点号/键数量x2
          obj.bond_number ---- 键数量/标量
          obj.bond_length ---- 键的长度/键数量x1
          修改：
          obj.bond_stretch --- 键的伸长率/键数量x1
        %}
        end
    end
end
