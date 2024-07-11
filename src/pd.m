classdef pd < handle
% 均匀离散网格
    properties
        node_index   % 存储pd点在所有网格节点中的序号/键数量x1
        bond_table   % 存储所有键的节点号/键数量x2
        bond_length  % 键的长度/键数量x1
        bond_stretch % 键的伸长率/键数量x1
        bond_status  % 键连接状态(true表示完好)/键数量x1
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
            if nargin == 6
                % 若有inter_crack函数输入
                obj.bond_table = obj.bonds_gen(node_index_pd, nodes, delta, inter_crack);
            else
                %若无需inter_crack函数输入
                obj.bond_table = obj.bonds_gen(node_index_pd, nodes, delta);
                end
            obj.bond_number = size(obj.bond_table, 1);
            obj.bond_status = true(obj.bond_number, 1);
            obj.bond_stretch = zeros(obj.bond_number, 1);
            obj.bond_length_gen(nodes);
            obj.node_index = unique(obj.bond_table);
        end
        function bond_table = bonds_gen(~, node_index_pd, nodes, delta, inter_crack)
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
             m = 0;
             if nargin == 5
                 % 若有inter_crack函数输入
                 for i = 1:length(idx)
                     for j = 2:length(idx{i, 1}) %idx{i, 1}的第一个值总是i自己
                         if ~inter_crack(nodes(node_index_pd(i), :), nodes(idx{i, 1}(j), :))
                             m = m + 1;
                             bond_table(m, 1) = node_index_pd(i);
                             bond_table(m, 2) = idx{i, 1}(j);
                         end
                     end
                 end
             else
                 %若无需inter_crack函数输入
                 for i = 1:length(idx)
                     for j = 2:length(idx{i, 1}) %idx{i, 1}的第一个值总是i自己
                         m = m + 1;
                         bond_table(m, 1) = node_index_pd(i);
                         bond_table(m, 2) = idx{i, 1}(j);
                     end
                 end
             end
             bond_table = bond_table(1:m, :);
             bond_table=sort(bond_table,2);
             bond_table=unique(bond_table,'rows');
        end
        function bond_length_gen(obj, nodes)
        %{
          生成所有初始键长
          输入：
          nodes ------------- 所有网格节点的坐标/节点数x2
          obj.bond_table ---- 存储所有键的节点号/键数量x2
          obj.bond_number --- 键数量/标量
          修改：
          obj.bond_length --- 键的长度/键数量x1
        %}

            obj.bond_length = zeros(obj.bond_number, 1);% 预分配键长向量
            % 获取键表中的节点索引
            nodei = obj.bond_table(:, 1);
            nodej = obj.bond_table(:, 2);

            % 提取节点坐标
            nodeix = nodes(nodei, 1);
            nodeiy = nodes(nodei, 2);
            nodejx = nodes(nodej, 1);
            nodejy = nodes(nodej, 2);

            % 计算节点坐标差值
            dx = nodejx - nodeix;
            dy = nodejy - nodeiy;

            % 计算键的长度
            obj.bond_length = sqrt(dx.^2 + dy.^2);
        end
        function K = stiffness_assemble(obj, nodes, dofs)
        %{
          组装pd整体刚度矩阵
          输入：
          nodes -------------------- 所有网格节点的坐标/节点数x2
          dofs --------------------- 所有网格(包括有限元部分)的总自由度数目/标量
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
        function K = stiffness_update(obj, K, fail_index, nodes, dofs)
        %{
          断键之后更新整体刚度矩阵
          输入：
          K ------------ 更新前整刚/dofs x dofs
          fail_index --- 断键的键的序号/断键数量x1
          nodes -------- 所有网格节点的坐标/节点数x2
          dofs  -------- 所有网格(包括有限元部分)的总自由度数目/标量
          输出：
          K ------------ 更新后的整刚/dofs x dofs
        %}
            fail_number = numel(fail_index);
            Is = zeros(16*fail_number, 1);
            Js = zeros(16*fail_number, 1);
            Vs = zeros(16*fail_number, 1);
            for i = 1:fail_number
                [I, J, V] = obj.bond_stiffness_gen(fail_index(i), nodes);
                Is(16*(i-1)+1:16*i) = I;
                Js(16*(i-1)+1:16*i) = J;
                Vs(16*(i-1)+1:16*i) = V;
            end
            K = K - sparse(Is, Js, Vs, dofs, dofs);
        end
        function [I, J, V] = bond_stiffness_gen(obj, bond_index, nodes)
        %{
          生成单键刚度矩阵
          输入：
          bond_index -------- 键序号/标量
          nodes ------------- 所有网格节点的坐标/节点数x2
          obj.bond_table ---- 存储所有键的节点号/键数量x2
          obj.bond_length --- 键的长度/键数量x1
          obj.delta --------- 近场域半径/标量
          obj.r ------------- pd点半径/标量
          obj.volume -------- pd点所占据体积/标量
          obj.c ------------- 键刚度/标量
          输出：
          I ----------------- 键刚度矩阵在整刚中位置的横坐标/列向量
          J ----------------- 键刚度矩阵在整刚中位置的纵坐标/列向量
          V ----------------- 键刚度矩阵的值/列向量
        %}
            nodei=obj.bond_table(bond_index, 1);   %键端点的节点数
            nodej=obj.bond_table(bond_index, 2);

            nodeix=nodes(nodei, 1);    %节点的坐标
            nodeiy=nodes(nodei, 2);
            nodejx=nodes(nodej, 1);
            nodejy=nodes(nodej, 2);

            %% 进行体积修正
            volFac=0;
            if obj.bond_length(bond_index)<=(obj.delta-obj.r)
                volFac=1;
            elseif (obj.bond_length(bond_index)>(obj.delta-obj.r)) && (obj.bond_length(bond_index)<=obj.delta)
                volFac=(obj.delta+obj.r-obj.bond_length(bond_index))/(2*obj.r);
            else
                volFac=0;
            end
            %% 求三角函数值
            cosine=(nodejx-nodeix)/obj.bond_length(bond_index);
            sine=(nodejy-nodeiy)/obj.bond_length(bond_index);

            V = [cosine^2 cosine*sine -cosine^2 -cosine*sine
                 cosine*sine sine^2 -cosine*sine -sine^2
                 -cosine^2 -cosine*sine cosine^2 cosine*sine
                 -cosine*sine -sine^2 cosine*sine sine^2]*obj.c*volFac* obj.volume^2/obj.bond_length(bond_index);

            %% 输出结果
            fod_nodeij=zeros(4,1);
            fod_nodeij(1)=2*nodei-1;
            fod_nodeij(2)=2*nodei;
            fod_nodeij(3)=2*nodej-1;
            fod_nodeij(4)=2*nodej;

            I=repmat(fod_nodeij,4,1);
            J=repelem(fod_nodeij,4);
            V =V(:);
        end
        function bond_stretch_gen(obj, nodes, u)
        %{
          计算所有键伸长率
          输入：
          nodes -------------- 所有网格节点的坐标/节点数x2
          u ------------------ 所有节点自由度的位移/总自由度数目x1
          obj.bond_table ----- 存储所有键的节点号/键数量x2
          obj.bond_number ---- 键数量/标量
          obj.bond_length ---- 键的长度/键数量x1
          修改：
          obj.bond_stretch --- 键的伸长率/键数量x1
        %}
            for i = 1:obj.bond_number
                node1 = obj.bond_table(i, 1); %节点1
                node2 = obj.bond_table(i, 2); %节点2
                n1 = (nodes(node2, 1) - nodes(node1, 1)) / obj.bond_length(i); %单位向量第一个值
                n2 = (nodes(node2, 2) - nodes(node1, 2)) / obj.bond_length(i); %单位向量第二个值
                ux1 = u(2*node1-1);
                uy1 = u(2*node1);
                ux2 = u(2*node2-1);
                uy2 = u(2*node2);
                obj.bond_stretch(i) = ((ux2-ux1)*n1+(uy2-uy1)*n2)/obj.bond_length(i);
            end
        end
        function fail_index = bond_fail_judge(obj, sc, max_fail_number)
        %{
          通过键的伸长率判断是否断键
          输入：
          sc ----------------- 键极限伸长率/标量
          max_fail_number ---- 最大断键数量/标量
          obj.bond_stretch --- 键的伸长率/键数量x1
          obj.bond_number ---- 键数量/标量
          输出：
          fail_index --------- 本回合断键的键的序号/断键数量(不超过最大断键数量)x1
          修改：
          obj.bond_status ---- 键连接状态(true表示完好)/键数量x1
        %}
        % 创建包含键序号和伸长率的矩阵
            bond_stretch = [(1:obj.bond_number)', obj.bond_stretch];

            % 找到未断键的键
            valid_bonds = bond_stretch(obj.bond_status, :);

            % 找到可能断键的键
            potential_failures = valid_bonds(valid_bonds(:, 2) > sc, :);

            % 根据可能断键的数量确定断键序号
            if size(potential_failures, 1) <= max_fail_number
                fail_index = potential_failures(:, 1);
            else
                % 按伸长率降序排序，取前 max_fail_number 个键
                sorted_failures = sortrows(potential_failures, -2);
                fail_index = sorted_failures(1:max_fail_number, 1);
            end

            % 更新断键键的连接状态为 false
            obj.bond_status(fail_index) = false;
        end
        function d = damage_gen(obj, node_number)
            %{
              计算pd点损伤值
              输入：
              node_number ------- 所有节点数量
              obj.node_index ---- 存储pd点在所有网格节点中的序号/键数量x1
              obj.bond_table ---- 存储所有键的节点号/键数量x2
              obj.bond_status --- 键连接状态(true表示完好)/键数量x1
              输出：
              d ----------------- 模型所有点的损伤值/总节点数x1
            %}
            d = zeros(node_number, 1);

            % 获取节点索引和键表
            node_index = obj.node_index;
            bond_table = obj.bond_table;
            bond_status = obj.bond_status;

            % 先找出所有断键相关点
            bond_failed = bond_table(~bond_status, :);
            node_failed = unique(bond_failed);

            % 计算节点的损伤值
            for k = 1:length(node_failed)
                i = node_failed(k);  % 获取节点编号

                % 找到所有与节点 i 相关的键的索引
                idx = (bond_table(:, 1) == i | bond_table(:, 2) == i);

                % 计算与节点 i 相关的键的连接状态
                relevant_bonds = bond_status(idx);

                % 计算节点 i 的损伤值
                if ~isempty(relevant_bonds)
                    d(i) = 1 - sum(relevant_bonds) / sum(idx);
                else
                    d(i) = 1;  % 如果没有相关键，则质点完好
                end
            end
        end
    end
end
