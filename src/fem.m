classdef  fem < handle
    properties
        elements             % 有限元单元集合/单元个数x单元结点数
        dof_per_element      % 单个单元的自由度数目/标量
        number_of_quadrature % 等参单元积分点数量/标量
        ref_coords           % 等参单元积分点坐标/等参单元积分点数量x2
        weights              % 等参单元积分点权重/等参单元积分点数量x1
        E                    % 杨氏模量/标量
        nu                   % 泊松比/标量
        thick                % 厚度/标量
    end
    methods
        function obj = fem(elements, filter, number1, number2, number3)
        %{
          生成fem类的一个对象
          输入：
          elements --- 有限元单元集合/单元个数x单元结点数
          filter ----- 代表单元类型/字符
          number1 ---- x方向积分点数量/标量
          number2 ---- y方向积分点数量/标量
          number3 ---- z方向积分点数量/标量
        %}
            obj.elements = elements;
            obj.dof_per_element = 2*size(elements, 2);
            if strcmp(filter, 'q')
                [obj.number_of_quadrature, obj.ref_coords, obj.weights] = obj.quadrature_gauss(number1, number2);
            elseif strcmp(filter, 'H')
                [obj.number_of_quadrature, obj.ref_coords, obj.weights] = obj.quadrature_gauss(number1, number2, number3);
            elseif strcmp(filter, 't')
                [obj.number_of_quadrature, obj.ref_coords, obj.weights] = obj.quadrature_triangle(number1);
            elseif strcmp(filter, 'T')
                [obj.number_of_quadrature, obj.ref_coords, obj.weights] = obj.quadrature_tetrahedron(number1);
            end
        end
        function [number_of_quadrature, ref_coords, weights] = quadrature_gauss(obj, number_of_quadrature_x, number_of_quadrature_y, number_of_quadrature_z)
        %{
          生成四边形及六面体网格高斯积分点信息
          输入：
          number_of_quadrature_x --- x方向积分点数量/标量
          number_of_quadrature_y --- y方向积分点数量/标量
          number_of_quadrature_z --- z方向积分点数量/标量
          输出：
          number_of_quadrature ----- 等参单元积分点数量
          ref_coords --------------- 等参单元积分点坐标
          weights -------------------等参单元积分点权重
        %}
            if nargin == 2 %这种情况是不是要删掉
                number_of_quadrature = number_of_quadrature_x;
                [ref_coords, weights] = GaussLegendre(number_of_quadrature_x);
            elseif nargin == 3
                number_of_quadrature = number_of_quadrature_x * number_of_quadrature_y;
                [coords_x, weights_x] = GaussLegendre(number_of_quadrature_x);
                [coords_y, weights_y] = GaussLegendre(number_of_quadrature_y);
                ref_coords = [repelem(coords_x, numel(coords_y)) repmat(coords_y, numel(coords_x), 1)];
                weights = repelem(weights_x, numel(weights_y)) .* repmat(weights_y, numel(weights_x), 1);
            elseif nargin == 4
                number_of_quadrature = number_of_quadrature_x * number_of_quadrature_y * number_of_quadrature_z;
                [coords_x, weights_x] = GaussLegendre(number_of_quadrature_x);
                [coords_y, weights_y] = GaussLegendre(number_of_quadrature_y);
                [coords_z, weights_z] = GaussLegendre(number_of_quadrature_z);
                ref_coords = [repelem(coords_x, numel(coords_y)*numel(coords_z)) repmat(repelem(coords_y, numel(coords_z)), numel(coords_x), 1) repmat(coords_z, numel(coords_x)*numel(coords_y), 1)];
                weights = repelem(weights_x, numel(weights_y)*numel(weights_z)) .* repmat(repelem(weights_y, numel(weights_z)), numel(weights_x), 1) .* repmat(weights_z, numel(weights_x)*numel(weights_y), 1);
            end
            function [x, w] = GaussLegendre(n)
            % 网上抄的[pasuka](http://bbs.fcode.cn/thread-97-1-1.html)
            % As you can see, the code consists of 2 blocks:
            % 1: construct a symmetrical companion matrix
            % 2: determine the (real) eigenvalues (i.e. the roots of the polynomial).
            % It can produce the correct abscissas and weights, for any value n>=2.
            %
            % input: n ---- number of intergrate points
            %
            % output: x --- GaussLegendre intergrate point
            %         w --- GaussLegendre intergrate weight
            %
                i = 1:n-1;
                a = i./sqrt(4*i.^2-1);
                CM = diag(a,1) + diag(a,-1);
                [V, L] = eig(CM);
                [x, ind] = sort(diag(L));
                V = V(:,ind)';
                w = 2 * V(:,1).^2;
                return
            end
        end
        function [number_of_quadrature, ref_coords, weights] = quadrature_triangle(obj, rank)
        %{
          生成三角形单元积分点信息，数据来自朱伯芳《有限单元法原理与应用（第四版）》P168
          输入：
          rank ------------------- 积分阶数
          输出：
          number_of_quadrature --- 等参单元积分点数量
          ref_coords ------------- 等参单元积分点坐标
          weights ---------------- 等参单元积分点权重
        %}
            if rank == 1
                number_of_quadrature = 1;
                ref_coords = [1/3  1/3  1/3];
                weights = 0.5;
            elseif rank == 2
                number_of_quadrature = 3;
                ref_coords = [0.5 0   0.5
                              0.5 0.5 0
                              0   0.5 0.5];
                weights = [1/6 1/6 1/6];
            elseif rank == 3
                number_of_quadrature = 7;
                ref_coords = [1/3 0.5 0   0.5 1 0 0
                              1/3 0.5 0.5 0   0 1 0
                              1/3 0   0.5 0.5 0 0 1];
                weights = [0.225 1/15 1/15 1/15 0.025 0.025 0.025];
            end
        end
        function [number_of_quadrature, ref_coords, weights] = quadrature_tetrahedron(obj, rank)
        %{
          生成四面体单元积分点信息，数据来自朱伯芳《有限单元法原理与应用（第四版）》P168
          输入：
          rank ------------------- 积分阶数
          输出：
          number_of_quadrature --- 等参单元积分点数量
          ref_coords ------------- 等参单元积分点坐标
          weights ---------------- 等参单元积分点权重
        %}
            if rank == 1
                number_of_quadrature = 1;
                ref_coords = [0.25 0.25 0.25 0.25];
                weights = 1;
            elseif rank == 2
                number_of_quadrature = 4;
                ref_coords = [0.5854102 0.1381966 0.1381966 0.1381966
                              0.1381966 0.5854102 0.1381966 0.1381966
                              0.1381966 0.1381966 0.5854102 0.1381966
                              0.1381966 0.1381966 0.1381966 0.5854102];
                weights = [0.25 0.25 0.25 0.25];
            elseif rank == 3
                number_of_quadrature = 5;
                ref_coords = [0.25 1/3 1/6 1/6 1/6
                              0.25 1/6 1/3 1/6 1/6
                              0.25 1/6 1/6 1/3 1/6
                              0.25 1/6 1/6 1/6 1/3];
                weights = [0.8 0.45 0.45 0.45 0.45];
            end
        end
        function K = stiffness_assemble(obj, nodes, dofs)
        %{
          组装fem整体刚度矩阵
          输入：
          nodes ----------------------- 所有网格节点的坐标/节点数x2
          dofs ------------------------ 所有网格(包括有限元部分)的总自由度数目/标量
          obj.dof_per_element --------- 单个单元的自由度数目/标量
          obj.element_stiffness_gen --- 生成单元刚度矩阵/函数
          输出：
          K --------------------------- fem整体刚度矩阵/dofs x dofs
        %}
            n = obj.dof_per_element^2;
            element_number = size(obj.elements, 1);
            Is = zeros(n*element_number, 1);
            Js = zeros(n*element_number, 1);
            Vs = zeros(n*element_number, 1);
            for i = 1:element_number
                [I, J, V] = obj.element_stiffness_gen(i, nodes);
                Is(n*(i-1)+1:n*i) = I;
                Js(n*(i-1)+1:n*i) = J;
                Vs(n*(i-1)+1:n*i) = V;
            end
            K = sparse(Is, Js, Vs, dofs, dofs);
        end
        function [I, J, V] = element_stiffness_gen(obj, element_index, nodes)
        %{
          生成单元刚度矩阵
          输入：
          element_index -------------- 单元序号/标量
          nodes ---------------------- 所有网格节点的坐标/节点数x2
          obj.E ---------------------- 杨氏模量/标量
          obj.nu --------------------- 泊松比/标量
          obj.number_of_quadrature --- 等参单元积分点数量/标量
          obj.thick ------------------ 单元厚度
          obj.matrix_B_gen ----------- 生成B矩阵和Jacobi行列式/函数
          输出：
          I -------------------------- 单元刚度矩阵在整刚中位置的横坐标/列向量
          J -------------------------- 单元刚度矩阵在整刚中位置的纵坐标/列向量
          V -------------------------- 单元刚度矩阵的值/列向量
        %}
            node_index = obj.elements(element_index, :);
            coords = nodes(node_index, :);
            D = [1 obj.nu 0
                 obj.nu 1 0
                 0 0 (1-obj.nu)/2]*obj.E/(1-obj.nu^2);
            V = zeros(8, 8);
            for i = 1:obj.number_of_quadrature
                [B, Jacobi_det] = obj.matrix_B_gen(i, coords);
                V = V + obj.weights(i) * B' * D * B * Jacobi_det * obj.thick;
            end
            dof(2:2:8) = 2 * node_index;
            dof(1:2:8) = 2 * node_index - 1;
            dof = dof';
            I =repmat(dof, 8, 1);
            J = repelem(dof, 8);
            V = V(:);
        end
        function [B, Jacobi_det] = matrix_B_gen(obj, quadrature_index, coords)
        %{
          生成B矩阵和Jacobi行列式(平面四边形单元)
          输入：
          quadrature_index --- 积分点序号/标量
          coords ------------- 积分点所处单元的节点坐标/4x2
          obj.ref_coords ----- 等参单元积分点坐标/等参单元积分点数量x2
          输出：
          B ------------------ B矩阵/3x8
          Jacobi_det --------- Jacobi矩阵行列式
        %}
            xi = obj.ref_coords(quadrature_index, 1);
            eta = obj.ref_coords(quadrature_index, 2);
            N_xieta = [eta-1  1-eta 1+eta -1-eta
                       xi-1 -1-xi  1+xi   1-xi] / 4;
            J = N_xieta * coords;
            N_xy = J \ N_xieta;
            B = [N_xy(1,1) 0 N_xy(1,2) 0 N_xy(1,3) 0 N_xy(1,4) 0
                 0 N_xy(2,1) 0 N_xy(2,2) 0 N_xy(2,3) 0 N_xy(2,4)
                 N_xy(2,1) N_xy(1,1) N_xy(2,2) N_xy(1,2) N_xy(2,3) N_xy(1,3) N_xy(2,4) N_xy(1,4)];
            Jacobi_det = det(J);
        end
    end
end
