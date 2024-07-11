classdef mesh < handle
    properties
        elements    %单元节点号/单元个数x单元结点数
        nodes       %节点坐标（模型所有节点）/结点个数x单点自由度数
        node_number %节点数量/标量
        dimension    %2或3/标量
        dofs        %number of degrees of freedom/标量
    end
    methods
        function obj = mesh(elements, nodes)
        %{
          生成mesh类的一个对象
          输入：
          elements --- 单元节点号/单元个数x单元结点数
          nodes ------ 节点坐标/结点个数x单点自由度数
        %}
            obj.elements = elements;
            obj.nodes = nodes;
            [obj.node_number, obj.dimension] = size(nodes);
            obj.dofs = obj.node_number * obj.dimension;
        end
        function plot_mesh(obj, fn)
        %{
          绘制网格图像
          输入：
          fn --- 输出文件名(不包含后缀)/字符串
          输出：
          fn.png
        %}
            figure
            patch('Vertices',obj.nodes,'Faces',obj.elements,'FaceVertexCData',zeros(obj.node_number, 1),'FaceColor','[0.6 0.6 0.7]','EdgeColor','blue');
            view(0, 90)
            axis equal
            axis on %显示坐标系
            saveas(gcf,[fn,'.png']);
        end
        function write_vtk(obj, fn, u, d, is_binary)

            fid = fopen(fn, 'w');% 打开文件

            % 写入VTK文件头
            fprintf(fid, '# vtk DataFile Version 3.0\n');
            fprintf(fid, 'VTK from MATLAB\n');
            fprintf(fid, 'BINARY\n');
            fprintf(fid, 'DATASET UNSTRUCTURED_GRID\n');

            % 写入点数据
            num_nodes = size(obj.nodes, 1);
            fprintf(fid, 'POINTS %d float\n', num_nodes);

            % 将点数据转换为大端字节顺序并写入文件
            fwrite(fid, swapbytes(single([obj.nodes, zeros(num_nodes, 1)]')),'float32');

            % 写入单元数据
            num_elements = size(obj.elements, 1);
            num_indices = numel(obj.elements);
            fprintf(fid, 'CELLS %d %d\n', num_elements, num_indices + num_elements);

            % 准备单元数据并写入文件
            cell_data = [repmat(uint32(4), num_elements, 1), uint32(obj.elements - 1)]';
            fwrite(fid, swapbytes(cell_data), 'uint32');

            % 写入单元类型 (9代表四边形)
            fprintf(fid, 'CELL_TYPES %d\n', num_elements);
            cell_types = repmat(uint32(9), num_elements, 1);
            fwrite(fid, swapbytes(cell_types), 'uint32');

            % 写入点数据 - 位移
            fprintf(fid, 'POINT_DATA %d\n', num_nodes);
            fprintf(fid, 'VECTORS displacements float\n');

            % 准备位移数据并写入文件
            displacements = [u(1:2:end), u(2:2:end), zeros(num_nodes, 1)]';
            fwrite(fid, swapbytes(single(displacements)), 'float32');

            % 写入点数据 - 损伤
            fprintf(fid, 'SCALARS damage float 1\n');
            fprintf(fid, 'LOOKUP_TABLE default\n');

            % 将损伤数据转换为大端字节顺序并写入文件
            fwrite(fid, swapbytes(single(d')), 'float32');

            % 关闭文件
            fclose(fid);

            disp('二进制VTK文件已成功生成');

        end
        function write_vtk_binary_time_series(obj)
        end
        function plot_variable(obj, d, fn)
        %{
          绘制场变量d的图像
          输入：
          d ---- 定义在节点上的标量场变量d
          fn --- 输出文件名(不包含后缀)/字符串
          输出：
          fn.png
        %}
            figure
            patch('Vertices',obj.nodes,'Faces',obj.elements,'FaceVertexCData',d,'FaceColor','[0.6 0.6 0.7]','EdgeColor','none');
            view(0, 90)
            shading('interp');
            colormap(jet)
            colorbar;
            axis equal
            axis off %不显示坐标系
            saveas(gcf,[fn,'.png']);
        end
    end
end
