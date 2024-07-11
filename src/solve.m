classdef  solve < handle
    properties
        free_dofs  %自由的自由度/mx1
        fixed_dofs %约束自由度/nx1
        load_dofs  %加载自由度/vx1
        load_value %加载自由度上的荷载值/vx1
        dofs       %总自由度数目/标量(m+n+v)
        u          % displacement vector/(m+n+v)x1
        f_ext      % external force vector/(m+n+v)x1
        f_int      % internal force vector/(m+n+v)x1
        R          % residual/(m+n+v)x1
        K          % global stiffness matrix/(m+n+v)x1
    end
    methods
        function obj = solve(fixed_dofs, load_dofs, load_value, dofs)
        %{
          生成solve类的一个对象
          输入：
          free_dofs ---- 自由的自由度/mx1
          fixed_dofs --- 约束自由度/nx1
          load_dofs ---- 加载自由度/vx1
          load_value --- 加载自由度上的荷载值/vx1
          dofs --------- 总自由度数目/标量(m+n+v)
        %}
            obj.fixed_dofs = fixed_dofs;
            obj.load_dofs = load_dofs;
            obj.free_dofs = setdiff((1:dofs)', union(fixed_dofs,load_dofs));
            obj.load_value = load_value;
            obj.dofs = dofs;
            obj.u = zeros(dofs, 1);
            obj.f_ext = zeros(dofs, 1);
            obj.f_int = zeros(dofs, 1);
            obj.R = zeros(dofs, 1);
        end
        function du = displacement_load(obj, nincr)
        %{
          位移荷载分nincr次加载，加载一步
          输入：
          nincr --- 荷载分成nincr次加载/标量
          输出：
          du ------ 位移增量/总自由度x1
        %}
            if (nargin == 1)
                %如果没有输入加载步数，默认分一步加载
                nincr = 1;
            end
            du = zeros(obj.dofs, 1);
            du(obj.free_dofs) = obj.K(obj.free_dofs,obj.free_dofs)\...
                (obj.R(obj.free_dofs) - obj.K(obj.free_dofs, obj.load_dofs) * obj.load_value / nincr);
            du(obj.load_dofs) = obj.load_value / nincr;
        end
    end
end
