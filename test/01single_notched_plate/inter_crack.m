function o = inter_crack(n1, n2)
%若两点连线穿过预制裂纹，inter_crack函数返回true
% 抄师兄的，能用
    ee=1e-4;
    % point=[];   %没有孔缝
    % 裂纹坐标  x1        y1        x2     y2
    point=[     0       0        500+ee   0] * 1e-3;
    ax = n1(1);
    ay = n1(2);
    bx = n2(1);
    by = n2(2);
    o = false;
    %判断连结是否和长方形边界相交
    if size(point,1)~=0
        for k=1:size(point,1)
            cx=point(k,1);
            cy=point(k,2);
            dx=point(k,3);
            dy=point(k,4);
            u=(cx-ax)*(by-ay)-(bx-ax)*(cy-ay);
            v=(dx-ax)*(by-ay)-(bx-ax)*(dy-ay);
            w=(ax-cx)*(dy-cy)-(dx-cx)*(ay-cy);
            z=(bx-cx)*(dy-cy)-(dx-cx)*(by-cy);
            if u*v<0 && w*z<0
                o = true;
                break;
            end
        end
    end
end
