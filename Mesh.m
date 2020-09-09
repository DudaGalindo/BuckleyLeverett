% Arquivo criado para gerar a malha (gerada manualmente) - só funciona para
% malhas 1-D.
classdef Mesh
    properties
        n_el 
        L 
    end
    methods 
        function conec = faces_conec(obj)
            conec = zeros(obj.n_el-1,2);
            conec(:,1) = linspace(1,obj.n_el-1, obj.n_el-1);
            conec(:,2) = linspace(2,obj.n_el, obj.n_el-1);
        end
        function vols_ID = all_vols_ID(obj)
            vols_ID = linspace(1,obj.n_el, obj.n_el);
        end
        function contour_vols_ID = contour_vols_ID(obj)
            contour_vols_ID = [1 obj.n_el];
        end
        function vols_coord = vols_coord(obj)
            vols_coord = linspace(0+obj.L/(2*obj.n_el), 1 - obj.L/(2*obj.n_el),obj.n_el); 
        end
        function faces_ID = internal_faces_ID(obj)
            faces_ID = linspace(1,obj.n_el);
        end
    end
end