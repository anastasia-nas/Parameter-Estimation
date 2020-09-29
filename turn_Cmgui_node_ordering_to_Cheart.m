function [Reordered_Elements_Cubic_ala_cheart,error]=turn_Cmgui_node_ordering_to_Cheart(Elements_Cubic_ala_cmgui,Nodes_per_elem_dir)

% the function reorders cubic lagrange elemnts so that they have the same
% node ordering as in cmgui (no priority to corner nodes )

%Elements_Cubic: ola ta Elements tou mesh sou: rows: element no, columns:
%nodes as ordered in element

% Nodes_per_Elem_dir: nodes per ksi directions in each element: e.g. 4 for
% cubic interpolation..
% cheart_order=[1,4,13,16,49,52,61,64,2,3,5:12,14,15,17:48,50,51,53:60,62,63];
% I was using cheart_order input explicitly as I was always using cub
% lagrange meshes, now trying to use some Quad Lagrange for computational
% economy (cheart crashes) and I do this iteratively. I checked with linear
% and cubic meshes and my "iterative" cheart order is correct..

% first store the corner nodes
cheart_order=[];
% corner nodes go first::
for n_ksi3=1:Nodes_per_elem_dir
    for n_ksi2=1:Nodes_per_elem_dir
        for n_ksi1=1:Nodes_per_elem_dir
            if (n_ksi1==1 || n_ksi1==Nodes_per_elem_dir) && (n_ksi2==1 || n_ksi2==Nodes_per_elem_dir) && (n_ksi3==1 || n_ksi3==Nodes_per_elem_dir)
                cheart_order=[cheart_order,(n_ksi3-1)*Nodes_per_elem_dir^2+(n_ksi2-1)*Nodes_per_elem_dir+n_ksi1];
            end
            
        end
    end
end
% rest of nodes later::
for n_ksi3=1:Nodes_per_elem_dir
    for n_ksi2=1:Nodes_per_elem_dir
        for n_ksi1=1:Nodes_per_elem_dir
            if (n_ksi1==1 || n_ksi1==Nodes_per_elem_dir) && (n_ksi2==1 || n_ksi2==Nodes_per_elem_dir) && (n_ksi3==1 || n_ksi3==Nodes_per_elem_dir)
                % do nothing
            else
                cheart_order=[cheart_order,(n_ksi3-1)*Nodes_per_elem_dir^2+(n_ksi2-1)*Nodes_per_elem_dir+n_ksi1];
            end
            
        end
    end
end

% no need for this any more
% if Nodes_per_elem_dir~=4
%     error=1;
% else
%     error=0;
% end
error=0; % I leave it like this for compatibility with earlier scripts.
Reordered_Elements_Cubic_ala_cheart=Elements_Cubic_ala_cmgui(:,cheart_order);