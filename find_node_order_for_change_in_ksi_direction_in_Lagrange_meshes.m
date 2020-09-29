function New_node_order=find_node_order_for_change_in_ksi_direction_in_Lagrange_meshes(Nodes_per_elem_dir,ksi_direction_ID)

% created by Anastasia on 11-10-2019
%-------------------------------------------------------------------------------------------------------------------------
% Purpose: 
%-------------------------------------------------------------------------------------------------------------------------
% This script provides the new node order in elements when of the ksi
% directions is reversed. This was needed because the cubic hermite meshes
% from Pablo's script had non-right-hand ordered xi directions.
%-------------------------------------------------------------------------------------------------------------------------
% Notes on methodology:
%-------------------------------------------------------------------------------------------------------------------------
% I provide a vector with the new order that existing nodes should appear
% in the element matrix rows
%-------------------------------------------------------------------------------------------------------------------------
% Input:
%-------------------------------------------------------------------------------------------------------------------------
% Nodes_per_elem_dir: Lagrange interpolation order
% ksi_direction_ID = 1,2 or 3 (choose if you switch direction in ksi1, ksi2
% or ksi3

%-------------------------------------------------------------------------------------------------------------------------
% Output:
%-------------------------------------------------------------------------------------------------------------------------
% New_node_order: column vector with new node ordering in element e.g.
% 4,3,2,1,8,7,6,5,...,64,63,62,61 for a cubic lagrange interpolation with a
% change in ksi1 direction (Nodes_per_elem_dir=4, ksi_direction_ID=1)
%-------------------------------------------------------------------------------------------------------------------------
% SOS points:
%-------------------------------------------------------------------------------------------------------------------------
% meshes must be ordered ala cmgui (NO priority over corner nodes!)

New_node_order=zeros(Nodes_per_elem_dir^3,1);

            
%            % find new order
for n_ksi1=1:Nodes_per_elem_dir
    for n_ksi2=1:Nodes_per_elem_dir
        for n_ksi3=1:Nodes_per_elem_dir
            n_ksi1_new=n_ksi1; n_ksi2_new=n_ksi2; n_ksi3_new=n_ksi3;  
            node_ind=(n_ksi3-1)*Nodes_per_elem_dir^2+(n_ksi2-1)*Nodes_per_elem_dir+n_ksi1;
            if ksi_direction_ID==1
                n_ksi1_new=Nodes_per_elem_dir-n_ksi1+1;                
            elseif ksi_direction_ID==2
                n_ksi2_new=Nodes_per_elem_dir-n_ksi2+1;                
            elseif ksi_direction_ID==3
                n_ksi3_new=Nodes_per_elem_dir-n_ksi3+1;
            else
                disp('SOS! error in function find_node_order_for_change_in_ksi_direction_in_Lagrange_meshes - ksi_direction_ID should be 1, 2, or 3')
                return
            end
            
            New_node_order(node_ind)=(n_ksi3_new-1)*Nodes_per_elem_dir^2+(n_ksi2_new-1)*Nodes_per_elem_dir+n_ksi1_new; % to New_order prepei gia cubic lagrange na einai: [4,3,2,1,8,7,6,5,12,11,10,9...,64,63,62,61]
        end
    end
end
        

