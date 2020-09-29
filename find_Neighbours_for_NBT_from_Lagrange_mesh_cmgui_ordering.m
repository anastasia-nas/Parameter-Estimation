function Neighbours=find_Neighbours_for_NBT_from_Lagrange_mesh_cmgui_ordering(Lagrange_Elements_CMGUI_ordering,Nodes_per_elem_dir)

% created by Anastasia on 03-06-2019
%-------------------------------------------------------------------------------------------------------------------------
% Purpose: 
%-------------------------------------------------------------------------------------------------------------------------
% create Neighbours matrix (ready to be printed into *.NBT file)

%-------------------------------------------------------------------------------------------------------------------------
% Notes on methodology:
%-------------------------------------------------------------------------------------------------------------------------
% neighbours file ordering: 
% [-ksi3,-ksi2,-ksi1,ksi1,ksi2,ksi3] -- according to the NBT file in the
% capped ellipsoid meshes.

%-------------------------------------------------------------------------------------------------------------------------
% Input:
%-------------------------------------------------------------------------------------------------------------------------
% Lagrange_Nodes: Nodes matrix (rows: node_ID, cols, nodal coordinate value) - size=[nodes_number,3] - corresponds 
%                 to Lagrange meshes of order indicated by the
%                 Nodes_per_elem_dir input (for cubic =4, for linear =2
%                 etc) -- I don't need this! I only need the element
%                 connectivity matrix (Lagrange_Elements)
% Lagrange_Elements: Elements matrix in CMGUI ordering size=[element_number, Nodes_per_elem_dir^3]- no priority over corner nodes 
%                    (rows: element_ID, cols: local node ID corresponding to global node id)

%-------------------------------------------------------------------------------------------------------------------------
% Output:
%-------------------------------------------------------------------------------------------------------------------------
% Neighbours: matrix size=[elements_number,6] -- (rows: element_ID, cols: neighbouring element in this order: [-ksi3,-ksi2,-ksi1,ksi1,ksi2,ksi3]

%-------------------------------------------------------------------------------------------------------------------------
% SOS points:
%-------------------------------------------------------------------------------------------------------------------------
% Nodes matrix must be in CMGUI ordering (no priority over nodes) and
% elements must be of Lagrange type. NBT ordering is
% [-ksi3,-ksi2,-ksi1,ksi1,ksi2,ksi3].
%===================================================================================================================================
Neighbours=zeros(size(Lagrange_Elements_CMGUI_ordering,1),6);

[indices_plane_ksi1eq0,indices_plane_ksi1eq1,indices_plane_ksi2eq0,indices_plane_ksi2eq1,indices_plane_ksi3eq0,indices_plane_ksi3eq1,...
    indices_plane_diag_EDGE1_ksi1ksi2eq0_EDGE2_ksi1ksi2eq1,indices_plane_diag_EDGE1_ksi1eq1ksi2eq0_EDGE2_ksi1eq0ksi2eq1,...
    indices_plane_diag_EDGE1_ksi2ksi3eq0_EDGE2_ksi2ksi3eq1,indices_plane_diag_EDGE1_ksi2eq1ksi3eq0_EDGE2_ksi2eq0ksi3eq1,...
    indices_plane_diag_EDGE1_ksi3ksi1eq0_EDGE2_ksi3ksi1eq1,indices_plane_diag_EDGE1_ksi3eq0ksi1eq1_EDGE2_ksi3eq1ksi1eq0]=internalNodeIDs_corresponding_natural_coord_planes_library(Nodes_per_elem_dir);


cell_boundary={[1:6].',[indices_plane_ksi3eq0;indices_plane_ksi2eq0;indices_plane_ksi1eq0;indices_plane_ksi1eq1;indices_plane_ksi2eq1;indices_plane_ksi3eq1]}; % this assumes the order in NBT is [-ksi3,-ksi2,-ksi1,ksi1,ksi2,ksi3]. If you need to change it, change the ordering of indices_plane_ksi*eq* in the second matrix within the cell
for m_el=1:size(Lagrange_Elements_CMGUI_ordering,1)
    for n_bound=1:length(cell_boundary{1})
        element_boundary_Nodes=Lagrange_Elements_CMGUI_ordering(m_el,cell_boundary{2}(n_bound,:));
%     ksi1_negative_Nodes=Lagrange_Elements(m_el,indices_plane_ksi1eq0);
%     ksi1_positive_Nodes=Lagrange_Elements(m_el,indices_plane_ksi1eq1);
%     ksi2_negative_Nodes=Lagrange_Elements(m_el,indices_plane_ksi2eq0);
%     ksi2_positive_Nodes=Lagrange_Elements(m_el,indices_plane_ksi2eq1);
%     ksi3_negative_Nodes=Lagrange_Elements(m_el,indices_plane_ksi3eq0);
%     ksi3_positive_Nodes=Lagrange_Elements(m_el,indices_plane_ksi3eq1);
    
        element_indices=ismember(Lagrange_Elements_CMGUI_ordering,element_boundary_Nodes); 
        neigh_candidates_at_boundary=find(sum(element_indices,2)==length(element_boundary_Nodes)); % there should only be up to 2 elements containing the boundary nodes: the current element and its neighbour
        neigh_ind=find(neigh_candidates_at_boundary~=m_el);
        if isempty(neigh_ind)==0 % if there is another element containing these nodes other than the current one assign its id at the right position; else do nothing (by preallocation this position is already zero
             Neighbours(m_el,cell_boundary{1}(n_bound))=neigh_candidates_at_boundary(neigh_ind);         
        end
    end       
        
    
    
end

