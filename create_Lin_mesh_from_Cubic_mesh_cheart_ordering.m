function [Linear_mesh_Nodes,Linear_mesh_Elements,Linear_mesh_Boundaries]=create_Lin_mesh_from_Cubic_mesh_cheart_ordering(Cubic_mesh_nodes,Cubic_mesh_elements_ala_CHeart,Cubic_mesh_Boundaries,Nodes_per_elem_dir)

% % % for debugging: Cubic_mesh_nodes=Lagrange_mesh_Nodes_Ref; Cubic_mesh_elements_ala_CHeart=Lagrange_mesh_elements_ala_CHeart_Ref; Cubic_mesh_Boundaries=Lagrange_mesh_Boundaries;

% I think this can also work for quadratic meshes
% Elements should be ordered the Cheart way (priority over corner nodes)!!!


Linear_mesh_Elem_pre=Cubic_mesh_elements_ala_CHeart(:,1:8);
cubic_mesh_nodes_in_use=unique(Linear_mesh_Elem_pre);
Linear_mesh_Nodes=Cubic_mesh_nodes(cubic_mesh_nodes_in_use,:);
cub2lin_map=cubic_mesh_nodes_in_use; %index in cub2lin_map (1:length(cub2lin_map) is mapped to node id in cubic mesh
cub2lin_map_inverse=zeros(1,size(Cubic_mesh_nodes,1)); % row vector
cub2lin_map_inverse(cubic_mesh_nodes_in_use)=1:length(cubic_mesh_nodes_in_use); % this vector will have zeros in all places except at cubic_mesh_nodes_in_use and will point to the position of this cubic mesh node in the Linear_mesh_Nodes matrix
Linear_mesh_Elements=Linear_mesh_Elem_pre;
for n_nod=1:length(cub2lin_map) % for each node in Linear_mesh_Elem_pre that needs to be replaced by proper index in Linear_mesh_nodes using the cub2lin_map
    
    indices=find(Linear_mesh_Elements==cub2lin_map(n_nod));
    Linear_mesh_Elements(indices)=ones(size(indices))*n_nod;
end

% find Boundaries
Linear_mesh_Boundaries=zeros(size(Cubic_mesh_Boundaries,1),2^2+2); % preallocate (Nodes_per_elem_dir=2 for linear mesh)
Linear_mesh_Boundaries(:,1)=Cubic_mesh_Boundaries(:,1); % element ID stays same in Lin mesh
Linear_mesh_Boundaries(:,2^2+2)=Cubic_mesh_Boundaries(:,Nodes_per_elem_dir^2+2); % boundary ID stays same in Lin mesh
for n_b=1:size(Cubic_mesh_Boundaries,1)
    boundary_nodes_cubic=Cubic_mesh_Boundaries(n_b,2:(Nodes_per_elem_dir^2+1));
    boundary_nodes_cubic_in_use=boundary_nodes_cubic(ismember(boundary_nodes_cubic,cubic_mesh_nodes_in_use));% find which of these are corner nodes
    Linear_mesh_Boundaries(n_b,2:(2^2+1))=cub2lin_map_inverse(boundary_nodes_cubic_in_use);
end
 
% Linear_mesh_Neighbours=Cubic_mesh_Neighbours; - no need to generate extra
% variables.. Neighbours  matrix is the same in both cubic and linear
% meshes

% plot this on top of cubic mesh to check this is ok and then print -- 
