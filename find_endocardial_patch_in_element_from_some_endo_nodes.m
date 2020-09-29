function Endocardial_boundary_nodes_of_elem=find_endocardial_patch_in_element_from_some_endo_nodes(some_element_endocardial_nodes,element_ID,mesh_Nodes,mesh_Elements_cmgui_order,Nodes_per_elem_dir,endo_patch_ID)


% created by Anastasia on 17-10-2019
%-------------------------------------------------------------------------------------------------------------------------
% Purpose: 
%-------------------------------------------------------------------------------------------------------------------------
% this function identifies the endocardial nodes of an element given a set of endocardial nodes (>=1),a specific element ID that contains these nodes and the 
% matrices of mesh nodes and element connectivity.
% This function is part of a script that finds the boundaries from
% structured meshes.
%-------------------------------------------------------------------------------------------------------------------------
% Notes on methodology:
%-------------------------------------------------------------------------------------------------------------------------
% Will not work on collapsed elements.
% Find patches containing these endocardial nodes. Depending on the size of
% "some_element_endocardial_nodes", there can be up to 3 patches containing
% these nodes (for example if only one endocardial node is provided then
% this node can belong to up to 3 patches if it sits at an element corner).
% To find which is the patch with the endocardial nodes find the one that
% has the biggest absolute (I use the absolute so that any definition of the
% surface normal positive direction is irrelevant) dot product between the
% surface normal of the patch and the vector joining the average
% endocardial node to the centre of the mesh (approximately equal to centre
% of blood cavity) and that distance vector serves as a surrogate for the
% ground truth surface normal of the endocardial patch.
% But check among these candidate patches, which is the one that has no
% neighbour.
%-------------------------------------------------------------------------------------------------------------------------
% Input:
%-------------------------------------------------------------------------------------------------------------------------
% some_element_endocardial_nodes: nodes IDs that lie on the endocarial patch of the element (usually these are the neighbouring endocardial nodes with an
%                                   already identified endocardial element, or just the node closest to the cavity centre, and which helps the initialisation
%                                   of finding the endocardial element)
% element_ID: element ID in mesh
% mesh_Nodes: Nodes coordinates in mesh
% mesh_Elements_cmgui_order: Element connectivity in mesh
% Nodes_per_elem_dir: mesh interpolation order (only works with Lagrange interpolation)
% endo_patch_ID: the boundary ID assigned to endocardium in the boundaries
% matrix.
%-------------------------------------------------------------------------------------------------------------------------
% Output:
%-------------------------------------------------------------------------------------------------------------------------
% Endocardial_boundary_nodes_of_elem: row vector (size: [1,Nodes_per_elem_dir^2]) containing Nodes_per_elem_dir^2 nodes exactly as they would be ordered in the
% Boundaries matrix (ordered first along ksi1 and then along ksi2 in terms of the ksi coords of the patch surface.
%-------------------------------------------------------------------------------------------------------------------------
% SOS points:
%-------------------------------------------------------------------------------------------------------------------------
% Not for collapsed elements!
% SOS!!: applies only to LV as it assumes that the endocardial nodes in the
% first element will be at endo only (not endo-base or endo-hole) so that
% the only patch containing the node and not being connected to a
% neighbouring element is endo. Otherwise it could be base or hole.
% 
% So this will break in meshes where the LV looks like an oblate spheroid
% where shortest distance from cavity mesh is not strictly endo (endo
% midwall) but can be endo-base or endo-apical hole.

%% cavity centre - will only be used if the surface normal will be used to find which boundary patch is endo
cavity_centre=mean(mesh_Nodes,1);

%% find which element patch the nodes can belong to.
[indices_plane_ksi1eq0,indices_plane_ksi1eq1,indices_plane_ksi2eq0,indices_plane_ksi2eq1,indices_plane_ksi3eq0,indices_plane_ksi3eq1,...
    indices_plane_diag_EDGE1_ksi1ksi2eq0_EDGE2_ksi1ksi2eq1,indices_plane_diag_EDGE1_ksi1eq1ksi2eq0_EDGE2_ksi1eq0ksi2eq1,...
    indices_plane_diag_EDGE1_ksi2ksi3eq0_EDGE2_ksi2ksi3eq1,indices_plane_diag_EDGE1_ksi2eq1ksi3eq0_EDGE2_ksi2eq0ksi3eq1,...
    indices_plane_diag_EDGE1_ksi3ksi1eq0_EDGE2_ksi3ksi1eq1,indices_plane_diag_EDGE1_ksi3eq0ksi1eq1_EDGE2_ksi3eq1ksi1eq0]=internalNodeIDs_corresponding_natural_coord_planes_library(Nodes_per_elem_dir);

cell_boundary_for_neighbours={(1:6).',[indices_plane_ksi3eq0;indices_plane_ksi2eq0;indices_plane_ksi1eq0;indices_plane_ksi1eq1;indices_plane_ksi2eq1;indices_plane_ksi3eq1]}; % this assumes the order in NBT is [-ksi3,-ksi2,-ksi1,ksi1,ksi2,ksi3]. If you need to change it, change the ordering of indices_plane_ksi*eq* in the second matrix within the cell

unit_indices_in_elem=ismember(mesh_Elements_cmgui_order(element_ID,:),some_element_endocardial_nodes);
endo_nodes_indices_in_elem=find(unit_indices_in_elem==1);

identify_cand_patch_pre=ismember(cell_boundary_for_neighbours{2},endo_nodes_indices_in_elem);
identify_cand_patch=find(sum(identify_cand_patch_pre,2)==length(some_element_endocardial_nodes)); % ta indices twn patches sto cell_boundary_for_neighbours pou periexoun olous tous endocardial nodes

% find out of these patches which one doesn't share its notes with a
% neighbouring element.
local_endo_patch_ID=[];
for n_p=1:length(identify_cand_patch)
    cand_neigh_elem_pre=[]; cand_neigh_elem=[];
    cand_neigh_elem_pre_ones=ismember(mesh_Elements_cmgui_order,mesh_Elements_cmgui_order(element_ID,cell_boundary_for_neighbours{2}(identify_cand_patch(n_p),:))); % poia elements exoun koina ola ta nodes tou candidate patch
    cand_neigh_elem_pre=find(sum(cand_neigh_elem_pre_ones,2)==Nodes_per_elem_dir^2);
    cand_neigh_elem=cand_neigh_elem_pre(find(cand_neigh_elem_pre~=element_ID));
    if isempty(cand_neigh_elem)==1
        local_endo_patch_ID=[local_endo_patch_ID,n_p];
    end
end

if isempty(local_endo_patch_ID)==1 
    disp('IN FUNCTION find_endocardial_patch_in_element_from_some_endo_nodes.m');
    disp('wtf at least one patch containing the some_element_endocardial_nodes, and not sharing them with a neighbouring element should have been identified!')
elseif length(local_endo_patch_ID)>1 % this is is highly improbable, as I am using as a criterion that the patch should include all endocardial nodes, so we should always be between an endo face and a face that is neighbouring with the previously identified endo element
    
    disp (' this should not be happening -- write code using collinearity to distance from cavity mesh as an additional criterion, although this should not work as you move away from the midwall towards apex or base...');
    
else
    Endocardial_boundary_nodes_of_elem=[element_ID,mesh_Elements_cmgui_order(element_ID,cell_boundary_for_neighbours{2}(identify_cand_patch(local_endo_patch_ID),:)),endo_patch_ID];
end


%
