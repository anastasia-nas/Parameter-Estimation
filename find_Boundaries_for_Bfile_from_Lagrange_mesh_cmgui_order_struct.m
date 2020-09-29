function Boundaries=find_Boundaries_for_Bfile_from_Lagrange_mesh_cmgui_order_struct(Lagrange_Nodes,Lagrange_Elements_CMGUI_ordering,Neighbours,Nodes_per_elem_dir,Endo_patch_ID,Epi_patch_ID,Base_patch_ID,Apex_hole_patch_ID)

% % % for debugging:: Lagrange_Nodes=Lagrange_mesh_Nodes_Ref; Lagrange_Elements_CMGUI_ordering=Lagrange_mesh_Elements_Cmgui; Neighbours=Lagrange_mesh_Neighbours; 


% modified 27/9/2019 from
% find_Boundaries_for_Bfile_from_Lagrange_mesh_cmgui_ordering_new.m to
% simplify the search because in quadratic meshes there were problems
% modified by Anastasia on 13-06-2019 from
% find_Boundaries_for_Bfile_from_Lagrange_mesh_cmgui_ordering.m (deleted it
% as that was not working)
%-------------------------------------------------------------------------------------------------------------------------
% Purpose: 
%-------------------------------------------------------------------------------------------------------------------------
% create Boundaries matrix (ready to be printed into *.B file)

%-------------------------------------------------------------------------------------------------------------------------
% Notes on methodology:
%-------------------------------------------------------------------------------------------------------------------------
% Like in the previous version I will find the nodes closest to the average
% node in the mesh, will find the element containing it and this has to be
% an endocardial element. That node might be an edge node, so find element
% faces containing this node. For each face find average of surf normals at
% GPs and identify which one points toward the midpoint (by finding the
% mean_surf_normal that has  the smallest angle (biggest absolute dot
% product) with the vector between central node and closest node. After you
% find the endocardial nodes of the first elements, find which of the
% neighbours contain one of the endocardial nodes. - but you still need to
% find which of these elements that will belong to edges are endo so you
% still need to go the way you're going--no need for finding the surface
% normal! Just find which of these patches don't have a neighbour
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
% Neighbours: matrix containing neighbours as given from function 'find_Neighbours_for_NBT_from_Lagrange_mesh_cmgui_ordering' (size=[elements_number,6]) - ordering [-ksi3,-ksi2,-ksi1,ksi1,ksi2,ksi3]
% Endo_patch_ID:
% Epi_patch_ID
% Base_patch_ID
% Apex_hole_patch_ID
% Apex_hole_patch_ID: leave empty if there is no apical hole!!! (Apex_hole_patch_ID=[])
%-------------------------------------------------------------------------------------------------------------------------
% Output:
%-------------------------------------------------------------------------------------------------------------------------
% Boundaries: matrix size=[elements per boundary,Nodes_per_elem_dir^2+2] -- (rows: element_ID, cols: [element_ID,nodes_in_patch_following_ksi_priority,Patch_ID])

%-------------------------------------------------------------------------------------------------------------------------
% SOS points:
%-------------------------------------------------------------------------------------------------------------------------
%! for structured meshes only
% NOT SUITABLE FOR COLLAPSED MESHES.
% Element matrix must be in CMGUI ordering (no priority over nodes) and
% elements must be of Lagrange type. NBT ordering is
% [-ksi3,-ksi2,-ksi1,ksi1,ksi2,ksi3]. 
% Automatically assumes that existence of 4th boundary corresponds to apical hole ! (see end of code- apical_hole_elements)
%===================================================================================================================================

% these are checked and are correct:
[indices_plane_ksi1eq0,indices_plane_ksi1eq1,indices_plane_ksi2eq0,indices_plane_ksi2eq1,indices_plane_ksi3eq0,indices_plane_ksi3eq1,...
    indices_plane_diag_EDGE1_ksi1ksi2eq0_EDGE2_ksi1ksi2eq1,indices_plane_diag_EDGE1_ksi1eq1ksi2eq0_EDGE2_ksi1eq0ksi2eq1,...
    indices_plane_diag_EDGE1_ksi2ksi3eq0_EDGE2_ksi2ksi3eq1,indices_plane_diag_EDGE1_ksi2eq1ksi3eq0_EDGE2_ksi2eq0ksi3eq1,...
    indices_plane_diag_EDGE1_ksi3ksi1eq0_EDGE2_ksi3ksi1eq1,indices_plane_diag_EDGE1_ksi3eq0ksi1eq1_EDGE2_ksi3eq1ksi1eq0]=internalNodeIDs_corresponding_natural_coord_planes_library(Nodes_per_elem_dir);

cell_boundary_for_neighbours={(1:6).',[indices_plane_ksi3eq0;indices_plane_ksi2eq0;indices_plane_ksi1eq0;indices_plane_ksi1eq1;indices_plane_ksi2eq1;indices_plane_ksi3eq1]}; % this assumes the order in NBT is [-ksi3,-ksi2,-ksi1,ksi1,ksi2,ksi3]. If you need to change it, change the ordering of indices_plane_ksi*eq* in the second matrix within the cell




% find elements that belong to boundaries (elements that have zeros in Neighbours matrix)
elements_in_boundaries=find(sum(ismember(Neighbours,0),2)>0);

% find elements that belong to 2 or more boundaries:
elements_in_two_or_more_patches=find(sum(ismember(Neighbours,0),2)>=2);
if length(elements_in_two_or_more_patches)==size(Lagrange_Elements_CMGUI_ordering,1) % if the mesh consists of 1 element in the transmural direction, then the elements that lie on the base and apex will belong to 3 boundaries and not just 1.
    max_common_boundaries_in_elem=max(sum(ismember(Neighbours,0),2));
    elements_at_endo_epi_base_hole_intesections=find(sum(ismember(Neighbours,0),2)==max_common_boundaries_in_elem);
else
    elements_at_endo_epi_base_hole_intesections=elements_in_two_or_more_patches;
end
% sanity check:
elements_in_over_2_patches=find(sum(ismember(Neighbours,0),2)>2); % this should be zero
if isempty(elements_in_over_2_patches)==0
    disp('elements should belong in up to two patches what is going on???????')
end

% find central point of mesh domain:
mesh_centre=mean(Lagrange_Nodes,1);


%% ========================================================ENDO========================================================================================================================================================================
%% Starting point - endo: find mesh node closest to element
[mesh_centre_mat,~]=meshgrid(mesh_centre,ones(size(Lagrange_Nodes,1),1));
[min_dist_from_centre,min_dist_from_centre_Node_ID]=min(sqrt(sum((Lagrange_Nodes-mesh_centre_mat).^2,2))); % central assumption is that this node will belong to endocardium
% [max_dist_from_centre,max_dist_from_centre_Node_ID]=max(sqrt(sum((Lagrange_Nodes-mesh_centre_mat).^2,2))); % this node might be epi-base or epi-apex, so not really useful for boundary identification


% Find element containing min distance from centre node (endocardium):
new_elements_contain_nodes=[];
new_elements_contain_nodes=find(sum(ismember(Lagrange_Elements_CMGUI_ordering,min_dist_from_centre_Node_ID),2)==1);
some_element_endocardial_nodes=min_dist_from_centre_Node_ID; % node closest to the cavity mesh (should be endocardial only - not lying in any other boundary- for an elongated LV where cavity mesh average lies in LV blood pool.)
Endocardial_Elements_Boundary_mat=[]; endocardial_nodes=[]; endocardial_elements=[];

% now write a while loop that ends when there are no new potential
% endocardial elements

while isempty(new_elements_contain_nodes)==0
    
    % use the first element in the list of new_elements_contain_nodes 
    disp(['finding endo nodes for element :', num2str(new_elements_contain_nodes(1))]);
    Endocardial_boundary_nodes_of_elem=find_endocardial_patch_in_element_from_some_endo_nodes(some_element_endocardial_nodes,new_elements_contain_nodes(1),Lagrange_Nodes,Lagrange_Elements_CMGUI_ordering,Nodes_per_elem_dir,Endo_patch_ID);
    


    Endocardial_Elements_Boundary_mat=[Endocardial_Elements_Boundary_mat;Endocardial_boundary_nodes_of_elem];
    endocardial_nodes=unique([endocardial_nodes,Endocardial_boundary_nodes_of_elem(2:Nodes_per_elem_dir^2+1)]);
    endocardial_elements=[endocardial_elements,new_elements_contain_nodes(1)];
    
    % find new element to identify its nodes:
    new_elements_contain_nodes=[]; new_elements_contain_nodes_pre=[];
    
    new_elements_contain_nodes_pre=find(sum(ismember(Lagrange_Elements_CMGUI_ordering,endocardial_nodes),2)>0); % find elements that contain at least an edge of endocardial nodes --could contain from Nodes_per_elem_dir(if only one neighbouring element has been identified) up to 4*Nodes_per_elem_dir endocardial nodes (nodes equivalent to 4 edges of endo patch if all 4 neighbouring endocardial elements have already been identified as endocardial elements 
    % of all the elements that contain endocardial nodes choose those that
    % haven't already been processed (aren't part of endocardial_elements
    % list yet)
    new_elements_contain_nodes=new_elements_contain_nodes_pre(~ismember(new_elements_contain_nodes_pre,endocardial_elements));
    if isempty(new_elements_contain_nodes)==0        
        some_element_endocardial_nodes=Lagrange_Elements_CMGUI_ordering(new_elements_contain_nodes(1),ismember(Lagrange_Elements_CMGUI_ordering(new_elements_contain_nodes(1),:),endocardial_nodes)); % find which nodes of element new_elements_contain_nodes(1) correspond to already identified endocadrial nodes and can thus serve as a starting point in identifying the endocardial nodes of the element 
    end
end




% % kane sort to boundary:
endocardial_elements=sort(endocardial_elements);
Endocardial_Elements_Boundary_mat=sortrows(Endocardial_Elements_Boundary_mat,1); % sort elements in ascending order

%% ========================================================EPI========================================================================================================================================================================
% Find Epicardial elements based on the structure of the mesh (going from
% -ksi to +ksi within element)
Epicardial_Elements_Boundary_mat=zeros(size(Endocardial_Elements_Boundary_mat));
epicardial_elements=zeros(size(endocardial_elements));
cell_boundary_new={(1:6).',[indices_plane_ksi1eq0;indices_plane_ksi2eq0;indices_plane_ksi3eq0;indices_plane_ksi1eq1;indices_plane_ksi2eq1;indices_plane_ksi3eq1]}; % I reorder this here because this makes use of modular arithmetic easier than using the ordering that was used for matching to the NBT file ordering (see cell_boundary ordering above)
Local_indices_in_elem=1:Nodes_per_elem_dir^3;
clear curr_elem;
for n_endo=1:size(Endocardial_Elements_Boundary_mat,1)
    curr_elem=Endocardial_Elements_Boundary_mat(n_endo,1);
    endocardial_nodes_in_elem=Endocardial_Elements_Boundary_mat(n_endo,2:Nodes_per_elem_dir^2+1);
    neigh_elem_towards_epi=curr_elem; % initialisation of while loop
    % psakse mexri na min vriskeis neighbour pou na exei ta antikrysta
    % nodes sto element patch pou einai closer to the endocardium
    while isempty(neigh_elem_towards_epi)==0
        curr_elem=neigh_elem_towards_epi;
        % find indices of endocardial elements in element:   
        Local_node_IDs_in_curr_patch=Local_indices_in_elem(ismember(Lagrange_Elements_CMGUI_ordering(curr_elem,:),endocardial_nodes_in_elem)); % Local nodes IDs: {1,64} of endocardial nodes in element
        for n_patch=1:6
            if sum(ismember(Local_node_IDs_in_curr_patch,cell_boundary_new{2}(n_patch,:)))==Nodes_per_elem_dir^2
                n_patch_oppos=mod(n_patch+3,6);
                if n_patch_oppos==0
                    n_patch_oppos=6;
                end
                Local_node_IDs_in_OPOSITE_patch=cell_boundary_new{2}(n_patch_oppos,:);
                endocardial_nodes_in_elem=Lagrange_Elements_CMGUI_ordering(curr_elem,Local_node_IDs_in_OPOSITE_patch); % these are the epicardial  nodes in current element that will serve as endocardial nodes in the element I'm searching for (neighbour of curr_elem towards epicardium)
            end
        end
        % find neighbouring element to curr_elem towards epi and name it curr_elem
        neigh_elem_towards_epi_cand=[];
        neigh_elem_towards_epi_cand=find(sum(ismember(Lagrange_Elements_CMGUI_ordering,endocardial_nodes_in_elem),2)==Nodes_per_elem_dir^2); % vaftizw tora curr_elem to geitoniko element pros epi
        neigh_elem_towards_epi=neigh_elem_towards_epi_cand(neigh_elem_towards_epi_cand~=curr_elem);
    end
    Epicardial_Elements_Boundary_mat(n_endo,:)=[curr_elem,Lagrange_Elements_CMGUI_ordering(curr_elem,ismember(Lagrange_Elements_CMGUI_ordering(curr_elem,:),endocardial_nodes_in_elem)),Epi_patch_ID]; % I use the "ismember" and don't just use endocardial_nodes_in_elem in cols 2:Nodes_per_elem_dir^2+1 to make sure the ordering is correct in epicardial element, as endocardial_nodes_in_elem ordering is based on the node ordering of the previous element to the epicardial one in the while loop                
    epicardial_elements(n_endo)=curr_elem;
end
epicardial_nodes=sort(unique(reshape(Epicardial_Elements_Boundary_mat(:,2:(Nodes_per_elem_dir^2+1) ),1,[])));
%% ========================================================BASE +APICAL HOLE========================================================================================================================================================================

% apo ta endocardial+epicardial afairese ta elements pou anikoun se 2
% patches. To apotelesma afairese to apo ta elements pou anikoun genika se
% patches kai exeis ta basal+apex_hole elements. Gia na ta ksexwriseis,
% vres ayta pou vriskontai panw apo to mesh centre (sti dieythynsi tou
% apicobasal direction) kai ayta pou einai stin anititheti dieythynsi
% anikoun sto apex_hole.
% Gia na vreis tous komvous tous pigaine sto Neighbours kai vres se poies
% dieythynseis den exoun neighbour. An den exoun mono se mia, tote ayti
% einai to basal boundary kai vres ta nodes, an einai 2 apokleise tin endo
% ws aytin pou periexei 16 theseis gemates endocardial nodes, enw i alli
% dieythynsi tha exei mono 4 (Nodes_per_elem_dir) along the endo-basal
% edge (dioti ta basal nodes den exoun akoma oristei afou den kserw poia elements einai sto basal boundary -ayto psaxnw)..

Endo_and_epi_elements=unique([endocardial_elements;epicardial_elements]);
non_base_apex_hole_elem=Endo_and_epi_elements(~ismember(Endo_and_epi_elements,elements_at_endo_epi_base_hole_intesections)); % boundary elements that are not part of base or apex hole boundaries
base_and_apex_hole_elements=elements_in_boundaries(~ismember(elements_in_boundaries,non_base_apex_hole_elem)); %base and apex hole boundary elements

% ksexwrise ta metaksy tous me vasi to Neighbours (ta basal elements kai ta
% apical hole elements den exoun neighbours metaksy tous)
boundary1_elements=base_and_apex_hole_elements(1); % initialisation
elem_list=boundary1_elements; % initialisation
while isempty(elem_list)==0
    new_elem=[];
    for nn=1:length(elem_list)
        curr_elem=elem_list(nn);
        neighbours_elem=Neighbours(curr_elem,:);
        new_elem=[new_elem,neighbours_elem(ismember(neighbours_elem,base_and_apex_hole_elements))];    % ayta tha ginoun to elem list sto epomeno iteration afou ginei ena ksekatharisma   
    end
    elem_list_pre=[];
    elem_list_pre=reshape(new_elem,1,[]);
    elem_list=[];
    elem_list=unique(elem_list_pre(~ismember(elem_list_pre,boundary1_elements))); % kratas mono ayta pou den exoun idi mpei sti lista boundary1_elements etsi wste na liksei kapoia stigmi to while loop, alliws that epaneksetazei synexeia ta idia kai ta idia elements  
    boundary1_elements=[boundary1_elements,elem_list];
end
boundary1_elements=sort(boundary1_elements);
boundary2_elements=base_and_apex_hole_elements(~ismember(base_and_apex_hole_elements,boundary1_elements)); % the remaining elements belong to the other boundary (apex or base)
boundary2_elements=sort(boundary2_elements);
if isempty(boundary2_elements)==0 % an yparxei kai tetarto boundary - ypothetw oti einai to apical hole
    % meta ta base einai ayta me to maximum radius kai to apical hole ayta me to minimum
    mean_boundary1=mean(Lagrange_Nodes(Lagrange_Elements_CMGUI_ordering(boundary1_elements,:),:));
    mean_boundary2=mean(Lagrange_Nodes(Lagrange_Elements_CMGUI_ordering(boundary2_elements,:),:));

    [mean_boundary1_mat,~]=meshgrid(mean_boundary1,ones(1,size(Lagrange_Nodes(Lagrange_Elements_CMGUI_ordering(boundary1_elements,:),:),1)));
    max_rad_boundary1=max(sqrt(sum((Lagrange_Nodes(Lagrange_Elements_CMGUI_ordering(boundary1_elements,:),:)-mean_boundary1_mat).^2,2)));


    mean_boundary2_mat=meshgrid(mean_boundary2,ones(1,size(Lagrange_Nodes(Lagrange_Elements_CMGUI_ordering(boundary2_elements,:),:),1)));
    max_rad_boundary2=max(sqrt(sum((Lagrange_Nodes(Lagrange_Elements_CMGUI_ordering(boundary2_elements,:),:)-mean_boundary2_mat).^2,2)));

    if max_rad_boundary1>max_rad_boundary2
        basal_elements=boundary1_elements;
        apical_hole_elements=boundary2_elements;
    elseif max_rad_boundary1<max_rad_boundary2
        basal_elements=boundary2_elements;
        apical_hole_elements=boundary1_elements;
    end
else
    basal_elements=boundary1_elements;
    apical_hole_elements=[];
end



%% ========================================================BASE========================================================================================================================================================================
% now build the boundary patch matrix for basal elements
Basal_Elements_Boundary_mat=zeros(length(basal_elements),Nodes_per_elem_dir^2+2); % preallocate
for n_el=1:length(basal_elements)
    candindate_face_indices=[]; 
    candindate_face_indices=find(Neighbours(basal_elements(n_el),:)==0); %element faces that lie on boundary patches
    for n_c=1:length(candindate_face_indices)
        cand_nodes_IDs=Lagrange_Elements_CMGUI_ordering(basal_elements(n_el),cell_boundary_for_neighbours{2}(candindate_face_indices(n_c),:)); %cell_boundary_for_neighbours{2}(candindate_face_indices(n_c),:): einai ta local nodes IDs twn nodes pou anikoun sto en logw patch
        % to face pou anoikei sto basal boundary tha einai ayto pou exei ligotera apo Nodes_per_elem_dir^2 nodes pou na anikoun sto
        % endocardial i epicardial nodes (mporei omws na exei Nodes_per_elem_dir nodes (osa diladi anoikoun se ena element edge) pou na anikoun sta endocardial_nodes i
        % epicardial_nodes ean to basal element anikei kai sto endo i epi ektos apo base)
        if sum(ismember(cand_nodes_IDs,endocardial_nodes))<=Nodes_per_elem_dir && sum(ismember(cand_nodes_IDs,epicardial_nodes))<=Nodes_per_elem_dir
            Basal_Elements_Boundary_mat(n_el,:)=[basal_elements(n_el),cand_nodes_IDs,Base_patch_ID];
        end
    end
end
%% ========================================================APICAL HOLE========================================================================================================================================================================
if isempty(Apex_hole_patch_ID)==0
    Apical_hole_Elements_Boundary_mat=zeros(length(apical_hole_elements),Nodes_per_elem_dir^2+2); % preallocate
    for n_el=1:length(apical_hole_elements)
        candindate_face_indices=[]; 
        candindate_face_indices=find(Neighbours(apical_hole_elements(n_el),:)==0); %element faces that lie on boundary patches
        for n_c=1:length(candindate_face_indices)
            cand_nodes_IDs=Lagrange_Elements_CMGUI_ordering(apical_hole_elements(n_el),cell_boundary_for_neighbours{2}(candindate_face_indices(n_c),:)); %cell_boundary_for_neighbours{2}(candindate_face_indices(n_c),:): einai ta local nodes IDs twn nodes pou anikoun sto en logw patch
            % to face pou anoikei sto apical hole boundary tha einai ayto pou exei ligotera apo Nodes_per_elem_dir^2 nodes pou na anikoun sto
            % endocardial i epicardial nodes (mporei omws na exei Nodes_per_elem_dir nodes (osa diladi anoikoun se ena element edge) pou na anikoun sta endocardial_nodes i
            % epicardial_nodes ean to basal element anikei kai sto endo i epi ektos apo base)
            if sum(ismember(cand_nodes_IDs,endocardial_nodes))<=Nodes_per_elem_dir && sum(ismember(cand_nodes_IDs,epicardial_nodes))<=Nodes_per_elem_dir
                Apical_hole_Elements_Boundary_mat(n_el,:)=[apical_hole_elements(n_el),cand_nodes_IDs,Apex_hole_patch_ID];
            end
        end
    end
else
    Apical_hole_Elements_Boundary_mat=[];
end

%% Bring them all together in one "Boundaries" matrix to be used for Cheart .B file:
Boundaries=[Endocardial_Elements_Boundary_mat;Epicardial_Elements_Boundary_mat;Basal_Elements_Boundary_mat;Apical_hole_Elements_Boundary_mat];
    
