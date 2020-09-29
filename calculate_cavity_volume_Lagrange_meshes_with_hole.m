function [cavity_volume, cav_volume_mesh_Elements_cmgui_mult, cav_volume_mesh_Nodes_mult]=calculate_cavity_volume_Lagrange_meshes_with_hole(Lagrange_Elements_Cmgui,Lagrange_Nodes,Lagrange_Boundaries,Endo_patch_ID,Epi_patch_ID,Base_patch_ID,Apex_hole_patch_ID,Nodes_per_elem_dir,GP_per_elem_dir,fig_ind)

%% debug::
% Lagrange_Elements_Cmgui=Ref_mesh_Elements_ala_Cmgui; Lagrange_Nodes=Ref_mesh_Nodes;Lagrange_Boundaries=Ref_mesh_Boundaries;
% Endo_patch_ID=Endo_Boundary_ID;Epi_patch_ID=Epi_Boundary_ID;Base_patch_ID=Base_Boundary_ID; Apex_hole_patch_ID=ApexHole_Boundary_ID;

% created by Anastasia on 14-04-2020
%-------------------------------------------------------------------------------------------------------------------------
% Purpose: 
%-------------------------------------------------------------------------------------------------------------------------
% for p-V based cost function from meshes with hole
%-------------------------------------------------------------------------------------------------------------------------
% Notes on methodology:
%-------------------------------------------------------------------------------------------------------------------------
% A volumetric mesh will be created, consisting of collapsed elements.
% Endocardial nodes will be organised by layers and the central node at
% each layer will be given as the average of nodes in that layer. And then
% the collapsed element mesh is formed.
%-------------------------------------------------------------------------------------------------------------------------
% Input:
%-------------------------------------------------------------------------------------------------------------------------
% Endo_Boundary_ID: patch IDs belonging to endocardium, can be a vector (in case multiple endocardial patches are defined).
% fig_ind=[] if no plot of cavity mesh required, else specify figure number
%-------------------------------------------------------------------------------------------------------------------------
% Output:
%-------------------------------------------------------------------------------------------------------------------------
% cav_volume_mesh_Elements_cmgui: created cavity volume mesh element connectivity in Cmgui ordering
% cav_volume_mesh_Nodes: created cavity volume mesh nodes
% cavity_volume: the size of the cavity volume estimated this way
%-------------------------------------------------------------------------------------------------------------------------
% SOS points:
%-------------------------------------------------------------------------------------------------------------------------
% verified successfully against analytical solution for ellipsoid


[indicesEDGE_ksi1ksi2eq0, indicesEDGE_ksi1ksi2eq1, indicesEDGE_ksi1eq1ksi2eq0, indicesEDGE_ksi1eq0ksi2eq1,...
    indicesEDGE_ksi2ksi3eq0, indicesEDGE_ksi2ksi3eq1, indicesEDGE_ksi2eq1ksi3eq0, indicesEDGE_ksi2eq0ksi3eq1,...
    indicesEDGE_ksi3ksi1eq0, indicesEDGE_ksi3ksi1eq1, indicesEDGE_ksi3eq0ksi1eq1, indicesEDGE_ksi3eq1ksi1eq0]=internalNodeIDs_corresponding_natural_coord_EDGES_library(Nodes_per_elem_dir);


Endocardial_Elements=Lagrange_Boundaries(ismember(Lagrange_Boundaries(:,Nodes_per_elem_dir^2+2),Endo_patch_ID));
Endocardial_Nodes=unique(reshape(Lagrange_Boundaries(ismember(Lagrange_Boundaries(:,Nodes_per_elem_dir^2+2),Endo_patch_ID),2:(Nodes_per_elem_dir^2+1)),[],1));
Epicardial_Elements=Lagrange_Boundaries(ismember(Lagrange_Boundaries(:,Nodes_per_elem_dir^2+2),Epi_patch_ID));
Epicardial_Nodes=unique(reshape(Lagrange_Boundaries(ismember(Lagrange_Boundaries(:,Nodes_per_elem_dir^2+2),Epi_patch_ID),2:(Nodes_per_elem_dir^2+1)),[],1));
Basal_Elements=Lagrange_Boundaries(ismember(Lagrange_Boundaries(:,Nodes_per_elem_dir^2+2),Base_patch_ID));
Basal_Nodes=unique(reshape(Lagrange_Boundaries(ismember(Lagrange_Boundaries(:,Nodes_per_elem_dir^2+2),Base_patch_ID),2:(Nodes_per_elem_dir^2+1)),[],1));
Apical_hole_Elements=Lagrange_Boundaries(ismember(Lagrange_Boundaries(:,Nodes_per_elem_dir^2+2),Apex_hole_patch_ID));
Apical_hole_Nodes=unique(reshape(Lagrange_Boundaries(ismember(Lagrange_Boundaries(:,Nodes_per_elem_dir^2+2),Apex_hole_patch_ID),2:(Nodes_per_elem_dir^2+1)),[],1));

Endo_basal_Elements=Endocardial_Elements(ismember(Endocardial_Elements,Basal_Elements));
Endo_basal_Nodes=Endocardial_Nodes(ismember(Endocardial_Nodes,Basal_Nodes));
Endo_apical_Elements=Endocardial_Elements(ismember(Endocardial_Elements,Apical_hole_Elements));
Endo_apical_Nodes=Endocardial_Nodes(ismember(Endocardial_Nodes,Apical_hole_Nodes));
%% inverse node order along ksi3 to be used later for ensuring right hand system of local xi coordinates. (only chance of error is along ksi3 which is the constructed dimension because otherwise I maintain the endocardial surface element local ordering (and only build nodes layer per layer until central node is reached)
inv_ksi3=Nodes_per_elem_dir:-1:1;
for n_ksi3=1:Nodes_per_elem_dir
    n_ksi3_inv=inv_ksi3(n_ksi3);
    new_order(((n_ksi3-1)*Nodes_per_elem_dir^2+1):n_ksi3*Nodes_per_elem_dir^2)=((n_ksi3_inv-1)*Nodes_per_elem_dir^2+1):n_ksi3_inv*Nodes_per_elem_dir^2;
end
if mod(length(Endocardial_Elements),length(Endo_basal_Elements))==0
    surf_elem_nod_indices=1:Nodes_per_elem_dir^2;
    endo_elem_rows=length(Endocardial_Elements)/length(Endo_basal_Elements);
    cav_volume_mesh_Nodes_mult=zeros(length(Endocardial_Nodes)*Nodes_per_elem_dir,3); % col vector -- cavity volume mesh nodes where each central node is inserted multiple times
    cav_volume_mesh_Nodes_mult(1:length(Endocardial_Nodes),:)=Lagrange_Nodes(Endocardial_Nodes,:);
    
    % create nodes: first layer of nodes closer to the endocardial layer,
    % then moving layer by layer towards center
    Endocardial_nodes_rows=zeros((endo_elem_rows*(Nodes_per_elem_dir-1)+1),length(Endo_basal_Nodes)); % matrix with rows: node IDs at each endocardial node row (parallel circles of nodes starting from endo-basal nodes till endo-apical nodes.
    curr_Endo_nodes_row=[];
    curr_Endo_nodes_row=Endo_basal_Nodes.';% initialisation
    for n_row=1:size(Endocardial_nodes_rows,1)
        Endocardial_nodes_rows(n_row,:)=curr_Endo_nodes_row;
        Central_nodes(n_row,:)=mean(Lagrange_Nodes(curr_Endo_nodes_row,:)); % node at the center of each node layer
        
        %% move to next node row:
        nodes_to_check=[];
        nodes_to_check=Endocardial_Nodes(~ismember(Endocardial_Nodes,Endocardial_nodes_rows)).'; % row
        if isempty(nodes_to_check)==0 
            if n_row<size(Endocardial_nodes_rows,1)
                nodes_indices=[];
                nodes_indices=knnsearch(Lagrange_Nodes(nodes_to_check,:),Lagrange_Nodes(curr_Endo_nodes_row,:));
                curr_Endo_nodes_row=nodes_to_check(nodes_indices);
            else 
                disp('n_row is final but I can still find nodes to check - this shouldnot happen (nodes_to_check should be empty)')
            end
        else
            if n_row<size(Endocardial_nodes_rows,1)
                disp(['nodes_to_check appears empty at iteration: ',num2str(n_row), ' < ',num2str(length(Endocardial_nodes_rows)),'. This shouldnot happen!'])
            end
        end
    end  
    % find which central node each endocardial node corresponds to (by finding to which row in "Endocardial_nodes_rows" the endocardial node belongs to:
    Endo_node_row=zeros(length(Endocardial_Nodes),1); % col vector
    for n_nod=1:length(Endocardial_Nodes)
        Endo_node_row(n_nod,:)=find(sum(ismember(Endocardial_nodes_rows,Endocardial_Nodes(n_nod)),2)==1); % gia kathe node to row pou anikei sto Endocardial_nodes_rows        
        if isempty(find(sum(ismember(Endocardial_nodes_rows,Endocardial_Nodes(n_nod)),2)>1))==0 % an yparxei endocardial node pou anikei se panw apo 2 rows tote oli i thewrisi pou eftiaxa oti ta endoacardial nodes katatassontai se rows apo base pros apex katarriptetai kai prepei na skefteis pws tha to ftiakseis.. 
            
            warning('each endocardial node should correspond to exactly one row in Endocardial_nodes_rows --what is going on??')
            warning(' is your mesh with collapsed elements or capped apex where apical singularity exists? if so rethink this script cause it may break..')
            warning('exiting function calculate_cavity_volume_Lagrange_meshes_with_hole.m')
            return;
        end
    end
    %%   
    Dist_nodes=Lagrange_Nodes(Endocardial_Nodes,:)-Central_nodes(Endo_node_row,:);
    for n_nod=1:(Nodes_per_elem_dir-1) % create intermediate nodes between endocardial layer and including the central node
        cav_volume_mesh_Nodes_mult((n_nod*length(Endocardial_Nodes)+1):((n_nod+1)*length(Endocardial_Nodes)),:)=Lagrange_Nodes(Endocardial_Nodes,:)-Dist_nodes*(n_nod)/(Nodes_per_elem_dir-1); % moving from endocardium towards central nodes - first layer is taken by Endocardial nodes
    end
    % build cavity volume mesh elements
    % first you need a map between endocardial nodes in "Endocardial Nodes"
    % and in Ref_Nodes ID numbering:
    Inverse_map=zeros(max(Endocardial_Nodes)-min(Endocardial_Nodes)+1,1); % preallocate- the Inverse map will be an array whose rows correspond to the Node ID in Ref nodes and its entry is the ID of the node in "Endocardial_Nodes" (and has only size equal to max(Endocardial_Nodes) cause I only need to find this for the endocardial nodes so no reason to make the array bigger
    Inverse_map(Endocardial_Nodes-(min(Endocardial_Nodes)-1),:)=(1:1:length(Endocardial_Nodes)).'; % afairw to (min(Endocardial_Nodes)-1) gia na glitwsw xwro sto vector Inverse_Map (min exei tzampa theseis me midenika) alla isws mperdevei, tha to dw..
    
    for m_el=1:length(Endocardial_Elements)
        % map apo reference nodes IDs sto element, sta endocardial node IDs
        % (pou einai to layer of nodes sto cav_volume_mesh_Nodes_mult.
        % kai meta ftiakse kai ta pollaplasia tous IDs
        for n_nod=1:Nodes_per_elem_dir
            cav_volume_mesh_Elements_cmgui_mult(m_el,((n_nod-1)*Nodes_per_elem_dir^2+1):n_nod*Nodes_per_elem_dir^2)=(n_nod-1)*length(Endocardial_Nodes)+Inverse_map(Lagrange_Elements_Cmgui(Endocardial_Elements(m_el),ismember(Lagrange_Elements_Cmgui(Endocardial_Elements(m_el),:),Endocardial_Nodes))-(min(Endocardial_Nodes)-1),:);            
        end
        % check element is ordered in such a way that ksi1,ksi2,ksi3 local
        % directions form right hand system::
        ksi3_dir=cav_volume_mesh_Nodes_mult(cav_volume_mesh_Elements_cmgui_mult(m_el,indicesEDGE_ksi1ksi2eq0(end)),:)-cav_volume_mesh_Nodes_mult(cav_volume_mesh_Elements_cmgui_mult(m_el,indicesEDGE_ksi1ksi2eq0(1)),:);
        ksi1_dir=cav_volume_mesh_Nodes_mult(cav_volume_mesh_Elements_cmgui_mult(m_el,indicesEDGE_ksi2ksi3eq0(end)),:)-cav_volume_mesh_Nodes_mult(cav_volume_mesh_Elements_cmgui_mult(m_el,indicesEDGE_ksi2ksi3eq0(1)),:);
        ksi2_dir=cav_volume_mesh_Nodes_mult(cav_volume_mesh_Elements_cmgui_mult(m_el,indicesEDGE_ksi3ksi1eq0(end)),:)-cav_volume_mesh_Nodes_mult(cav_volume_mesh_Elements_cmgui_mult(m_el,indicesEDGE_ksi3ksi1eq0(1)),:);
        power_of_10_pre=norm(ksi3_dir)*norm(ksi2_dir)*norm(ksi1_dir);
        power_of_10 = round(abs(log10(power_of_10_pre)));
        power_of_10_TU=power_of_10+2; % so that you are taking issues of round-offs and binary representation into account
        if round(10^power_of_10_TU*dot(cross(ksi1_dir,ksi2_dir),ksi3_dir))<0
            cav_volume_mesh_Elements_cmgui_mult(m_el,:)=cav_volume_mesh_Elements_cmgui_mult(m_el,new_order);
            
        elseif round(10^power_of_10_TU*dot(cross(ksi1_dir,ksi2_dir),ksi3_dir))==0
            disp('possible error here these should not be exactly perpendicular --what''s going on?');
            disp('is it that your meshes are collapsed and one of ksi1_dir, ksi2_dir or ksi3_dir are =0?');s
            cav_volume_mesh_Elements_cmgui_mult(m_el,:)=cav_volume_mesh_Elements_cmgui_mult(m_el,new_order);
        end
    end
    if isempty(fig_ind)==0
        figure(fig_ind);
        plot_cmgui_ordered_mesh(cav_volume_mesh_Elements_cmgui_mult,cav_volume_mesh_Nodes_mult,Nodes_per_elem_dir,3)
    end
    
    
    %% volume calculation: 
    %% Cartesian coords basis vectors:
    e_x=[1,0,0]; e_y=[0,1,0]; e_z=[0,0,1];
    %% Gauss points in 3D element and shapefunctions and derivatives at GPs:
    [GPcoords, GPweights_mat]=Gauss_Points_In_3D_Element(GP_per_elem_dir,GP_per_elem_dir,GP_per_elem_dir);
    N_at_ksi_per_GP_mat=zeros(Nodes_per_elem_dir^3,size(GPcoords,1));
    thetaN_thetaksi_per_GP_mat=zeros(Nodes_per_elem_dir^3,3,size(GPcoords,1));
    for n_GP=1:size(GPcoords,1)
        N_at_ksi_per_GP_mat(:,n_GP)=Shapefunction_lagrange(GPcoords(n_GP,:),Nodes_per_elem_dir); % N_at_ksi_per_GP(:,n_GP): col vector --each col is Nat ksi for each GP
        thetaN_thetaksi_per_GP_mat(:,:,n_GP)=Shapefunction_derivative(GPcoords(n_GP,:),Nodes_per_elem_dir); % thetaN_thetaksi_per_GP(:,:,n_GP): [Nodes_per_elem_dir^3,3] size matrix (nodes in element, ksi directions)
    end
    
    %preallocations:
    Total_weights_per_GP_mat=zeros(size(cav_volume_mesh_Elements_cmgui_mult,1)*size(GPcoords,1),1);
    dV_per_GP_mat=zeros(size(cav_volume_mesh_Elements_cmgui_mult,1)*size(GPcoords,1),1);
    thetaX_thetaksi_per_GP_mat=zeros(size(cav_volume_mesh_Elements_cmgui_mult,1)*size(GPcoords,1),9);
    for m_el=1:size(cav_volume_mesh_Elements_cmgui_mult,1)
        X_at_Element_Nodes=cav_volume_mesh_Nodes_mult(cav_volume_mesh_Elements_cmgui_mult(m_el,:),:);
        for n_GP=1:size(GPcoords,1)
            GP_ID_total=(m_el-1)*size(GPcoords,1)+n_GP;
            thetaN_thetaksi=squeeze(thetaN_thetaksi_per_GP_mat(:,:,n_GP)); % [Nodes_per_elem_dir^3,3]
            N_atksi=N_at_ksi_per_GP_mat(:,n_GP); % col [Nodes_per_elem_dir^3,1]
            thetaX_thetaksi=X_at_Element_Nodes.'*thetaN_thetaksi;
            thetaX_thetaksi_per_GP_mat(GP_ID_total,:)=[thetaX_thetaksi(1,:),thetaX_thetaksi(2,:),thetaX_thetaksi(3,:)]; %thetaX/thetaksi1,thetaXthetaksi2,thetaXthetaksi3,thetaY/thetaksi1,thetaYthetaksi2,thetaYthetaksi3,thetaZ/thetaksi1,thetaZthetaksi2,thetaZthetaksi3]
            dV_per_GP_mat(GP_ID_total,:)=det(thetaX_thetaksi);
            Total_weights_per_GP_mat(GP_ID_total,:)=GPweights_mat(n_GP,1)*GPweights_mat(n_GP,2)*GPweights_mat(n_GP,3);
        end
    end
    cavity_volume=sum(dV_per_GP_mat.*Total_weights_per_GP_mat);
    %%
%     for n_ser=1:endo_elem_rows % I could use a while loop but prefer the for to avoid infinite loops in case of errors.
%         curr_Endo_nodes_row=Endo_basal_Nodes;
%         new_Central_node=mean(Lagrange_Nodes(curr_Endo_nodes_row,:));
%         Endo_elements_layer=Endo_basal_Elements;
%         % create new nodes from endo to cavity centre  per endocardial
%         % nodes layer:
%         
%         for m_el=1:length(Endo_elements_layer)
%             % create new layer of volume mesh elements
%             Endo_surf_elem=Lagrange_Elements_Cmgui(Endo_elements_layer(m_el),ismember(Lagrange_Elements_Cmgui(Endo_elements_layer(m_el),:),Endocardial_Nodes)); % ta node IDs tou endocardial element
%             endo_layer_nodes_in_Endo_surf_elem=Endo_surf_elem(ismember(Endo_surf_elem,curr_Endo_nodes_row));
%             
%             current_layer_indices=surf_elem_nod_indices(ismember(Endo_surf_elem,curr_Endo_nodes_row));% indices of nodes belonging to current endocardial node row in endocardial surface element
%             
%             cav_volume_mesh_Elements_cmgui((n_ser-1)*length(Endo_elements_layer)+m_el,current_layer_indices)=endo_layer_nodes_in_Endo_surf_elem;
%             
%             % find next layer of endocardial elements:
%             
%         end
%     end
else
    
    disp('this works on meshes where circumferential number of elements doesnot change');
end

% check volumetric mesh elements are ordered following a right hand system