function [Wext_endo_per_DF, virt_Wext_endo_per_DF,virt_Wext_epi_per_DF,virt_Wext_epi_over_press_per_DF]=calc_Wext_endo_and_virtWext_endo_per_frame_for_CF_verbose_epi(diastolic_frames_list,positive_pressures_at_diastolic_frames,Elements_cmgui,Boundaries,endo_patch_ID,epi_patch_ID,Def_Nodes_per_Frame_struct,virt_disp_Nodes_per_Frame_struct,GP_per_elem_dir,Nodes_per_elem_dir,flag_calc_virt_work_epi)
% % % % % % for debug:
% % endo_patch_ID=Endo_patch_ID; epi_patch_ID=Epi_patch_ID; Boundaries=Ref_mesh_Boundaries; Elements_cmgui=Ref_mesh_Elements_ala_Cmgui; positive_pressures_at_diastolic_frames=diast_pressures_Pa_TU; diastolic_frames_list=Diastolic_Frames_Indices_TU;
% created by Anastasia on 24-02-2020
%-------------------------------------------------------------------------------------------------------------------------
% Purpose: 
%-------------------------------------------------------------------------------------------------------------------------

%-------------------------------------------------------------------------------------------------------------------------
% Notes on methodology:
%-------------------------------------------------------------------------------------------------------------------------
% Wext is only Wext at endo
%-------------------------------------------------------------------------------------------------------------------------
% Input:
%-------------------------------------------------------------------------------------------------------------------------
% diastolic_frames_list: this needs to be all frames from reference frame on e.g. [27,28,1,2,3] - EXCLUDING the reference frame
% positive_pressures_at_diastolic_frames: POSITIVE for INFLATION! corresponding pressure to each frame: e.g. [1500,1600,1700,1900,2000]-excluding pressure 
%                                 at reference frame which is assumed to be zero!
% Ref_Nodes: mesh nodes at reference frame
% Elements_cmgui:elements with with cmgui ordering (no priority at corner nodes)
% Def_Nodes_per_Frame_struct: a structure containing the deformed mesh at each frame (e.g. Def_Nodes_per_Frame_struct.Frame_25,Def_Nodes_per_Frame_struct.Frame_27, Def_Nodes_per_Frame_struct.Frame_1 etc..
% virt_disp_Nodes_per_Frame_struct: virtual field used for each deformed state, again in a structure: virt_disp_Nodes_per_Frame_struct.Frame_25, virt_disp_Nodes_per_Frame_struct.Frame_27, virt_disp_Nodes_per_Frame_struct.Frame_1 etc
%                                  if =[] then all virtual field related variables will be left empty
% GP_per_elem_dir: number of GPs per element dimension
% Nodes_per_elem_dir: Lagrange mesh nodes number per element dimension indicating the interpolation order (=4 for Cubic).


%-------------------------------------------------------------------------------------------------------------------------
% Output:
%-------------------------------------------------------------------------------------------------------------------------
% Wext_endo_per_DF: vector of Wext at endo per diastolic frame - column vector - size(length(diastolic_frames_list),1)
% virt_Wext_endo_per_DF: vector of virtual Wext at endo per diastolic frame - column vector - size(length(diastolic_frames_list),1)
%-------------------------------------------------------------------------------------------------------------------------
% SOS points:
%-------------------------------------------------------------------------------------------------------------------------

%% ==============================================START==============================================================
local_nodes_IDs_in_elem=1:Nodes_per_elem_dir^3;
%% Gauss points in 3D element and shapefunctions and derivatives at GPs:
[surfGPcoords, surfGPweights]=Gauss_Points_In_2D_Element(GP_per_elem_dir,GP_per_elem_dir);
Nsurf_at_ksi_per_GP_mat=zeros(Nodes_per_elem_dir^2,size(surfGPcoords,1));
thetaNsurf_thetaksi_per_GP_mat=zeros(Nodes_per_elem_dir^2,2,size(surfGPcoords,1));
for n_GP=1:size(surfGPcoords,1)
    Nsurf_at_ksi_per_GP_mat(:,n_GP)=Shapefunction_lagrange_Surf(surfGPcoords(n_GP,:),Nodes_per_elem_dir); % N_at_ksi_per_GP(:,n_GP): col vector --each col is Nat ksi for each GP
    thetaNsurf_thetaksi_per_GP_mat(:,:,n_GP)=Surface_Shapefunction_derivative(surfGPcoords(n_GP,:),Nodes_per_elem_dir);% thetaN_thetaksi_per_GP(:,:,n_GP): [Nodes_per_elem_dir^3,3] size matrix (nodes in element, ksi directions)    
end
surfGPweights_total=surfGPweights(:,1).*surfGPweights(:,2);

%% find endocardial elements and endocardial nodes:
endocardial_elements=Boundaries(ismember(Boundaries(:,Nodes_per_elem_dir^2+2),endo_patch_ID),1);
endocardial_nodes=sort(unique(reshape(Boundaries(ismember(Boundaries(:,Nodes_per_elem_dir^2+2),endo_patch_ID),2:Nodes_per_elem_dir^2+1),1,[])));

epicardial_elements=Boundaries(ismember(Boundaries(:,Nodes_per_elem_dir^2+2),epi_patch_ID),1);
epicardial_nodes=sort(unique(reshape(Boundaries(ismember(Boundaries(:,Nodes_per_elem_dir^2+2),epi_patch_ID),2:Nodes_per_elem_dir^2+1),1,[])));

%% calculate Wext at each frame
% preallocations:
virt_Wext_epi_per_DF=[]; virt_Wext_epi_over_press_per_DF=[];
p_da_per_GP_and_frame=zeros(length(endocardial_elements)*size(surfGPcoords,1),3,length(diastolic_frames_list));
x_per_GP_and_frame=zeros(length(endocardial_elements)*size(surfGPcoords,1),3,length(diastolic_frames_list));
p_da_virt_u_w_per_GP_and_frame=zeros(length(endocardial_elements)*size(surfGPcoords,1),length(diastolic_frames_list));
p_da_virt_u_w_per_GP_and_frame_epi=zeros(length(epicardial_elements)*size(surfGPcoords,1),length(diastolic_frames_list));
da_virt_u_w_per_GP_and_frame_epi=zeros(length(epicardial_elements)*size(surfGPcoords,1),length(diastolic_frames_list));
Global_GP_weights=zeros(length(endocardial_elements)*size(surfGPcoords,1),1);
Wext_endo_incr_per_DF=zeros(length(diastolic_frames_list),1); % external work per increment: (\sum_GP dot([(p*da)_n+(p*da)_n-1],[x_n-x_n-1])*W_GP )
Wext_endo_per_DF=zeros(length(diastolic_frames_list),1);
for n_DF=1:length(diastolic_frames_list)
    Def_Nodes=Def_Nodes_per_Frame_struct.(['Frame_',num2str(diastolic_frames_list(n_DF))]);
    if isempty(virt_disp_Nodes_per_Frame_struct)==0
        Virtual_Disp_Nodes=virt_disp_Nodes_per_Frame_struct.(['Frame_',num2str(diastolic_frames_list(n_DF))]);
    else
        Virtual_Disp_Nodes=[];
    end
    press=-positive_pressures_at_diastolic_frames(n_DF); % I reverse sign here because I want positive pressure in input to correspond to inflation (opposite sign to da vector!)
    for m_el=1:length(endocardial_elements)
        Def_Nodes_at_surfElem=Def_Nodes(Elements_cmgui(endocardial_elements(m_el),ismember(Elements_cmgui(endocardial_elements(m_el),:),endocardial_nodes)),:);        
        Local_node_IDs_in_curr_patch=local_nodes_IDs_in_elem(ismember(Elements_cmgui(endocardial_elements(m_el),:),endocardial_nodes));
        Local_node_IDs_in_OPPOSITE_patch=find_local_nodeIDs_of_opposite_patch_given_curr_elem_patch_IDs(Local_node_IDs_in_curr_patch,Nodes_per_elem_dir);

        for n_GP=1:size(surfGPcoords,1)
            GP_ID_global=(m_el-1)*size(surfGPcoords,1)+n_GP;
            thetaNsurf_thetaksi=squeeze(thetaNsurf_thetaksi_per_GP_mat(:,:,n_GP));
            thetax_thetaksisurf=Def_Nodes_at_surfElem.'*thetaNsurf_thetaksi;
            da_vec=cross(thetax_thetaksisurf(:,1),thetax_thetaksisurf(:,2)).';
            x_at_GP=Nsurf_at_ksi_per_GP_mat(:,n_GP).'*Def_Nodes_at_surfElem; % position vector at GP - I use it only as correction for da_vec
            
            
            [~,min_ind_nod]=min(sum((Def_Nodes(Elements_cmgui(endocardial_elements(m_el),Local_node_IDs_in_OPPOSITE_patch),:)-[x_at_GP(1)*ones(length(Local_node_IDs_in_OPPOSITE_patch),1),x_at_GP(2)*ones(length(Local_node_IDs_in_OPPOSITE_patch),1),x_at_GP(3)*ones(length(Local_node_IDs_in_OPPOSITE_patch),1)]).^2,2));
            nearest_nod_opposite_patch=Def_Nodes(Elements_cmgui(endocardial_elements(m_el),Local_node_IDs_in_OPPOSITE_patch(min_ind_nod)),:);
            
            if dot(da_vec,x_at_GP-nearest_nod_opposite_patch)<0 % if da points opposite to the vector pointing from the nearest node at the opposite (epicardial)patch to the endocardial surface GP
                %the above control will only be compromised if the element
                %is very long and thin, so that the two vectors may become
                %perpendicular and the dot product zero, so that numerical
                %roundoffs may make it a very small negative number..-
                %that's why I use the round (10^5*..) extra level of
                %security
                da_vec=-da_vec; % then change direction so da points towards blood cavity
            end
            p_da_per_GP_and_frame(GP_ID_global,:,n_DF)=press*da_vec;
            if isempty(Virtual_Disp_Nodes)==0
                Virtual_Disp_Nodes_at_surfElem=Virtual_Disp_Nodes(Elements_cmgui(endocardial_elements(m_el),ismember(Elements_cmgui(endocardial_elements(m_el),:),endocardial_nodes)),:);
                virt_u=Nsurf_at_ksi_per_GP_mat(:,n_GP).'*Virtual_Disp_Nodes_at_surfElem;
                p_da_virt_u_w_per_GP_and_frame(GP_ID_global,n_DF)=press*dot(virt_u,da_vec)*surfGPweights_total(n_GP);%p*da \cdot delta_u*w_GP - just sum this and you've got \deltaWext
            end
            x_per_GP_and_frame(GP_ID_global,:,n_DF)=Nsurf_at_ksi_per_GP_mat(:,n_GP).'*Def_Nodes_at_surfElem;
            Global_GP_weights(GP_ID_global,:)=surfGPweights_total(n_GP);            
        end                
    end
    %% estimate the virtual work on epicardium to verify the intervention works and to what extent:
    for m_el=1:length(epicardial_elements)
        Def_Nodes_at_surfElem=Def_Nodes(Elements_cmgui(epicardial_elements(m_el),ismember(Elements_cmgui(epicardial_elements(m_el),:),epicardial_nodes)),:);
        if isempty(Virtual_Disp_Nodes)==0
            Virtual_Disp_Nodes_at_surfElem=Virtual_Disp_Nodes(Elements_cmgui(epicardial_elements(m_el),ismember(Elements_cmgui(epicardial_elements(m_el),:),epicardial_nodes)),:);
        end
        Local_node_IDs_in_curr_patch=local_nodes_IDs_in_elem(ismember(Elements_cmgui(epicardial_elements(m_el),:),epicardial_nodes));
        Local_node_IDs_in_OPPOSITE_patch=find_local_nodeIDs_of_opposite_patch_given_curr_elem_patch_IDs(Local_node_IDs_in_curr_patch,Nodes_per_elem_dir);

        for n_GP=1:size(surfGPcoords,1)
            GP_ID_global=(m_el-1)*size(surfGPcoords,1)+n_GP;
            thetaNsurf_thetaksi=squeeze(thetaNsurf_thetaksi_per_GP_mat(:,:,n_GP));
            thetax_thetaksisurf=Def_Nodes_at_surfElem.'*thetaNsurf_thetaksi;
            da_vec=cross(thetax_thetaksisurf(:,1),thetax_thetaksisurf(:,2)).';
            x_at_GP=Nsurf_at_ksi_per_GP_mat(:,n_GP).'*Def_Nodes_at_surfElem; % position vector at GP - I use it only as correction for da_vec
            
            
            [~,min_ind_nod]=min(sum((Def_Nodes(Elements_cmgui(epicardial_elements(m_el),Local_node_IDs_in_OPPOSITE_patch),:)-[x_at_GP(1)*ones(length(Local_node_IDs_in_OPPOSITE_patch),1),x_at_GP(2)*ones(length(Local_node_IDs_in_OPPOSITE_patch),1),x_at_GP(3)*ones(length(Local_node_IDs_in_OPPOSITE_patch),1)]).^2,2));
            nearest_nod_opposite_patch=Def_Nodes(Elements_cmgui(epicardial_elements(m_el),Local_node_IDs_in_OPPOSITE_patch(min_ind_nod)),:);
            
            if dot(da_vec,x_at_GP-nearest_nod_opposite_patch)<0 % if da points opposite to the vector pointing from the nearest node at the opposite (epicardial)patch to the endocardial surface GP
                %the above control will only be compromised if the element
                %is very long and thin, so that the two vectors may become
                %perpendicular and the dot product zero, so that numerical
                %roundoffs may make it a very small negative number..-
                %that's why I use the round (10^5*..) extra level of
                %security
                da_vec=-da_vec; % then change direction so da points towards blood cavity
            end
% % %       p_da_per_GP_and_frame_epi(GP_ID_global,:,n_DF)=press*da_vec; % don't need this
            if isempty(Virtual_Disp_Nodes)==0
                virt_u=Nsurf_at_ksi_per_GP_mat(:,n_GP).'*Virtual_Disp_Nodes_at_surfElem;
                p_da_virt_u_w_per_GP_and_frame_epi(GP_ID_global,n_DF)=press*dot(virt_u,da_vec)*surfGPweights_total(n_GP);%p*da \cdot delta_u*w_GP - just sum this and you've got \deltaWext
                da_virt_u_w_per_GP_and_frame_epi(GP_ID_global,n_DF)=dot(virt_u,da_vec)*surfGPweights_total(n_GP);%p*da \cdot delta_u*w_GP - just sum this and you've got \deltaWext                                       
            end
        end                
    end
    %%
    % now that p*da and x have been estimated at the previous and the
    % current frame (at n_DF==1 the previous one is zero so no estimation
    % is needed you can calculate the parts of Wext:
    if n_DF==1
        deltax_mat=squeeze(x_per_GP_and_frame(:,:,n_DF)); % row :GP col: x,y,z
        p_da_average_mat=0.5*squeeze(p_da_per_GP_and_frame(:,:,n_DF)); % row :GP col: x,y,z
        Wext_endo_incr_per_DF(n_DF,:)=sum(sum(p_da_average_mat.*deltax_mat,2).*Global_GP_weights);
    else
        deltax_mat=(squeeze(x_per_GP_and_frame(:,:,n_DF))-squeeze(x_per_GP_and_frame(:,:,n_DF-1))); % row :GP col: x,y,z
        p_da_average_mat=0.5*(squeeze(p_da_per_GP_and_frame(:,:,n_DF))+squeeze(p_da_per_GP_and_frame(:,:,n_DF-1))); % row :GP col: x,y,z
        Wext_endo_incr_per_DF(n_DF,:)=sum(sum(p_da_average_mat.*deltax_mat,2).*Global_GP_weights);
    end
    Wext_endo_per_DF(n_DF,:)=sum(Wext_endo_incr_per_DF(1:n_DF,:));
end
if isempty(Virtual_Disp_Nodes)==0
    virt_Wext_endo_per_DF=sum(p_da_virt_u_w_per_GP_and_frame,1);
    virt_Wext_epi_per_DF=sum(p_da_virt_u_w_per_GP_and_frame_epi,1);
    virt_Wext_epi_over_press_per_DF=sum(da_virt_u_w_per_GP_and_frame_epi,1);
    if length(virt_Wext_endo_per_DF)~=length(diastolic_frames_list)
        disp('size of virt Wext is not what it should be -did you sum along the right direction?');
    end
else
    virt_Wext_endo_per_DF=[];
    virt_Wext_epi_per_DF=[];
    virt_Wext_epi_over_press_per_DF=[];
end




