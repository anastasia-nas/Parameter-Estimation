function [alpha_EB_CF,alpha_EB_pV_CF,alpha_VWB_CF,C1_EB_CF,C1_EB_pV_CF,C1_VWB_CF,EB_CF,EB_pV_CF,VWB_CF]=calculations_for_EB_and_VWB_functional(reference_frame,Diastolic_Frames_Indices_TU,diast_pressures_Pa_TU,Ref_mesh_Nodes,Ref_mesh_Elements_ala_Cmgui,Ref_mesh_Boundaries,Fib_mesh_Nodes,Fib_mesh_Elements_ala_Cmgui,alpha_PS_vector,r_ff,r_fs,r_sn,Def_Nodes_per_Frame_struct,virt_disp_Nodes_per_Frame_struct,Ref_cavity_volume,Def_cavity_volume_struct,Endo_patch_ID,Epi_patch_ID,flag_calc_virt_work_epi,GP_per_elem_dir,Nodes_per_elem_dir)

% created by Anastasia on 18-04-2020
%-------------------------------------------------------------------------------------------------------------------------
% Purpose: 
%-------------------------------------------------------------------------------------------------------------------------
% to have a solid piece of code that doesn't repeat in order to avoid
% inadvertent bugs..
%-------------------------------------------------------------------------------------------------------------------------
% Notes on methodology:
%-------------------------------------------------------------------------------------------------------------------------
% virtual work and energy based functional virt and regular: Wext is just \int_a p*cdot(da,\deltau) or \int_a\int_u pcdot(da,du) -- so just endo
% also include a p*dV based Wext estimation which more robust for clinical images.

%-------------------------------------------------------------------------------------------------------------------------
% Input:
%-------------------------------------------------------------------------------------------------------------------------
% reference_frame: reference frame number
% Diastolic_Frames_Indices_TU
% diast_pressures_Pa_TU
% Ref_mesh_Nodes
% Ref_mesh_Elements_ala_Cmgui
% Ref_mesh_Boundaries
% Fib_mesh_Nodes
% Fib_mesh_Elements_ala_Cmgui: for continuous fiber meshes this is equal to Ref_mesh_Elements_ala_Cmgui
% alpha_PS_vector
% r_ff
% r_fs
% r_sn
% Def_Nodes_per_Frame_struct
% virt_disp_Nodes_per_Frame_struct: if left empty no VF calculations will take place and the respective results based on VWB_CF will return empty variables
% Ref_cavity_volume
% Def_cavity_volume_struct
% Endo_patch_ID
% Epi_patch_ID
% flag_calc_virt_work_epi
% GP_per_elem_dir
% Nodes_per_elem_dir

%-------------------------------------------------------------------------------------------------------------------------
% Output:
%-------------------------------------------------------------------------------------------------------------------------
% structures containing the alpha, C1 and functionals at the current
% reference frame
% alpha_EB_CF:
% alpha_EB_pV_CF
% alpha_VWB_CF
% C1_EB_CF
% C1_EB_pV_CF
% C1_VWB_CF
% EB_CF
% EB_pV_CF
% VWB_CF
%-------------------------------------------------------------------------------------------------------------------------
% SOS points:
%-------------------------------------------------------------------------------------------------------------------------


%% start main body::
% Wint_wrt_alpha: matrix of Wint per alpha in alpha_PS_vector per frame (size: [length(alpha_PS_vector),length(diastolic_frames_list)])
% virtual_Wint_wrt_alpha: matrix of \deltaWint per alpha in alpha_PS_vector per frame (size: [length(alpha_PS_vector),length(diastolic_frames_list)])
[Wint_wrt_alpha_per_frame, virtual_Wint_wrt_alpha_per_frame]=calc_Wint_and_virtWint_per_frame_and_a_for_CF_Guccione(alpha_PS_vector,r_ff,r_fs,r_sn,Diastolic_Frames_Indices_TU,Ref_mesh_Nodes,Ref_mesh_Elements_ala_Cmgui,Fib_mesh_Nodes,Fib_mesh_Elements_ala_Cmgui,Def_Nodes_per_Frame_struct,virt_disp_Nodes_per_Frame_struct,GP_per_elem_dir,Nodes_per_elem_dir);

Wext_pV_per_DF=calc_Wext_with_pV(Diastolic_Frames_Indices_TU,diast_pressures_Pa_TU,Ref_cavity_volume,Def_cavity_volume_struct);

% Wext_endo_per_DF: vector of Wext at endo per diastolic frame - column vector - size(length(diastolic_frames_list),1)
% virt_Wext_endo_per_DF: vector of virtual Wext at endo per diastolic frame - column vector - size(length(diastolic_frames_list),1)
%         [Wext_endo_per_DF, virtual_Wext_endo_per_DF]=calc_Wext_endo_and_virtWext_endo_per_frame_for_CF(Diastolic_Frames_Indices_TU,diast_pressures_Pa_TU,Ref_mesh_Elements_ala_Cmgui,Ref_mesh_Boundaries,Endo_patch_ID,Def_Nodes_per_Frame_struct,virt_disp_Nodes_per_Frame_struct,GP_per_elem_dir,Nodes_per_elem_dir);
[Wext_endo_per_DF, virtual_Wext_endo_per_DF,virtual_Wext_epi_per_DF,virtual_Wext_epi_over_press_per_DF]=calc_Wext_endo_and_virtWext_endo_per_frame_for_CF_verbose_epi(Diastolic_Frames_Indices_TU,diast_pressures_Pa_TU,Ref_mesh_Elements_ala_Cmgui,Ref_mesh_Boundaries,Endo_patch_ID,Epi_patch_ID,Def_Nodes_per_Frame_struct,virt_disp_Nodes_per_Frame_struct,GP_per_elem_dir,Nodes_per_elem_dir,flag_calc_virt_work_epi);
for n_DF2=2:length(Diastolic_Frames_Indices_TU)
    DF2=Diastolic_Frames_Indices_TU(n_DF2);
    for n_DF1=1:(n_DF2-1)
        DF1=Diastolic_Frames_Indices_TU(n_DF1);

        %energy based CF -Wext endo:
        EBfunctional=abs(Wint_wrt_alpha_per_frame(:,n_DF1)./Wint_wrt_alpha_per_frame(:,n_DF2)-ones(size(Wint_wrt_alpha_per_frame,1),1)*Wext_endo_per_DF(n_DF1)/Wext_endo_per_DF(n_DF2));
        EB_CF.(['Ref_Fr_',num2str(reference_frame),'_DF1_',num2str(DF1),'_DF2_',num2str(DF2)])=EBfunctional;
        [~,min_ind_func]=min(EBfunctional);
        alpha_EB_CF.(['Ref_Fr_',num2str(reference_frame),'_DF1_',num2str(DF1),'_DF2_',num2str(DF2)])=alpha_PS_vector(min_ind_func);
        C1_EB_CF.(['Ref_Fr_',num2str(reference_frame),'_DF1_',num2str(DF1),'_DF2_',num2str(DF2)])=Wext_endo_per_DF(n_DF2)/Wint_wrt_alpha_per_frame(min_ind_func,n_DF2);

        %energy based CF -pV:
        EBfunctional_pV=abs(Wint_wrt_alpha_per_frame(:,n_DF1)./Wint_wrt_alpha_per_frame(:,n_DF2)-ones(size(Wint_wrt_alpha_per_frame,1),1)*Wext_pV_per_DF(n_DF1)/Wext_pV_per_DF(n_DF2));
        EB_pV_CF.(['Ref_Fr_',num2str(reference_frame),'_DF1_',num2str(DF1),'_DF2_',num2str(DF2)])=EBfunctional_pV;
        [~,min_ind_func]=min(EBfunctional_pV);
        alpha_EB_pV_CF.(['Ref_Fr_',num2str(reference_frame),'_DF1_',num2str(DF1),'_DF2_',num2str(DF2)])=alpha_PS_vector(min_ind_func);
        C1_EB_pV_CF.(['Ref_Fr_',num2str(reference_frame),'_DF1_',num2str(DF1),'_DF2_',num2str(DF2)])=Wext_pV_per_DF(n_DF2)/Wint_wrt_alpha_per_frame(min_ind_func,n_DF2);

        
        %virtual works based CF:
        if isempty(virt_disp_Nodes_per_Frame_struct)==0
            VWBfunctional=abs(virtual_Wint_wrt_alpha_per_frame(:,n_DF1)./virtual_Wint_wrt_alpha_per_frame(:,n_DF2)-ones(size(Wint_wrt_alpha_per_frame,1),1)*virtual_Wext_endo_per_DF(n_DF1)/virtual_Wext_endo_per_DF(n_DF2));
            VWB_CF.(['Ref_Fr_',num2str(reference_frame),'_DF1_',num2str(DF1),'_DF2_',num2str(DF2)])=VWBfunctional;
            [~,min_ind_func]=min(VWBfunctional);
            alpha_VWB_CF.(['Ref_Fr_',num2str(reference_frame),'_DF1_',num2str(DF1),'_DF2_',num2str(DF2)])=alpha_PS_vector(min_ind_func);
            C1_VWB_CF.(['Ref_Fr_',num2str(reference_frame),'_DF1_',num2str(DF1),'_DF2_',num2str(DF2)])=virtual_Wext_endo_per_DF(n_DF2)/virtual_Wint_wrt_alpha_per_frame(min_ind_func,n_DF2);
        else
            VWB_CF=[];
            alpha_VWB_CF=[];
            C1_VWB_CF=[];
        end

    end            
end