close all; clear all; clc; 

% created by Anastasia on 24-02-2020
%-------------------------------------------------------------------------------------------------------------------------
% Purpose: 
%-------------------------------------------------------------------------------------------------------------------------
% analyse deformation and pressure data to get C1-alpha parameters in
% Guccione. For original energy based CF formulation and VF based. C1 in
% this script is estimated from energy balance and this can serve as an
% index of the robustness..
%-------------------------------------------------------------------------------------------------------------------------
% Notes on methodology:
%-------------------------------------------------------------------------------------------------------------------------
% use create_cheart_prep_setups_new_cases_hole_new.m first to create the
% Virtual fields at each frame of the data. For this you can use the
% following options (still searching which is the best..)
% either : flag_create_VF_just_base_and_hole_NeoHook=1 (for neohookean material
% description)
% or: flag_create_VF_just_base_and_hole_ISO_Gucc=1 AND flag_fake_fibers=1 (for
% Guccione with isotropic parameters and "fake fibers" where ef=[1,0,0],es=[0,1,0] and en=[0,0,1]
%....................................
% also run create_cheart_prep_setups_new_cases_hole_new.m (with flag flag_prep_for_parameter_estimation=1)
% to create the folder and data structure that I will use here (Saving deformation in
% Cheart format and pressure from data (modified so that I stick to the
% convention I had made before of offseting pressure so that 0 sits at the
% reference frame, and reference frame is frame with minimum pressure).

%-------------------------------------------------------------------------------------------------------------------------
% Input:
%-------------------------------------------------------------------------------------------------------------------------


%-------------------------------------------------------------------------------------------------------------------------
% Output:
%-------------------------------------------------------------------------------------------------------------------------

%-------------------------------------------------------------------------------------------------------------------------
% SOS points:
%-------------------------------------------------------------------------------------------------------------------------

%================================================================================================================================
%% FLAGS:
flag_Joshuas_laptop=1; % Joshua leave this to =1..
%% simulation specification flags (they specify which Synthetic data set I use as input here)
flag_from_simulation=0; 
flag_prescr_basal_disps=0; 
flag_prescr_apicalHole_disps=0; 
flag_fix_apex=0; 
flag_press_epi=0;   % if you're analysing deformation from a synthetic data set - flag_prescr_basal_disps=1 if it's the simulation from prescribing basal disps you are focusing on (if flag_prescr_basal_disps==0 then it's the fixed base simulations)
%% or if you analyse data:
flag_data=1; % =1 if you are processing the actual clinical data
Lagrange_meshes_interp_order='Cubic' ;  % choose between 'Cubic' and 'Quadratic'

%% i_case:
i_case='JUL_new'; res_str='9124'; %this needs to be '341' or '9124'
flag_correct_pressure_data=1; % leave this =1 for now - this involves using corrected pressures (with no negatives) as in BMMB.
%% prep for data analysis for parameter estimation (print boundaries, neighbours and used pressure per frame after correction according to pipeline (where min pressure is at reference frame).
% flag_prep_for_parameter_estimation=1;% this is equivalent to the flag_data here.. % this option will override all others 
%% choose fiber fitting options if your material is anisotropic. Option flag_fake_fibers is to be used ONLY for VF generation
flag_my_fibers=1; flag_Daves_fibers=0; flag_fake_fibers=0;% flag_my_fibers=1 for using my tool-flag_Daves_fibers=1 for preping folder to run cheart fibre fitting tool, or both can be zero if you have an isotropic material
%% choose these for generating VF: (note flag_prescr_basal_disps must be ==1!)
%note if mesh has hole you also need to fix hole to avoid having work doneat the hole by the stresses:
flag_use_VF_just_base_and_hole_NeoHook=0; % =1 if this script is run to prepare VF generation with fixed base and neohookean material (for any other type of VF generating simulation - e.g VF for epi/Stokes/isotropic Guccione/ whatever- create a new flag describing the case!)- this didn't work well in Cheart (sth wrong with NeoHookean?)
flag_use_VF_just_base_and_hole_ISO_Gucc=0;
flag_create_VF_epi_NeoHookean_sin_mu3=0;

%% SPECS:
fig_ind=[]; % for numbering plots
Endo_patch_ID=1;
Epi_patch_ID=2;
Base_patch_ID=3;
Apex_hole_patch_ID=4; % =[] if no apical hole exists in mesh.
flag_disc_Fib=0; % =1 if discontinuous fibers are used and =0 otherwise. For meshes with hole 
fiber_epi_degrees=-60;fiber_endo_degrees=60;
% material parameters default:
r_ff=0.55; r_fs=0.25; r_sn=0.2; % change them further down if you need an isotropic material for example
alpha_PS_vector=[-0.2:0.05:0,1:1:100,105:5:500];
flag_calc_virt_work_epi=1; % to estimate work at epi to check if VF epi works..
sim_with_hole=0; sim_capped_apex=0; sim_collapsed=0; % preallocation - mesh description will be specified under each "i_case"
%=============================================================================================================================================================
%% catching quick errors about contradicting flag choices::
if sum([flag_from_simulation,flag_data])==0 || sum([flag_from_simulation,flag_data])>1
    warning('ERROR: these flags are mutually exclusive - either you import deformation from data or from a simulation')
    warning('Stopping script "Parameter_estimation_Guccione_energy_based_VW_based_CF" execution now!')
   return
end

if sum([flag_use_VF_just_base_and_hole_NeoHook,flag_use_VF_just_base_and_hole_ISO_Gucc,flag_create_VF_epi_NeoHookean_sin_mu3])>1
   warning('ERROR: these flags are mutually exclusive as they correspond to completely different actions. Set one to 1 and the rest to 0 to proceed!')
   warning('Stopping script "Parameter_estimation_Guccione_energy_based_VW_based_CF" execution now!')
   return
elseif sum([flag_use_VF_just_base_and_hole_NeoHook,flag_use_VF_just_base_and_hole_ISO_Gucc,flag_create_VF_epi_NeoHookean_sin_mu3])==0
    warning('!! No Virtual Field selected --proceeding without VW based analysis!!')
    warning('....................................................................')
end
% for new cases:
if any(strcmp(i_case,{'BAL','HAW','JER_new','JUL_new','MCP','OSU','ROL','SAN','WIL'}))==1
 %%
     

   
    % FOLDERS:
    if flag_Joshuas_laptop==1
        folder_pre_fix='C:\Users\fitzp\Documents\Meshes';
    else
        folder_pre_fix='C:\Users\an11\Desktop\Work\Teaching\Stiffness_Project_Joshua_Kader';
    end
    if strcmp(res_str,'341')==1
        gen_folder=[folder_pre_fix,'\Manav_Cases_Meshed_341']; % use your _341 or _9124 folders here
    elseif strcmp(res_str,'9124')==1
        gen_folder=[folder_pre_fix,'\Manav_Cases_Meshed_9124']; % use your _341 or _9124 folders here
    end
    pres_frame_data_folder_intro_new_cases=[folder_pre_fix,'\new_cases_pressures_JK'];% from laptop
    

%     % mesh/data
%     Lagrange_meshes_interp_order='Cubic' ;  % choose between 'Cubic' and 'Quadratic'
%     flag_rescale_mesh_to_meters=1; % if you use the meshes converted to m
%     flag_correct_pressure_data=1; % leave like that here
%     % FLAGS:
%%
%      % FOLDERS:
%     if flag_run_from_laptop==1
%         % from laptop::
%         gen_folder='C:\Users\fitzp\Documents\Meshes\Manav_Cases_Meshed_341'; % use from laptop    
%     else
%         % from bioeng139::
%         gen_folder='/staff/an11/new_cases_meshed'; % from bioeng-139
%     end

    switch Lagrange_meshes_interp_order
        case 'Cubic' % CUBIC LAGRANGE MESHES
            Lag_mesh_general_folder=[gen_folder,'/Cubic_Lagrange_Meshes_in_m'];
            VF_Simulations_general_folder=[gen_folder,'/Cheart_simulations_for_VF_generation/Cubic_Lagrange_meshes'];
            Data_analysis_general_folder=[gen_folder,'/Data_prep_for_param_estimation/Cubic_Lagrange_meshes'];
            Simulations_general_folder=[gen_folder,'/Cheart_simulations_setup_folders/Cubic_Lagrange_meshes'];
            Nodes_per_elem_dir=4;

        case 'Quadratic' % QUADRATIC LAGRANGE MESHES
            Lag_mesh_general_folder=[gen_folder,'/Quadratic_Lagrange_Meshes_in_m'];
            VF_Simulations_general_folder=[gen_folder,'/Cheart_simulations_for_VF_generation/Quadratic_Lagrange_meshes'];
            Data_analysis_general_folder=[gen_folder,'/Data_prep_for_param_estimation/Quadratic_Lagrange_meshes'];
            Simulations_general_folder=[gen_folder,'/Cheart_simulations_setup_folders/Quadratic_Lagrange_meshes'];
            Nodes_per_elem_dir=3;
    end
    GP_per_elem_dir=Nodes_per_elem_dir;
    % as I lost all my data in auto/nas/heart-nas, I use my bioeng003 saved data for this. I think the PAS_NF (passive inflation - no filter) is
    % the one I used in BMMB..
    % (e.g. /home/baddisk/an11/Desktop/new_cases_processing_folder_PAS_NF/case_BAL/Parameter_Estimation_plots/From_Matlab/chosen_offset_38/case_BAL)
    
    switch i_case
        case 'ADEGBES' % not used in BMMB
            disp('this case was not included in BMMB - are you sure you want to run it?');
            i_case_pres=i_case;
            press_data_folder='SA';
            chosen_offset=32; % stin tyxi to dialexa den kserw an einai swsto.. NOv 2018 - check with pv sync script-> vgazei 32 kai to xrismopooiw
            chosen_offset=46; % to 32 den mou to kanei print sto folder gia kapoio logo opote as kanw debug to script me to 46 kai epistrefw an thelw (allwste to case ADE den itan katholou kalo)        
            mesh_root_name='MNaNMesh9124';
            ref_frame_list=[];
            
        case 'BAL'
            i_case_pres=i_case;
%             chosen_offset=75;
            chosen_offset=38; % use this
            press_data_folder='MB';
            mesh_root_name='MNaNMesh9124';
            ref_frame_list=16; %[15:17]; % default reference frame: 16
            C1=5300; % put the identified values on for the forward simulation, will only be used for the Cheart Pfile generation
            alpha=61;  
            all_frames_ind=1:25;
            
        case 'HAW'
            i_case_pres=i_case;
%             chosen_offset=78;
            chosen_offset=39; % use this
            press_data_folder='GH';
            mesh_root_name='MNaNMesh9124';
            ref_frame_list=23; %[22:24]; % default reference frame is 23
            C1= 820; % put the identified values on for the forward simulation, will only be used for the Cheart Pfile generation
            alpha=61;
            all_frames_ind=1:30;

        case 'JER_new'
            i_case_pres='JER';
%             chosen_offset=41; % this threw error that file is missing - maybe this is missing because the copy of files was performed before i ran offset 41
            chosen_offset=46; % use this one
            press_data_folder='MJ';
            mesh_root_name='MNaNMesh9124';
            ref_frame_list=20; %[19:21]; % default reference frame is 20 
            C1= 4780; % put the identified values on for the forward simulation, will only be used for the Cheart Pfile generation
            alpha= 5;
            all_frames_ind=1:35;

        case 'JUL_new'
            i_case_pres='JUL';
%             chosen_offset=20;
%             chosen_offset=17;
            chosen_offset=19; % use this one gia syncing me JUL_new mesh pou xrisimopoiw
            press_data_folder='AJ';
            mesh_root_name='MNaNMesh9124';
            ref_frame_list=27; % [26:28]; % default reference frame is 27 
            C1= 1960; % put the identified values on for the forward simulation, will only be used for the Cheart Pfile generation
            alpha= 66;
            all_frames_ind=1:35;

        case 'MCP' % not used in BMMB
            disp('this case was not included in BMMB - are you sure you want to run it?');
            i_case_pres=i_case;
            chosen_offset=75; 
            chosen_offset=41; 
            press_data_folder='SMc';
            mesh_root_name='MNaNMesh9124';
            ref_frame_list=[];
            C1=[]; % put the identified values on for the forward simulation, will only be used for the Cheart Pfile generation
            alpha=[];
            all_frames_ind=1:35;

        case 'OSU' % not forward simulation could not converge for collapsed elem mesh
            i_case_pres=i_case;
            chosen_offset=22; 
            press_data_folder='EOS';
            mesh_root_name='MNaNMesh9124';
            ref_frame_list=[23:25]; % default reference frame is 24
            C1= 3140; % put the identified values on for the forward simulation, will only be used for the Cheart Pfile generation
            alpha=24; 
            all_frames_ind=1:35;

        case 'ROL' % not used in BMMB
            disp('this case was not included in BMMB - are you sure you want to run it?');
            i_case_pres=i_case;
            chosen_offset=60; 
            chosen_offset=52;
            press_data_folder='OR';
            mesh_root_name='MNaNMesh9124';
            ref_frame_list=[];
            C1=[]; % put the identified values on for the forward simulation, will only be used for the Cheart Pfile generation
            alpha=[];
            all_frames_ind=1:35;

        case 'SAN'
            i_case_pres=i_case;
            chosen_offset=32; 
            press_data_folder='SS';
            mesh_root_name='MNaNMesh9124';
            ref_frame_list=[17:19]; % default reference frame is 18 
            C1= 1460; % put the identified values on for the forward simulation, will only be used for the Cheart Pfile generation
            alpha= 29;
            all_frames_ind=1:35;

        case 'WIL'
            i_case_pres=i_case;
            chosen_offset=46; 
            press_data_folder='MW';
            mesh_root_name='MNaNMesh9124';
            ref_frame_list=23; %[22:24]; % default reference frame is 23 
            C1= 3300; % put the identified values on for the forward simulation, will only be used for the Cheart Pfile generation
            alpha= 66;
            all_frames_ind=1:35;
    end
    
    %% continue with building the folder paths:
    % VF folder path:
    if flag_use_VF_just_base_and_hole_NeoHook==1
        VF_mesh_folder_path_intro=[VF_Simulations_general_folder,'/Fixed_Base_and_ApexHole_NeoHookean'];
    elseif flag_use_VF_just_base_and_hole_ISO_Gucc==1
        VF_mesh_folder_path_intro=[VF_Simulations_general_folder,'/Fixed_Base_and_ApexHole_ISOTROPIC_Guccione'];
    elseif flag_create_VF_epi_NeoHookean_sin_mu3==1
        VF_mesh_folder_path_intro=[VF_Simulations_general_folder,'/Prescr_Tang_Epi_sinmu3_NeoHookean'];
    else
        VF_mesh_folder_path_intro=[];
    end
    if isempty(VF_mesh_folder_path_intro)==0
        VF_mesh_folder_path_intro=[VF_mesh_folder_path_intro,'/case_',i_case];
    end
    
    % Mesh and pressure data folder:
    Mesh_and_pres_folder=[Data_analysis_general_folder,'/Cheart_Mesh_and_fitted_fibers_at_frame'];    
    if any([flag_my_fibers, flag_Daves_fibers, flag_fake_fibers])==0  || sum([flag_my_fibers, flag_Daves_fibers, flag_fake_fibers])>1
        warning('exactly one fiber fitting method must be chosen - set either "flag_my_fibers" or "flag_Daves_fibers" or "flag_fake_fibers" to 1')
                disp('stopping script.....')
                return
    else % My_Fibers Daves_Fibers 
        if flag_my_fibers==1
            Mesh_and_pres_folder=[Mesh_and_pres_folder,'/My_Fibers/Endo_',num2str(fiber_endo_degrees),'_epi_',num2str(abs(fiber_epi_degrees))];           
        elseif flag_Daves_fibers==1 % I am currently removing this option from the script as it is faster to run my code, works fine and I wanted to remove reduntant flags
            Mesh_and_pres_folder=[Mesh_and_pres_folder,'/Daves_Fibers/Endo_',num2str(fiber_endo_degrees),'_epi_',num2str(abs(fiber_epi_degrees))];          
        elseif flag_fake_fibers==1 % 
            Mesh_and_pres_folder=[Mesh_and_pres_folder,'/Fake_Fibers_XYZ'];
        end
    end    
    Mesh_and_pres_folder=[Mesh_and_pres_folder,'/case_',i_case];
    
    %% also sort folders for the case that you analyse synthetic datasets
    % from simulations with prescr. basal disps for example:   
    if flag_prescr_basal_disps==1  
        if flag_prescr_apicalHole_disps==1
            Simulations_mesh_folder_path=[Simulations_general_folder,'/Prescribe_Disps_at_Base_Apex_Nodes'];
        elseif flag_fix_apex==1
            Simulations_mesh_folder_path=[Simulations_general_folder,'/Prescribe_Disps_at_Base_Nodes_Fix_Apex'];
        else
            Simulations_mesh_folder_path=[Simulations_general_folder,'/Prescribe_Disps_at_Base_Nodes'];            
        end           
    else % this corresponds to preparing Cheart sims for running like in " flag_prescr_basal_disps=1 " but with fixed base (no prescribed basal disps)
        if flag_prescr_apicalHole_disps==1
            Simulations_mesh_folder_path=[Simulations_general_folder,'/Fixed_Base_Prescr_Apex']; % this combination is not that useful as a fixed base is usually chosen for control and prescribing apex is not that important as is prescribing base, so no real point in trying this other than for numerical reasons (checking which boundary causes simulation to non-coverge)
        elseif flag_fix_apex==1
            Simulations_mesh_folder_path=[Simulations_general_folder,'/Fixed_Base_Apex'];
        else
            Simulations_mesh_folder_path=[Simulations_general_folder,'/Fixed_Base'];
        end        
    end
    if flag_press_epi==1 % this flag only applies to the simulations for creating fake data
        Simulations_mesh_folder_path=[Simulations_mesh_folder_path,'_PresEpi'];
    end
    
    if flag_my_fibers==1 || flag_Daves_fibers==1 || flag_fake_fibers==1
        if flag_my_fibers==1
            Simulations_mesh_folder_path=[Simulations_mesh_folder_path,'/My_Fibers/Endo_',num2str(fiber_endo_degrees),'_epi_',num2str(abs(fiber_epi_degrees))]; %folder paths are too close to max path length - economize!          
        elseif flag_Daves_fibers==1 % I am currently removing this option from the script as it is faster to run my code, works fine and I wanted to remove reduntant flags
            Simulations_mesh_folder_path=[Simulations_mesh_folder_path,'/Daves_Fibers/Endo_',num2str(fiber_endo_degrees),'_epi_',num2str(abs(fiber_epi_degrees))];
        elseif flag_fake_fibers==1 % 
            Simulations_mesh_folder_path=[Simulations_mesh_folder_path,'/Fake_Fibers_XYZ'];
        end
    end
    Simulations_mesh_folder_path=[Simulations_mesh_folder_path,'/case_',i_case];            
    
    for n_ref_fr=1:length(ref_frame_list)
        reference_frame=ref_frame_list(n_ref_fr);
        
        Ref_mesh_folder=[Mesh_and_pres_folder,'/Ref_Frame_',num2str(reference_frame)];
        if flag_correct_pressure_data==1
            Pressure_data_folder_path=[Ref_mesh_folder,'/Corrected_pressures_no_neg'];                    
        end
        press_data_file=[Pressure_data_folder_path,'/frame_and_press.data'];
        data_pre=[];
        fid=fopen(press_data_file);
        tot_diast_frames=fscanf(fid,'%d',[1,1]); % number of frames in total
        data_pre=fscanf(fid,'%d     %f',[2,Inf]).';
        fclose(fid);
        Diastolic_Frames_Indices_TU=data_pre(:,1); 
        diast_pressures_Pa_TU=data_pre(:,2);
        
        fid=fopen([Ref_mesh_folder,'/Ref_mesh.T']);
        Element_intro=fscanf(fid,'%d',[2,1]).';
        Ref_mesh_Elements_ala_CHeart=fscanf(fid,'%d',[Nodes_per_elem_dir^3,Inf]).';
        fclose(fid);
        
        [Ref_mesh_Elements_ala_Cmgui,error]=turn_Cheart_node_ordering_to_Cmgui(Ref_mesh_Elements_ala_CHeart,Nodes_per_elem_dir);
        
        fid=fopen([Ref_mesh_folder,'/Ref_mesh.X']);
        Nodes_intro=fscanf(fid,'%d',[2,1]).';
        Ref_mesh_Nodes=fscanf(fid,'%f',[3,Inf]).';
        fclose(fid);
        
        fid=fopen([Ref_mesh_folder,'/Ref_mesh.B']);
        Bound_intro=fscanf(fid,'%d',[1,1]).';
        Ref_mesh_Boundaries=fscanf(fid,'%d',[Nodes_per_elem_dir^2+2,Inf]).';
        fclose(fid);
        
        fid=fopen([Ref_mesh_folder,'/Ref_mesh.NBT']);
        Neigh_intro=fscanf(fid,'%d',[1,1]).';
        Ref_mesh_Neighbours=fscanf(fid,'%d',[6,Inf]).';
        fclose(fid);
                 
        if flag_disc_Fib==0
            fid=fopen([Ref_mesh_folder,'/FIBERS.X']);
            Fib_Nodes_intro=fscanf(fid,'%d',[2,1]).';
            Fib_mesh_Nodes=fscanf(fid,'%f',[9,Inf]).';
            fclose(fid);
        else
            disp('you need to build the discontinuous fiber field option- currently not available for meshes with hole');
        end
        

        for def_fr_ind=1:length(Diastolic_Frames_Indices_TU)
            frame_ID=Diastolic_Frames_Indices_TU(def_fr_ind);
            if flag_data==1
                % read Def_Nodes per frame in a structure
                fid=fopen([Ref_mesh_folder,'/Def_mesh-',num2str(frame_ID),'.X']);
                Nodes_intro=fscanf(fid,'%d',[2,1]).';
                Def_mesh_Nodes=fscanf(fid,'%f',[3,Inf]).';
                fclose(fid);
                
                % and virtual field Disp Nodes per frame in a structure
                if isempty(VF_mesh_folder_path_intro)==0
                    list_DispVF=dir([VF_mesh_folder_path_intro,'/Frame_',num2str(frame_ID),'/Disp-*.D']);
                    end_VF_Disp_file=[VF_mesh_folder_path_intro,'/Frame_',num2str(frame_ID),'/Disp-',num2str(length(list_DispVF)-2),'.D']; %using the output from the last increment in the simulation as that has the largest deformation - if simulation didn't converge though maybe it's best to go 1or 2 increments before the last one to avoid weird deformation occuring.
                    disp('using 2 frames before last converged increment as VF disp to avoid possible issues with bad elements')
                end
                
            elseif flag_from_simulation==1
                Sim_prep_folder=[Simulations_mesh_folder_path,'/Ref_Frame_',num2str(reference_frame),'/Def_Frame_',num2str(frame_ID)]; 
                % check reference frame and pressure is according to
                % pressure data
                fid=fopen([Sim_prep_folder,'/',Lagrange_meshes_interp_order,'_mesh_FE.X']);
                Nodes_intro=fscanf(fid,'%d',[2,1]).';
                Lagrange_mesh_Nodes_Ref=fscanf(fid,'%f',[3,Inf]).';
                fclose(fid);   
                
                if max(max(abs(Lagrange_mesh_Nodes_Ref-Ref_mesh_Nodes)))>10^(-6)
                    disp('wtf?? ref nodes in data and simulation should match!!!???')
                    disp('');
                end
                [dp_dt,max_incrs,flag_found_it_press]=read_dp_dt_from_pressure_expression_in_Pfile(Sim_prep_folder,'Ellipse.P'); % dp_dt is <0 for inflation           
                % find max increments in simulation
                list_Disps_sim=dir([Sim_prep_folder,'/Disp-*.D']);
                if abs(-diast_pressures_Pa_TU(def_fr_ind)-dp_dt*max_incrs)>10^(-5) % diast_pressures_Pa_TU is >0 and hence add the minus sign (-) to give it the same sign as dp_dt.
                    disp('wtf?? pressure in data and simulation should match!!!???')
                    disp([' from data: ',num2str(-diast_pressures_Pa_TU(def_fr_ind)),' from sim: ',num2str(dp_dt*length(list_Disps_sim))])
                    disp('');
                else
                    diast_pressures_Pa_TU(def_fr_ind)=-dp_dt*length(list_Disps_sim); % correct for pressure as the maximum pressure may not have been reached (simulation not converged)-- and the (-dp_dt) is because diast_pressures_Pa_TU must be >0
                end
                % read that disp and add to reference frame for creating
                % deformed state
                fid=fopen([Sim_prep_folder,'/Disp-',num2str(length(list_Disps_sim)),'.D']);
                Nodes_intro=fscanf(fid,'%d',[2,1]).';
                Disp_Nodes_Ref=fscanf(fid,'%f',[3,Inf]).';
                fclose(fid); 
                Def_mesh_Nodes=Disp_Nodes_Ref+Ref_mesh_Nodes;
                
                % and virtual field Disp Nodes per frame in a structure
                if flag_use_VF_just_base_and_hole_NeoHook==1
                    if exist([Sim_prep_folder,'/VF_Frame_',num2str(length(list_Disps_sim))],'dir')==1
                        VF_folder=[Sim_prep_folder,'/VF_Frame_',num2str(length(list_Disps_sim))]; % in this folder the last increment of the simulation with prescr disps is added to reference nodes and this deformed state is used as reference for generating the VF frame with fixed base&hole , so that the VF field satisfies the properties at the deformed configuration as is required                        
                    elseif exist([Sim_prep_folder,'/VFFixBasApNH_',num2str(length(list_Disps_sim))],'dir')==1
                        VF_folder=[Sim_prep_folder,'/VFFixBasApNH_',num2str(length(list_Disps_sim))];
                    else
                        disp('what is going on here? Is there no VF folder specified?');
                    end
                elseif flag_create_VF_epi_NeoHookean_sin_mu3==1                  
                    VF_folder=[Sim_prep_folder,'/VFEpiNH_',num2str(length(list_Disps_sim))];
                end
                
                if isempty(VF_mesh_folder_path_intro)==0
                    list_DispVF=dir([VF_folder,'/Disp-*.D']);
                    end_VF_Disp_file=[VF_folder,'/Disp-',num2str(length(list_DispVF)-2),'.D']; %using the output from the last increment in the simulation as that has the largest deformation - if simulation didn't converge though maybe it's best to go 1or 2 increments before the last one to avoid weird deformation occuring.
                    disp('using 2 frames before last converged increment as VF disp to avoid possible issues with bad elements');
                end
            end
            Def_Nodes_per_Frame_struct.(['Frame_',num2str(frame_ID)])=Def_mesh_Nodes;
            
            % calculate cavity volume at deformed configuration
            [Def_cavity_volume, ~, ~]=calculate_cavity_volume_Lagrange_meshes_with_hole(Ref_mesh_Elements_ala_Cmgui,Def_mesh_Nodes,Ref_mesh_Boundaries,Endo_patch_ID,Epi_patch_ID,Base_patch_ID,Apex_hole_patch_ID,Nodes_per_elem_dir,GP_per_elem_dir,fig_ind);
            Def_cavity_volume_struct.(['Frame_',num2str(frame_ID)])=Def_cavity_volume;
            
            if isempty(fig_ind)==0
                fig_ind=fig_ind+1;
            end
            if isempty(VF_mesh_folder_path_intro)==0
                fid=fopen(end_VF_Disp_file);
                VFNodes_intro=fscanf(fid,'%d',[2,1]).';
                VF_Disp_Nodes=fscanf(fid,'%f',[3,Inf]).';
                fclose(fid);
                virt_disp_Nodes_per_Frame_struct.(['Frame_',num2str(frame_ID)])=VF_Disp_Nodes;  
            else
                virt_disp_Nodes_per_Frame_struct=[];
            end
        end
        [Ref_cavity_volume, ~, ~]=calculate_cavity_volume_Lagrange_meshes_with_hole(Ref_mesh_Elements_ala_Cmgui,Ref_mesh_Nodes,Ref_mesh_Boundaries,Endo_patch_ID,Epi_patch_ID,Base_patch_ID,Apex_hole_patch_ID,Nodes_per_elem_dir,GP_per_elem_dir,fig_ind);
         
        %% start insert function 18/4/2020
%           [alpha_EB_CF_tu,alpha_EB_CF_pV,alpha_VWB_CF,C1_EB_CF_tu,C1_EB_CF_pV,C1_VWB_CF,EB_CF,EB_pV_CF,VWB_CF]=calculations_for_EB_and_VWB_functional(reference_frame,Diastolic_Frames_Indices_TU,diast_pressures_Pa_TU,Ref_mesh_Nodes,Ref_mesh_Elements_ala_Cmgui,Ref_mesh_Boundaries,Fib_mesh_Nodes,Ref_mesh_Elements_ala_Cmgui,alpha_PS_vector,r_ff,r_fs,r_sn,Def_Nodes_per_Frame_struct,virt_disp_Nodes_per_Frame_struct,Ref_cavity_volume,Def_cavity_volume_struct,Endo_patch_ID,Epi_patch_ID,flag_calc_virt_work_epi,GP_per_elem_dir,Nodes_per_elem_dir);
        %% end insert function -18/4/2020
        %% 4/6/2020: pedantic printing of the estimations: outputting Wext and E (elastic energy/internal work)
        [alpha_EB_CF_tu,alpha_EB_CF_pV,alpha_VWB_CF,C1_EB_CF_tu,C1_EB_CF_pV,C1_VWB_CF,EB_CF,EB_pV_CF,VWB_CF,Wint_wrt_alpha_per_frame,Wext_pV_per_DF,Wext_endo_per_DF,~,~,~,~]=calculations_for_EB_and_VWB_functional_pedantic(reference_frame,Diastolic_Frames_Indices_TU,diast_pressures_Pa_TU,Ref_mesh_Nodes,Ref_mesh_Elements_ala_Cmgui,Ref_mesh_Boundaries,Fib_mesh_Nodes,Ref_mesh_Elements_ala_Cmgui,alpha_PS_vector,r_ff,r_fs,r_sn,Def_Nodes_per_Frame_struct,virt_disp_Nodes_per_Frame_struct,Ref_cavity_volume,Def_cavity_volume_struct,Endo_patch_ID,Epi_patch_ID,flag_calc_virt_work_epi,GP_per_elem_dir,Nodes_per_elem_dir);
        
        % print strain energy for C1=1700 a=15 (assuming healthy subject
        % parameters (cite Anastasia BMMB) you get E=\int_V \Psi dV,
        % \Psi=0.5*C1*(exp(Q)-1)
        C1=1700;
        Internal_energy_per_DF=Wint_wrt_alpha_per_frame(20,:)*1700; % row vector of number of frames
        plot_output_mat=[Diastolic_Frames_Indices_TU;Internal_energy_per_DF.';Wext_pV_per_DF.'];
        %plot 3 cols: col1: Diastolic_Frames_Indices_TU; col2:
        %Internal_energy_per_DF; col3: Wext_pV_per_DF;
        
        
        
        
        
        
    end
    
    
    
    
    
    
    
else
    disp('modify script for use with other meshes except for new cases')
    %%
      % FOLDERS:
    if flag_run_from_laptop==1
        % from laptop::
        gen_dir='C:\Users\an11\Desktop\Work\Synthetic_data\For_VF_as_unknown_BCs_solution_method_verification';
    else
        % from bioeng139::
        gen_folder=''; % from bioeng-139
    end
    switch i_case
        case 'Fixed_base_and_hole_PIendo_grad_epi_sept_a_15_C1_1000'
            geometry='Ellipsoid';
            C1=1000; alpha=15;
            r_ff=0.55; r_fs=0.25; r_sn=0.2;
            pmax_endo=-1500; % <0 gia passive inflation sto LV
            pmax_epi=-500; % <0 gia na exw piesi antistoixi tou na antisteketai to epicaridum sto inflation tou LV
            pmax_sept=-750;
            septal_elements=[220:222,232:234,244:246,256:258,268:270,280:282];
            Septum_Boundary_ID=5;
            flag_dp_dZ_0_to_1_endo=1;
            TOT_incrs=30;
            flag_prescr_base_disps=0;
            material_law='Guccione';
            material_parameters=[C1,alpha*r_ff,alpha*r_sn,alpha*r_fs];
            fiber_endo_degrees=60; fiber_epi_degrees=-60;
            % specify the flags for the Cheart file generation.
            flag_prescr_bound_base=0; %=1 if you prescr disps 
            flag_prescr_bound_epi=0;
            flag_prescr_bound_endo=0; 
            flag_prescr_bound_apex_hole=0;
            flag_prescr_APEX_NODE=0;
            flag_disc_Fib=0;
            flag_fix_base=1;
            flag_fix_apex_hole=1; 
            %% VF specs:
            flag_also_create_VF_epi=1; n_incr_VF=30;
            flag_continuous_fibers=1;

        case 'Fixed_base_and_hole_PIendo_epi_sept_a_15_C1_1000'
            geometry='Ellipsoid';
            C1=1000; alpha=15;
            r_ff=0.55; r_fs=0.25; r_sn=0.2;
            pmax_endo=-1500; % <0 gia passive inflation sto LV
            pmax_epi=-500; % <0 gia na exw piesi antistoixi tou na antisteketai to epicaridum sto inflation tou LV    
            pmax_sept=-750;
            septal_elements=[220:222,232:234,244:246,256:258,268:270,280:282];
            Septum_Boundary_ID=5;
            TOT_incrs=30;
            flag_prescr_base_disps=0;
            material_law='Guccione';
            material_parameters=[C1,alpha*r_ff,alpha*r_sn,alpha*r_fs];
            fiber_endo_degrees=60; fiber_epi_degrees=-60;
            % specify the flags for the Cheart file generation.
            flag_prescr_bound_base=0; %=1 if you prescr disps 
            flag_prescr_bound_epi=0;
            flag_prescr_bound_endo=0; 
            flag_prescr_bound_apex_hole=0;
            flag_prescr_APEX_NODE=0;
            flag_disc_Fib=0;
            flag_fix_base=1;
            flag_fix_apex_hole=1;  
            %% VF specs:
            flag_also_create_VF_epi=1; n_incr_VF=30;
            flag_continuous_fibers=1;

        case 'Fixed_base_and_hole_PIendo_grad_epi_a_15_C1_1000'
            geometry='Ellipsoid';
            C1=1000; alpha=15;
            r_ff=0.55; r_fs=0.25; r_sn=0.2;
            pmax_endo=-1500; % <0 gia passive inflation sto LV
            pmax_epi=-500; % <0 gia na exw piesi antistoixi tou na antisteketai to epicaridum sto inflation tou LV
            flag_dp_dZ_0_to_1_endo=1;
            TOT_incrs=30;
            flag_prescr_base_disps=0;
            material_law='Guccione';
            material_parameters=[C1,alpha*r_ff,alpha*r_sn,alpha*r_fs];
            fiber_endo_degrees=60; fiber_epi_degrees=-60;
            % specify the flags for the Cheart file generation.
            flag_prescr_bound_base=0; %=1 if you prescr disps 
            flag_prescr_bound_epi=0;
            flag_prescr_bound_endo=0; 
            flag_prescr_bound_apex_hole=0;
            flag_prescr_APEX_NODE=0;
            flag_disc_Fib=0;
            flag_fix_base=1;
            flag_fix_apex_hole=1;  
            %% VF specs:
            flag_also_create_VF_epi=1; n_incr_VF=30;
            flag_continuous_fibers=1;

        case 'Fixed_base_and_hole_PIendo_epi_a_15_C1_1000'
            geometry='Ellipsoid';
            C1=1000; alpha=15;
            r_ff=0.55; r_fs=0.25; r_sn=0.2;
            pmax_endo=-1500; % <0 gia passive inflation sto LV
            pmax_epi=-500; % <0 gia na exw piesi antistoixi tou na antisteketai to epicaridum sto inflation tou LV        
            TOT_incrs=30;
            flag_prescr_base_disps=0;
            material_law='Guccione';
            material_parameters=[C1,alpha*r_ff,alpha*r_sn,alpha*r_fs];
            fiber_endo_degrees=60; fiber_epi_degrees=-60;
            % specify the flags for the Cheart file generation.
            flag_prescr_bound_base=0; %=1 if you prescr disps 
            flag_prescr_bound_epi=0;
            flag_prescr_bound_endo=0; 
            flag_prescr_bound_apex_hole=0;
            flag_prescr_APEX_NODE=0;
            flag_disc_Fib=0;
            flag_fix_base=1;
            flag_fix_apex_hole=1; 
            %% VF specs:
            flag_also_create_VF_epi=1; n_incr_VF=30;
            flag_continuous_fibers=1;

        case 'Fixed_base_and_hole_PIendo_a_15_C1_1000'
            geometry='Ellipsoid';
            C1=1000; alpha=15;
            r_ff=0.55; r_fs=0.25; r_sn=0.2;
            pmax_endo=-1500;
            TOT_incrs=30;
            flag_prescr_base_disps=0;
            material_law='Guccione';
            material_parameters=[C1,alpha*r_ff,alpha*r_sn,alpha*r_fs];
            fiber_endo_degrees=60; fiber_epi_degrees=-60;
            % specify the flags for the Cheart file generation.
            flag_prescr_bound_base=0; %=1 if you prescr disps 
            flag_prescr_bound_epi=0;
            flag_prescr_bound_endo=0; 
            flag_prescr_bound_apex_hole=0;
            flag_prescr_APEX_NODE=0;
            flag_disc_Fib=0;
            flag_fix_base=1;
            flag_fix_apex_hole=1;   
            %% VF specs:
            flag_also_create_VF_epi=1; n_incr_VF=30;
            flag_continuous_fibers=1;

        case 'Fixed_base_PIendo_a_15_C1_1000'
            geometry='Ellipsoid';
            C1=1000; alpha=15;
            r_ff=0.55; r_fs=0.25; r_sn=0.2;
            pmax_endo=-1500;
            TOT_incrs=30;
            flag_prescr_base_disps=0;
            material_law='Guccione';
            fiber_endo_degrees=60; fiber_epi_degrees=-60;
            material_parameters=[C1,alpha*r_ff,alpha*r_sn,alpha*r_fs];
            % specify the flags for the Cheart file generation.
            flag_prescr_bound_base=0; %=1 if you prescr disps 
            flag_prescr_bound_epi=0;
            flag_prescr_bound_endo=0; 
            flag_prescr_bound_apex_hole=0;
            flag_prescr_APEX_NODE=0;
            flag_disc_Fib=0;
            flag_fix_base=1;
            flag_fix_apex_hole=0;    
            %% VF specs:
            flag_also_create_VF_epi=1; n_incr_VF=30;
            flag_continuous_fibers=1;

        case 'Sphere_Fixed_base_PIendo_a_9_C1_1000'
            gen_dir='C:\Users\an11\Desktop\Work\Synthetic_data\For_VF_as_unknown_BCs_solution_method_verification';
            geometry='Spheroid';
            Lagrange_meshes_interp_order='Quadratic';
            C1=1000; alpha=9;
            r_ff=0.333; r_fs=0.333; r_sn=0.333;
            pmax_endo=-1000;
            TOT_incrs=10;
            flag_prescr_base_disps=0;
            material_law='Guccione';
            fiber_endo_degrees=60; fiber_epi_degrees=-60;
            material_parameters=[C1,alpha*r_ff,alpha*r_sn,alpha*r_fs];
            % specify the flags for the Cheart file generation.
            flag_prescr_bound_base=0; %=1 if you prescr disps 
            flag_prescr_bound_epi=0;
            flag_prescr_bound_endo=0; 
            flag_prescr_bound_apex_hole=0;
            flag_prescr_APEX_NODE=0;
            flag_prescr_BASE_NODES=1;
            flag_disc_Fib=0;
            flag_fix_base=1;
            flag_fix_apex_hole=0;    
            mu_base=0;
            ref_frame_list=0;
            %% VF specs:
            flag_also_create_VF_epi=1; n_incr_VF=30;
            flag_continuous_fibers=1;
            sim_with_hole=1
    end
    if sum([sim_with_hole,sim_capped_apex,sim_collapsed])~=1
        disp('choose if sim has hole, capped apex or collapsed elements! - exiting script')
        return
    end
    
    general_Imaging_folder=[gen_dir,'/',geometry,'_Flat_base_hole'];
    general_Imaging_folder=[general_Imaging_folder,'/',Lagrange_meshes_interp_order,'_mesh'];
    if any([flag_my_fibers, flag_Daves_fibers, flag_fake_fibers])==0  || sum([flag_my_fibers, flag_Daves_fibers, flag_fake_fibers])>1
        warning('exactly one fiber fitting method must be chosen - set either "flag_my_fibers" or "flag_Daves_fibers" or "flag_fake_fibers" to 1')
                disp('stopping script.....')
                return
    else % My_Fibers Daves_Fibers 
        if flag_my_fibers==1
            general_Imaging_folder=[general_Imaging_folder,'/My_Fibers/Endo_',num2str(fiber_endo_degrees),'_epi_',num2str(abs(fiber_epi_degrees))];   

        elseif flag_Daves_fibers==1 % I am currently removing this option from the script as it is faster to run my code, works fine and I wanted to remove reduntant flags
            general_Imaging_folder=[general_Imaging_folder,'/Daves_Fibers/Endo_',num2str(fiber_endo_degrees),'_epi_',num2str(abs(fiber_epi_degrees))];          
        elseif flag_fake_fibers==1 % 
            general_Imaging_folder=[general_Imaging_folder,'/Fake_Fibers_XYZ'];
        end
    end   

    % fix boundaries to account for existence of septal elements - DO THIS
    % AFTER FITTING THE FIBERS! (or else your endo-epi laplace is messed up
    % (since only the non septal elements will be recognised as epicardium!)

    general_Imaging_folder=[general_Imaging_folder,'/',i_case];
    mesh_folder=[general_Imaging_folder,'/Mesh_per_frame_for_E_from_data'];
    press_folder=[general_Imaging_folder,'/pressure_data/default'];
    if flag_use_VF_just_base_and_hole_NeoHook==1
        VF_mesh_folder_path_intro=[general_Imaging_folder,'/Fixed_Base_and_ApexHole_NeoHookean'];
    elseif flag_use_VF_just_base_and_hole_ISO_Gucc==1
        VF_mesh_folder_path_intro=[general_Imaging_folder,'/Fixed_Base_and_ApexHole_ISOTROPIC_Guccione'];
    elseif flag_create_VF_epi_NeoHookean_sin_mu3==1
        VF_mesh_folder_path_intro=[general_Imaging_folder,'/Prescr_Tang_Epi_sinmu3_NeoHookean'];
    else
        VF_mesh_folder_path_intro=[];
    end
    
    
%     Lag_mesh_general_folder=[gen_folder,'/new_cases_hole/Cubic_Lagrange_Meshes'];
%             VF_Simulations_general_folder=[gen_folder,'/new_cases_hole/Cheart_simulations_for_VF_generation/Cubic_Lagrange_meshes'];
%             Data_analysis_general_folder=[gen_folder,'/new_cases_hole/Data_prep_for_param_estimation/Cubic_Lagrange_meshes'];
%             Simulations_general_folder=[gen_folder,'/new_cases_hole/Cheart_simulations_setup_folders/Cubic_Lagrange_meshes'];
%             
    switch Lagrange_meshes_interp_order
        case 'Cubic' % CUBIC LAGRANGE MESHES
            Nodes_per_elem_dir=4;
        case 'Quadratic' % QUADRATIC LAGRANGE MESHES
            Nodes_per_elem_dir=3;
    end
    GP_per_elem_dir=Nodes_per_elem_dir;
    
    
    %% continue with building the folder paths:
    % VF folder path:
%     if flag_use_VF_just_base_and_hole_NeoHook==1
%         VF_mesh_folder_path_intro=[VF_Simulations_general_folder,'/Fixed_Base_and_ApexHole_NeoHookean'];
%     elseif flag_use_VF_just_base_and_hole_ISO_Gucc==1
%         VF_mesh_folder_path_intro=[VF_Simulations_general_folder,'/Fixed_Base_and_ApexHole_ISOTROPIC_Guccione'];
%     elseif flag_create_VF_epi_NeoHookean_sin_mu3==1
%         VF_mesh_folder_path_intro=[VF_Simulations_general_folder,'/Prescr_Tang_Epi_sinmu3_NeoHookean'];
%     end
%     VF_mesh_folder_path_intro=[VF_mesh_folder_path_intro,'/case_',i_case];
%     
%     % Mesh and pressure data folder:
%     Mesh_and_pres_folder=[Data_analysis_general_folder,'/Cheart_Mesh_and_fitted_fibers_at_frame'];    
    
    
    
    for n_ref_fr=1:length(ref_frame_list)
        reference_frame=ref_frame_list(n_ref_fr);
        
%         Ref_mesh_folder=[Mesh_and_pres_folder,'/Ref_Frame_',num2str(reference_frame)];
        Ref_mesh_folder=[mesh_folder,'/Frame_',num2str(reference_frame)];
        
        press_data_file=[press_folder,'/pressure_at_frame.txt']; % ! note: this file is called frame_and_press.data in the clinical cases..
        data_pre=[];
        fid=fopen(press_data_file);
        tot_diast_frames=fscanf(fid,'%d',[1,1]); % number of frames in total
        data_pre=fscanf(fid,'%d     %f',[2,Inf]).';
        fclose(fid);
        indices_TU=find(data_pre(:,1)==reference_frame):length(data_pre(:,1));
        Diastolic_Frames_Indices_TU=data_pre(indices_TU,1); 
        diast_pressures_Pa_TU=data_pre(indices_TU,2);
        
%         fid=fopen([Ref_mesh_folder,'/Ref_mesh.T']);
        fid=fopen([Ref_mesh_folder,'/',Lagrange_meshes_interp_order,'_mesh_FE.T']);
        Element_intro=fscanf(fid,'%d',[2,1]).';
        Ref_mesh_Elements_ala_CHeart=fscanf(fid,'%d',[Nodes_per_elem_dir^3,Inf]).';
        fclose(fid);
        
        [Ref_mesh_Elements_ala_Cmgui,error]=turn_Cheart_node_ordering_to_Cmgui(Ref_mesh_Elements_ala_CHeart,Nodes_per_elem_dir);
        
%         fid=fopen([Ref_mesh_folder,'/Ref_mesh.X']);
        fid=fopen([Ref_mesh_folder,'/',Lagrange_meshes_interp_order,'_mesh_FE.X']);
        Nodes_intro=fscanf(fid,'%d',[2,1]).';
        Ref_mesh_Nodes=fscanf(fid,'%f',[3,Inf]).';
        fclose(fid);
        
%         fid=fopen([Ref_mesh_folder,'/Ref_mesh.B']);
        fid=fopen([Ref_mesh_folder,'/',Lagrange_meshes_interp_order,'_mesh_FE.B']);
        Bound_intro=fscanf(fid,'%d',[1,1]).';
        Ref_mesh_Boundaries=fscanf(fid,'%d',[Nodes_per_elem_dir^2+2,Inf]).';
        fclose(fid);
        
%         fid=fopen([Ref_mesh_folder,'/Ref_mesh.NBT']);
        fid=fopen([Ref_mesh_folder,'/',Lagrange_meshes_interp_order,'_mesh_FE.NBT']);
        Neigh_intro=fscanf(fid,'%d',[1,1]).';
        Ref_mesh_Neighbours=fscanf(fid,'%d',[6,Inf]).';
        fclose(fid);
                 
        if flag_disc_Fib==0
            fid=fopen([Ref_mesh_folder,'/FIBERS.X']);
            Fib_Nodes_intro=fscanf(fid,'%d',[2,1]).';
            Fib_mesh_Nodes=fscanf(fid,'%f',[9,Inf]).';
            fclose(fid);
        else
            disp('you need to build the discontinuous fiber field option- currently not available for meshes with hole');
        end
        

        for def_fr_ind=1:length(Diastolic_Frames_Indices_TU)
            frame_ID=Diastolic_Frames_Indices_TU(def_fr_ind);
            if flag_data==1
                Def_mesh_folder=[mesh_folder,'/Frame_',num2str(frame_ID)];
                % read Def_Nodes per frame in a structure
%                 fid=fopen([Ref_mesh_folder,'/Def_mesh-',num2str(frame_ID),'.X']);
                fid=fopen([Def_mesh_folder,'/',Lagrange_meshes_interp_order,'_mesh_FE.X']);
                Nodes_intro=fscanf(fid,'%d',[2,1]).';
                Def_mesh_Nodes=fscanf(fid,'%f',[3,Inf]).';
                fclose(fid);
                
                % and virtual field Disp Nodes per frame in a structure
                if isempty(VF_mesh_folder_path_intro)==0
                    list_DispVF=dir([VF_mesh_folder_path_intro,'/Frame_',num2str(frame_ID),'/Disp-*.D']);
                    end_VF_Disp_file=[VF_mesh_folder_path_intro,'/Frame_',num2str(frame_ID),'/Disp-',num2str(length(list_DispVF)-2),'.D']; %using the output from the last increment in the simulation as that has the largest deformation - if simulation didn't converge though maybe it's best to go 1or 2 increments before the last one to avoid weird deformation occuring.
                    disp('using 2 frames before last converged increment as VF disp to avoid possible issues with bad elements')
                end
                
            elseif flag_from_simulation==1
                Sim_prep_folder=[Simulations_mesh_folder_path,'/Ref_Frame_',num2str(reference_frame),'/Def_Frame_',num2str(frame_ID)]; 
                % check reference frame and pressure is according to
                % pressure data
                fid=fopen([Sim_prep_folder,'/',Lagrange_meshes_interp_order,'_mesh_FE.X']);
                Nodes_intro=fscanf(fid,'%d',[2,1]).';
                Lagrange_mesh_Nodes_Ref=fscanf(fid,'%f',[3,Inf]).';
                fclose(fid);   
                
                if max(max(abs(Lagrange_mesh_Nodes_Ref-Ref_mesh_Nodes)))>10^(-6)
                    disp('wtf?? ref nodes in data and simulation should match!!!???')
                    disp('');
                end
                [dp_dt,max_incrs,flag_found_it_press]=read_dp_dt_from_pressure_expression_in_Pfile(Sim_prep_folder,'Ellipse.P'); % dp_dt is <0 for inflation           
                % find max increments in simulation
                list_Disps_sim=dir([Sim_prep_folder,'/Disp-*.D']);
                if abs(-diast_pressures_Pa_TU(def_fr_ind)-dp_dt*max_incrs)>10^(-5) % diast_pressures_Pa_TU is >0 and hence add the minus sign (-) to give it the same sign as dp_dt.
                    disp('wtf?? pressure in data and simulation should match!!!???')
                    disp([' from data: ',num2str(-diast_pressures_Pa_TU(def_fr_ind)),' from sim: ',num2str(dp_dt*length(list_Disps_sim))])
                    disp('');
                else
                    diast_pressures_Pa_TU(def_fr_ind)=-dp_dt*length(list_Disps_sim); % correct for pressure as the maximum pressure may not have been reached (simulation not converged)-- and the (-dp_dt) is because diast_pressures_Pa_TU must be >0
                end
                % read that disp and add to reference frame for creating
                % deformed state
                fid=fopen([Sim_prep_folder,'/Disp-',num2str(length(list_Disps_sim)),'.D']);
                Nodes_intro=fscanf(fid,'%d',[2,1]).';
                Disp_Nodes_Ref=fscanf(fid,'%f',[3,Inf]).';
                fclose(fid); 
                Def_mesh_Nodes=Disp_Nodes_Ref+Ref_mesh_Nodes;
                
                % and virtual field Disp Nodes per frame in a structure
                if isempty(VF_mesh_folder_path_intro)==0
                    if flag_use_VF_just_base_and_hole_NeoHook==1
                        if exist([Sim_prep_folder,'/VF_Frame_',num2str(length(list_Disps_sim))],'dir')==1
                            VF_folder=[Sim_prep_folder,'/VF_Frame_',num2str(length(list_Disps_sim))]; % in this folder the last increment of the simulation with prescr disps is added to reference nodes and this deformed state is used as reference for generating the VF frame with fixed base&hole , so that the VF field satisfies the properties at the deformed configuration as is required                        
                        elseif exist([Sim_prep_folder,'/VFFixBasApNH_',num2str(length(list_Disps_sim))],'dir')==1
                            VF_folder=[Sim_prep_folder,'/VFFixBasApNH_',num2str(length(list_Disps_sim))];
                        else
                            disp('what is going on here? Is there no VF folder specified?');
                        end
                    elseif flag_create_VF_epi_NeoHookean_sin_mu3==1                  
                        VF_folder=[Sim_prep_folder,'/VFEpiNH_',num2str(length(list_Disps_sim))];
                    end                
                
                    list_DispVF=dir([VF_folder,'/Disp-*.D']);
                    end_VF_Disp_file=[VF_folder,'/Disp-',num2str(length(list_DispVF)-2),'.D']; %using the output from the last increment in the simulation as that has the largest deformation - if simulation didn't converge though maybe it's best to go 1or 2 increments before the last one to avoid weird deformation occuring.
                    disp('using 2 frames before last converged increment as VF disp to avoid possible issues with bad elements');
                end
            end
            Def_Nodes_per_Frame_struct.(['Frame_',num2str(frame_ID)])=Def_mesh_Nodes;
            
            % calculate cavity volume at deformed configuration
            [Def_cavity_volume, ~, ~]=calculate_cavity_volume_Lagrange_meshes_with_hole(Ref_mesh_Elements_ala_Cmgui,Def_mesh_Nodes,Ref_mesh_Boundaries,Endo_patch_ID,Epi_patch_ID,Base_patch_ID,Apex_hole_patch_ID,Nodes_per_elem_dir,GP_per_elem_dir,fig_ind);
            Def_cavity_volume_struct.(['Frame_',num2str(frame_ID)])=Def_cavity_volume;
            
             
            if isempty(VF_mesh_folder_path_intro)==0
                fid=fopen(end_VF_Disp_file);
                VFNodes_intro=fscanf(fid,'%d',[2,1]).';
                VF_Disp_Nodes=fscanf(fid,'%f',[3,Inf]).';
                fclose(fid);
                virt_disp_Nodes_per_Frame_struct.(['Frame_',num2str(frame_ID)])=VF_Disp_Nodes;       
            else
                virt_disp_Nodes_per_Frame_struct=[];
            end
        end
        [Ref_cavity_volume, ~, ~]=calculate_cavity_volume_Lagrange_meshes_with_hole(Ref_mesh_Elements_ala_Cmgui,Ref_mesh_Nodes,Ref_mesh_Boundaries,Endo_patch_ID,Epi_patch_ID,Base_patch_ID,Apex_hole_patch_ID,Nodes_per_elem_dir,GP_per_elem_dir,fig_ind);
         
        %% start insert function 18/4/2020
        [alpha_EB_CF_tu,alpha_EB_CF_pV,alpha_VWB_CF,C1_EB_CF_tu,C1_EB_CF_pV,C1_VWB_CF,EB_CF,EB_pV_CF,VWB_CF]=calculations_for_EB_and_VWB_functional(reference_frame,Diastolic_Frames_Indices_TU,diast_pressures_Pa_TU,Ref_mesh_Nodes,Ref_mesh_Elements_ala_Cmgui,Ref_mesh_Boundaries,Fib_mesh_Nodes,Ref_mesh_Elements_ala_Cmgui,alpha_PS_vector,r_ff,r_fs,r_sn,Def_Nodes_per_Frame_struct,virt_disp_Nodes_per_Frame_struct,Ref_cavity_volume,Def_cavity_volume_struct,Endo_patch_ID,Epi_patch_ID,flag_calc_virt_work_epi,GP_per_elem_dir,Nodes_per_elem_dir);

% % %         %% end insert function -18/4/2020
% % %         
% % %         % Wint_wrt_alpha: matrix of Wint per alpha in alpha_PS_vector per frame (size: [length(alpha_PS_vector),length(diastolic_frames_list)])
% % %         % virtual_Wint_wrt_alpha: matrix of \deltaWint per alpha in alpha_PS_vector per frame (size: [length(alpha_PS_vector),length(diastolic_frames_list)])
% % %         [Wint_wrt_alpha_per_frame, virtual_Wint_wrt_alpha_per_frame]=calc_Wint_and_virtWint_per_frame_and_a_for_CF_Guccione(alpha_PS_vector,r_ff,r_fs,r_sn,Diastolic_Frames_Indices_TU,Ref_mesh_Nodes,Ref_mesh_Elements_ala_Cmgui,Fib_mesh_Nodes,Ref_mesh_Elements_ala_Cmgui,Def_Nodes_per_Frame_struct,virt_disp_Nodes_per_Frame_struct,GP_per_elem_dir,Nodes_per_elem_dir);
% % %         % Wext_endo_per_DF: vector of Wext at endo per diastolic frame - column vector - size(length(diastolic_frames_list),1)
% % %         % virt_Wext_endo_per_DF: vector of virtual Wext at endo per diastolic frame - column vector - size(length(diastolic_frames_list),1)
% % % %         [Wext_endo_per_DF, virtual_Wext_endo_per_DF]=calc_Wext_endo_and_virtWext_endo_per_frame_for_CF(Diastolic_Frames_Indices_TU,diast_pressures_Pa_TU,Ref_mesh_Elements_ala_Cmgui,Ref_mesh_Boundaries,Endo_patch_ID,Def_Nodes_per_Frame_struct,virt_disp_Nodes_per_Frame_struct,GP_per_elem_dir,Nodes_per_elem_dir);
% % %         [Wext_endo_per_DF, virtual_Wext_endo_per_DF,virtual_Wext_epi_per_DF,virtual_Wext_epi_over_press_per_DF]=calc_Wext_endo_and_virtWext_endo_per_frame_for_CF_verbose_epi(Diastolic_Frames_Indices_TU,diast_pressures_Pa_TU,Ref_mesh_Elements_ala_Cmgui,Ref_mesh_Boundaries,Endo_patch_ID,Epi_patch_ID,Def_Nodes_per_Frame_struct,virt_disp_Nodes_per_Frame_struct,GP_per_elem_dir,Nodes_per_elem_dir,flag_calc_virt_work_epi);
% % %         for n_DF2=2:length(Diastolic_Frames_Indices_TU)
% % %             DF2=Diastolic_Frames_Indices_TU(n_DF2);
% % %             for n_DF1=1:(n_DF2-1)
% % %                 DF1=Diastolic_Frames_Indices_TU(n_DF1);
% % %                 
% % %                 %energy based CF:
% % %                 EBfunctional=abs(Wint_wrt_alpha_per_frame(:,n_DF1)./Wint_wrt_alpha_per_frame(:,n_DF2)-ones(size(Wint_wrt_alpha_per_frame,1),1)*Wext_endo_per_DF(n_DF1)/Wext_endo_per_DF(n_DF2));
% % %                 EB_CF.(['Ref_Fr_',num2str(reference_frame),'_DF1_',num2str(DF1),'_DF2_',num2str(DF2)])=EBfunctional;
% % %                 [~,min_ind_func]=min(EBfunctional);
% % %                 alpha_EB_CF.(['Ref_Fr_',num2str(reference_frame),'_DF1_',num2str(DF1),'_DF2_',num2str(DF2)])=alpha_PS_vector(min_ind_func);
% % %                 C1_EB_CF.(['Ref_Fr_',num2str(reference_frame),'_DF1_',num2str(DF1),'_DF2_',num2str(DF2)])=Wext_endo_per_DF(n_DF2)/Wint_wrt_alpha_per_frame(min_ind_func,n_DF2);
% % %                 
% % %                 %virtual works based CF:
% % %                 VWBfunctional=abs(virtual_Wint_wrt_alpha_per_frame(:,n_DF1)./virtual_Wint_wrt_alpha_per_frame(:,n_DF2)-ones(size(Wint_wrt_alpha_per_frame,1),1)*virtual_Wext_endo_per_DF(n_DF1)/virtual_Wext_endo_per_DF(n_DF2));
% % %                 VWB_CF.(['Ref_Fr_',num2str(reference_frame),'_DF1_',num2str(DF1),'_DF2_',num2str(DF2)])=VWBfunctional;
% % %                 [~,min_ind_func]=min(VWBfunctional);
% % %                 alpha_VWB_CF.(['Ref_Fr_',num2str(reference_frame),'_DF1_',num2str(DF1),'_DF2_',num2str(DF2)])=alpha_PS_vector(min_ind_func);
% % %                 C1_VWB_CF.(['Ref_Fr_',num2str(reference_frame),'_DF1_',num2str(DF1),'_DF2_',num2str(DF2)])=virtual_Wext_endo_per_DF(n_DF2)/virtual_Wint_wrt_alpha_per_frame(min_ind_func,n_DF2);
% % %                 
% % %             end            
% % %         end
    end
end
