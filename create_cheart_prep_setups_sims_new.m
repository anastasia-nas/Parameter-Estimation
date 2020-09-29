close all; clear all; clc;

% modified from create_cheart_prep_setups_new_cases_hole_new.m by Anastasia on 12-03-2020
%-------------------------------------------------------------------------------------------------------------------------
% Purpose: 
%-------------------------------------------------------------------------------------------------------------------------
% this script is to be run after
% My_cubic_Hermite_to_cubic_Lagrange_script_home.m (where you obtain a Cubic_mesh_FE.X and Cubic_mesh_FE.T file) to create a simulation
% ready folder for cheart. So a Linear mesh has to be created from the
% cubic one and .B and .NBT files have to be created for both Cubic and
% Linear meshes.
%-------------------------------------------------------------------------------------------------------------------------
% Notes on methodology:
%-------------------------------------------------------------------------------------------------------------------------
% MODIFIED from: create_cheart_prep_setups_new_cases_hole.m. Here I made
% changes so that this script can be used for two contradictory purposes in
% the way the folders are generated:
% 1) for running parameter sweeps with

% Cheart files are used as input, hence cubic mesh is easily turned into
% linear by just selecting the corner nodes and then renaming the nodes so
% they have a consistent numbering (from 1 to max_Node_ID in Lin_mesh).
% Boundary file (.B) names nodes in ksi1-ksi2 ordering.
% NBT file has ordering : [-ksi3,-ksi2,-ksi1,ksi1,ksi2,ksi3] % I got this
% from the prolate spheroidal meshes - I hope it's correct..
%-------------------------------------------------------------------------------------------------------------------------
% Input:
%-------------------------------------------------------------------------------------------------------------------------
% Cub_mesh_folder: folder where the Cubic_mesh_FE.X and Cubic_mesh_FE.T files are stored.
% Apex_hole_patch_ID: leave empty if no apical hole in mesh, otherwise =4


%-------------------------------------------------------------------------------------------------------------------------
% Output:
%-------------------------------------------------------------------------------------------------------------------------
% printed Linear_mesh_FE.X Linear_mesh_FE.T files and Linear_mesh_FE.B and
% Linear_mesh_FE.NBT and Cubic_mesh_FE.B and Cubic_mesh_FE.NBT files
%-------------------------------------------------------------------------------------------------------------------------
% SOS points:
%-------------------------------------------------------------------------------------------------------------------------
% 1) NBT file has this ordering: [-ksi3,-ksi2,-ksi1,ksi1,ksi2,ksi3]
% 2) I use the corrected pressure values ('PAS-NF') version in parameter
% estimation which means I only use the frames from chosen reference frames till the end diastolic one (as identified from the p-V analysis)
% and from these I only use those where pressure is positive. Then I offset the pressure trace in those frames by the pressure at the
% first of these frames. Hence choosing a different reference frame means
% that the pressure data are modified.

% use then script to turn fibers into continuous mesh

%================================================================================================================================
flag_Joshuas_laptop=0; % Joshua leave this to =1..

i_case='OSU';
% FOLDERS:
if flag_Joshuas_laptop==1
    folder_pre_fix='C:\Users\fitzp\Documents\Josh_Documents\KCL\iBSc\iBSc_4_Stiffness_Project';
else
    folder_pre_fix='C:\Users\an11\Desktop\Work\Teaching\Stiffness_Project_Joshua_Kader';
end
gen_folder_new_cases=[folder_pre_fix,'\Manav_Cases_Meshed_341']; % use your _341 or _9124 folders here
pres_frame_data_folder_intro_new_cases=[folder_pre_fix,'\new_cases_pressures_JK'];% from laptop

% mesh/data
Lagrange_meshes_interp_order='Cubic' ;  % choose between 'Cubic' and 'Quadratic'
flag_rescale_mesh_to_meters=1; % if you use the meshes converted to m
flag_correct_pressure_data=1; % leave like that here
% FLAGS:

%% prep for data analysis for parameter estimation (print boundaries, neighbours and used pressure per frame after correction according to pipeline (where min pressure is at reference frame).
flag_prep_for_parameter_estimation=1; % this option will override all others 
%% for simulations with prescribed disps - for anatomical synthetic data generation / parameter sweeps
flag_prescr_basal_disps=0; % =1 if you prescribe basal disps, =0 otherwise. This option is for creating template folders for Cheart sims with prescribed base for parameter sweeps or for synthetic data calcuation or for verifying the L2 norm between sims with identified parameters and data.
flag_prescr_apicalHole_disps=0; % =1 for more realism to data - however this might lead to strain concentrations around the hole that may affect results..
flag_fix_apex=0; % only for synthetic data with anatomical meshes. this choice will help identify what is causing the discordance in the analysis. If it is the deformation at the hole that leads to divergence.
flag_press_epi=0; % =1 if you also apply pressure to epi for simulating effect of pericardium. (pressure is chosen to be 1/3 of endocardial pressure)
%% choose fiber fitting options if your material is anisotropic. Option flag_fake_fibers is to be used ONLY for VF generation
flag_my_fibers=1; flag_Daves_fibers=0; flag_fake_fibers=0;% flag_my_fibers=1 for using my tool-flag_Daves_fibers=1 for preping folder to run cheart fibre fitting tool, or both can be zero if you have an isotropic material
%% choose these for generating VF: (note flag_prescr_basal_disps must be ==1!)
%note if mesh has hole you also need to fix hole to avoid having work doneat the hole by the stresses:
flag_create_VF_just_base_and_hole_NeoHook=0; % use this for tackling unknown base tractions with VF! -- =1 if this script is run to prepare VF generation with fixed base and neohookean material (for any other type of VF generating simulation - e.g VF for epi/Stokes/isotropic Guccione/ whatever- create a new flag describing the case!)- this didn't work well in Cheart (sth wrong with NeoHookean?)-no I think it was the hole (that was not fixed and deforming, and also the steps were too big (too large pressure increments)
flag_create_VF_just_base_and_hole_ISO_Gucc=0;% this is slower and leads to reduced deformation  -- =1 if this script is run to prepare VF generation with fixed base and isotropic Guccione material (for any other type of VF generating simulation - e.g VF for epi/Stokes/whatever- create a new flag describing the case!)
flag_create_VF_epi_NeoHookean_sin_mu3=0; % prescribed tangential disp at epi maximum at midwall and zero at apex and base this is achieved by weighting the unit tangential vectors by the sin(mu^3) where mu is the scalar field obtained by solving the Dirichlet Laplace apex to base.
% SPECS:
fig_ind=[]; % for numbering plots
Endo_patch_ID=1;
Epi_patch_ID=2;
Base_patch_ID=3;
Apex_hole_patch_ID=4; % =[] if no apical hole exists in mesh.
Septum_patch_ID=[]; % =5 if there is hole, 4 otherwise.
simulation_increments=30;
pmax=-150; % <0 for passive inflation, >0 for deflation -- this is to be used only in case of building the VF
% existence of epi/septum (also set Septum_patch_ID to =5)
press_expr_sept=[]; press_expr_epi=[];
% fibers:
flag_disc_Fib=0; % print continuous fiber field - no need to waste space with discontinuous -- use script turn_discontinuous_fibers_to_continuous_new_cases_hole.m for the conversion from cheart fiber fitting tool.
fiber_epi_degrees=-60;fiber_endo_degrees=60; % fiber angle at epi must be negative!
% material parameters default:
r_ff=0.55; r_fs=0.25; r_sn=0.2; % change them further down if you need an isotropic material for example


flag_VF_sim_only__fit_fibers_to_def_configuration=0; % initialisation - will be turned to =1 if one of the VF building options has been chosen

%=============================================================================================================================================================
%% catching quick errors about contradicting flag choices::
if sum([flag_prep_for_parameter_estimation,flag_prescr_basal_disps,flag_create_VF_just_base_and_hole_NeoHook,flag_create_VF_just_base_and_hole_ISO_Gucc,flag_create_VF_epi_NeoHookean_sin_mu3])>1
   warning('ERROR: these flags are mutually exclusive as they correspond to completely different actions. Set one to 1 and the rest to 0 to proceed!')
   warning('Stopping script "create_cheart_prep_setups_new_cases_hole_new" execution now!')
   return
end

if sum([flag_prescr_apicalHole_disps,flag_fix_apex])>1 
    warning('ERROR: these flags are mutually exclusive - you either fix or prescribe disps at apical hole');
    return
end
if flag_create_VF_just_base_and_hole_NeoHook==1 || flag_create_VF_just_base_and_hole_ISO_Gucc==1 || flag_create_VF_epi_NeoHookean_sin_mu3==1 % || any other VF field generation
    disp('this is run to generate VF  fields- so you only care about generating simulations at each frame without caring about reference frames')
    disp('fibers_must_be_fitted_to_deformed_config');
    flag_VF_sim_only__fit_fibers_to_def_configuration=1;
    flag_prescr_basal_disps=0;     
    disp('forcing prescribed BCs to be fixed base!');  
    flag_prescr_apicalHole_disps=0; 
    flag_fix_apex=0; 
    disp('returning all other displacement prescription flags to default (=0)');
end
% safety check return epi pressure flag to zero:
if (flag_create_VF_just_base_and_hole_NeoHook==1 || flag_create_VF_just_base_and_hole_ISO_Gucc==1 || flag_create_VF_epi_NeoHookean_sin_mu3==1 ) || flag_prep_for_parameter_estimation==1
    flag_press_epi=0; % no need  for prescribing pressure at epi here, this is only for the creation of anatomical in silico cases...
end
if fiber_endo_degrees<0 || fiber_epi_degrees>0
    disp('ERROR! fiber angle at epi must negative and endo positive - see Streeter 1960 - exiting script')
    return
end

% for new cases:
if any(strcmp(i_case,{'BAL','HAW','JER_new','JUL_new','MCP','OSU','ROL','SAN','WIL'}))==1
    switch Lagrange_meshes_interp_order
        case 'Cubic' % CUBIC LAGRANGE MESHES
            Lag_mesh_general_folder=[gen_folder_new_cases,'/Cubic_Lagrange_Meshes'];
            if flag_rescale_mesh_to_meters==1
                Lag_mesh_general_folder=[Lag_mesh_general_folder,'_in_m'];
            end
            Simulations_general_folder=[gen_folder_new_cases,'/Cheart_simulations_setup_folders/Cubic_Lagrange_meshes'];
            VF_Simulations_general_folder=[gen_folder_new_cases,'/Cheart_simulations_for_VF_generation/Cubic_Lagrange_meshes'];
            Data_analysis_general_folder=[gen_folder_new_cases,'/Data_prep_for_param_estimation/Cubic_Lagrange_meshes'];
            Nodes_per_elem_dir=4;

        case 'Quadratic' % QUADRATIC LAGRANGE MESHES
            Lag_mesh_general_folder=[gen_folder_new_cases,'/Quadratic_Lagrange_Meshes'];
            if flag_rescale_mesh_to_meters==1
                Lag_mesh_general_folder=[Lag_mesh_general_folder,'_in_m'];
            end
            Simulations_general_folder=[gen_folder_new_cases,'/Cheart_simulations_setup_folders/Quadratic_Lagrange_meshes'];
            VF_Simulations_general_folder=[gen_folder_new_cases,'/Cheart_simulations_for_VF_generation/Quadratic_Lagrange_meshes'];
            Data_analysis_general_folder=[gen_folder_new_cases,'/Data_prep_for_param_estimation/Quadratic_Lagrange_meshes'];
            Nodes_per_elem_dir=3;
    end

    GP_per_elem_dir=Nodes_per_elem_dir;

    Lag_mesh_folder_path=[Lag_mesh_general_folder,'/case_',i_case]; % folder where Lagrange meshes are stored after processing from Cubic Hermite meshes.
   
    % for more choices copy them from Refine_Cub_Lagrange_mesh_for_Cheart_new_cases_hole.m
    
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
            ref_frame_list=[15:17]; % default reference frame: 16
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
            ref_frame_list=27;% [26:28]; % default reference frame is 27 
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
    
  
    
    %% if this about simulations for creating fake data (anisotropic material and prescribing BCs from data
    if flag_create_VF_just_base_and_hole_NeoHook==0 && flag_create_VF_just_base_and_hole_ISO_Gucc==0 && flag_create_VF_epi_NeoHookean_sin_mu3==0 && flag_prep_for_parameter_estimation==0 % && and all other VF flags zero: ( and all other flags regarding VF generating options ( so not real data processing) are set to =0 )
        C2=r_ff*alpha; C4=r_fs*alpha; C3=r_sn*alpha; % the order these are written in .P file
        material_law='Guccione'; material_parameters=[C1,C2,C3,C4];
        if flag_prescr_basal_disps==1  
            if flag_prescr_apicalHole_disps==1
                Simulations_mesh_folder_path=[Simulations_general_folder,'/Prescribe_Disps_at_Base_Apex_Nodes'];
            elseif flag_fix_apex==1
                Simulations_mesh_folder_path=[Simulations_general_folder,'/Prescribe_Disps_at_Base_Nodes_Fix_Apex'];
            else
                Simulations_mesh_folder_path=[Simulations_general_folder,'/Prescribe_Disps_at_Base_Nodes'];            
            end
            % .P file generation specs: 
            flag_base_P_file=1; flag_epi_P_file=0; flag_endo_P_file=0; flag_apex_P_file=0; flag_prescr_APEX_NODE_P_file=0; flag_fix_base_P_file=0; flag_fix_apex_hole_P_file=0;
        % else if % add flags if you want to prescribe alternative BCs
        else % this corresponds to preparing Cheart sims for running like in " flag_prescr_basal_disps=1 " but with fixed base (no prescribed basal disps)
            if flag_prescr_apicalHole_disps==1
                Simulations_mesh_folder_path=[Simulations_general_folder,'/Fixed_Base_Prescr_Apex']; % this combination is not that useful as a fixed base is usually chosen for control and prescribing apex is not that important as is prescribing base, so no real point in trying this other than for numerical reasons (checking which boundary causes simulation to non-coverge)
            elseif flag_fix_apex==1
                Simulations_mesh_folder_path=[Simulations_general_folder,'/Fixed_Base_Apex'];
            else
                Simulations_mesh_folder_path=[Simulations_general_folder,'/Fixed_Base'];
            end
            % .P file generation specs: 
            flag_base_P_file=0; flag_epi_P_file=0; flag_endo_P_file=0; flag_apex_P_file=0; flag_prescr_APEX_NODE_P_file=0; flag_fix_base_P_file=1; flag_fix_apex_hole_P_file=0;
        end
        if flag_press_epi==1 % this flag only applies to the simulations for creating fake data
            Simulations_mesh_folder_path=[Simulations_mesh_folder_path,'_PresEpi'];
        end
    else % add each VF generating type here:
        if flag_create_VF_just_base_and_hole_NeoHook==1
            simulation_increments=30; pmax=-150; % I chose this after comparing options for BAL Frame 24- see Zim notes for more. (pmax<0 for passive inflation, >0 for deflation )
            disp('')
            disp('SOS: overwriting pmax and simulations increments for the case of VF with fixed base & apex and Neohookean material')
            Simulations_mesh_folder_path=[VF_Simulations_general_folder,'/Fixed_Base_and_ApexHole_NeoHookean'];
            flag_Daves_fibers=0; flag_my_fibers=0; flag_fake_fibers=0; % no fibers required here.
            % .P file generation specs:
            material_law='NeoHookean'; material_parameters=1000; dp_dt=pmax/simulation_increments;
            flag_base_P_file=0; flag_epi_P_file=0; flag_endo_P_file=0; flag_apex_P_file=0; flag_prescr_APEX_NODE_P_file=0;  flag_fix_base_P_file=1; flag_fix_apex_hole_P_file=1;
        end  
        if flag_create_VF_just_base_and_hole_ISO_Gucc==1
            Simulations_mesh_folder_path=[VF_Simulations_general_folder,'/Fixed_Base_and_ApexHole_ISOTROPIC_Guccione'];
            
            % .P file generation specs:
            material_law='Guccione'; material_parameters=[1000,5,5,5]; dp_dt=pmax/simulation_increments;
            flag_base_P_file=0; flag_epi_P_file=0; flag_endo_P_file=0; flag_apex_P_file=0; flag_prescr_APEX_NODE_P_file=0;  flag_fix_base_P_file=1; flag_fix_apex_hole_P_file=1;
        end  
        if flag_create_VF_epi_NeoHookean_sin_mu3==1
            simulation_increments=30; pmax=-150; % I chose this after comparing options for BAL Frame 24- see Zim notes for more. (pmax<0 for passive inflation, >0 for deflation )
            max_disp=0.01; % max disp tangentially at epi occuring midwall
            disp('')
            disp('SOS: overwriting pmax and simulations increments for the case of VF with prescribed tangential disp at epi and fixed base & apex and Neohookean material')
            Simulations_mesh_folder_path=[VF_Simulations_general_folder,'/Prescr_Tang_Epi_sinmu3_NeoHookean'];
            flag_Daves_fibers=0; flag_my_fibers=0; flag_fake_fibers=0; % no fibers required here.
            % .P file generation specs:
            material_law='NeoHookean'; material_parameters=100; dp_dt=pmax/simulation_increments;
            flag_base_P_file=0; flag_epi_P_file=1; flag_endo_P_file=0; flag_apex_P_file=0; flag_prescr_APEX_NODE_P_file=0;  flag_fix_base_P_file=1; flag_fix_apex_hole_P_file=1;            
        end
        
        if flag_prep_for_parameter_estimation==1 % if you just prep folders for parameter estimation processing.
            % copy structure from
            % create_copies_of_sim_outputs_to_Mesh_per_frame_for_E_from_data3
            % (where pressure and deformation data are separated) to allow
            % to insert noise.
            
            Simulations_mesh_folder_path=[Data_analysis_general_folder,'/Cheart_Mesh_and_fitted_fibers_at_frame'];
%             Pressure_data_folder_path=[Data_analysis_general_folder,'/Pressure_data_per_frame'];
%             % won't use a separate folder for this will just store the
%             pressure under the reference frame, since they are reference
%             frame dependent (following the processing)
%             if flag_correct_pressure_data==1
%                 Pressure_data_folder_path=[Pressure_data_folder_path,'/Corrected_pressures_no_neg'];
%             end
        end
        
    end
    
    % check a fiber fitting option is selected if an isotropic material law
    % is chosen:
    if flag_prep_for_parameter_estimation==0 % I don't need I material law description other wise
        if strcmp(material_law,'Guccione') % !add any other anisotropic laws here!!!
            if any([flag_my_fibers, flag_Daves_fibers, flag_fake_fibers])==0 % the method for fitting fibers must be specified
                warning('no fiber fitting method specified in " " - set either "flag_my_fibers" or "flag_Daves_fibers" to 1')
                disp('stopping script.....')
                return
            end
            if sum([flag_my_fibers, flag_Daves_fibers, flag_fake_fibers])>1 % don't use contradictory flags for fitting fibers
                warning('more than one fiber fitting method found in " " - set either "flag_my_fibers" or "flag_Daves_fibers" or "flag_fake_fibers" to 1')
                disp('stopping script.....')
                return
            end
        end
    end
    % now specify the fiber fitting if any:
    if flag_my_fibers==1 || flag_Daves_fibers==1 || flag_fake_fibers==1
        if flag_my_fibers==1
            Simulations_mesh_folder_path=[Simulations_mesh_folder_path,'/My_Fibers/Endo_',num2str(fiber_endo_degrees),'_epi_',num2str(abs(fiber_epi_degrees))]; %folder paths are too close to max path length - economize!
            disp('describe how you will fit fibers here!!!!!!!!!!!')
        elseif flag_Daves_fibers==1 % I am currently removing this option from the script as it is faster to run my code, works fine and I wanted to remove reduntant flags
            Simulations_mesh_folder_path=[Simulations_mesh_folder_path,'/Daves_Fibers/Endo_',num2str(fiber_endo_degrees),'_epi_',num2str(abs(fiber_epi_degrees))];
            
        elseif flag_fake_fibers==1 % 
            Simulations_mesh_folder_path=[Simulations_mesh_folder_path,'/Fake_Fibers_XYZ'];


        end
    end
    Simulations_mesh_folder_path=[Simulations_mesh_folder_path,'/case_',i_case];
    mkdir(Simulations_mesh_folder_path);
    
    
    % pressure data folder:
    diast_frames_info_folder=[pres_frame_data_folder_intro_new_cases,'/',i_case_pres,'/',press_data_folder,'/process_data/sync_pV'];
    if isempty(regexp(i_case,'_new','match'))==0 % an exei _new to case (JUL_new,JER_new)
        diast_frames_info_folder=[pres_frame_data_folder_intro_new_cases,'/',i_case_pres,'/',press_data_folder,'/process_data_new_mesh/sync_pV'];
    end
    if isempty(regexp(i_case,'_iso','match'))==0
        diast_frames_info_folder=[pres_frame_data_folder_intro_new_cases,'/',i_case_pres,'/',press_data_folder,'/process_data_iso_meshes/sync_pV'];
    end
    
    if flag_VF_sim_only__fit_fibers_to_def_configuration==1
        % don't care about reference frames. Each simulation will be
        % generated by treating the deformed config as a reference frame.
        % So I only use 1 frame to read the Elements and Nodes and create
        % Boundaries, Neighbours and Lin_meshes for pressure (as is the case for the "proper processing" for fake data or parameter sweep simulation.
        ref_frame_list_TU=ref_frame_list(1);        
    else
        ref_frame_list_TU=ref_frame_list;
    end
    for n_ref_fr=1:length(ref_frame_list_TU)
        reference_frame=ref_frame_list_TU(n_ref_fr);
        
        fid=fopen([Lag_mesh_folder_path,'/Frame_',num2str(reference_frame),'/',Lagrange_meshes_interp_order,'_mesh_FE.X']);
        Nodes_intro=fscanf(fid,'%d',[2,1]).';
        Lagrange_mesh_Nodes_Ref=fscanf(fid,'%f',[3,Inf]).';
        fclose(fid);    
                                
        if n_ref_fr==1   % you only need to do this once if your element and node numbering is consistent across frames -check                                   
            fid=fopen([Lag_mesh_folder_path,'/Frame_',num2str(reference_frame),'/',Lagrange_meshes_interp_order,'_mesh_FE.T']);
            Element_intro=fscanf(fid,'%d',[2,1]).';
            Lagrange_mesh_elements_ala_CHeart_Ref=fscanf(fid,'%f',[Nodes_per_elem_dir^3,Inf]).';
            fclose(fid);

            [Lagrange_mesh_Elements_Cmgui,error]=turn_Cheart_node_ordering_to_Cmgui(Lagrange_mesh_elements_ala_CHeart_Ref,Nodes_per_elem_dir);
                            
            % Find cubic mesh Neighbours & create .NBT file:
            Lagrange_mesh_Neighbours=find_Neighbours_for_NBT_from_Lagrange_mesh_cmgui_ordering(Lagrange_mesh_Elements_Cmgui,Nodes_per_elem_dir);
            % Find cubic mesh Boundaries & create .B file:
%           old and wrong version of finding boundaries:  Cubic_mesh_Boundaries=find_Boundaries_for_Bfile_from_Lagrange_mesh_cmgui_ordering_new(Cubic_mesh_Nodes_Ref,Cubic_mesh_Elements_Cmgui,Cubic_mesh_Neighbours,Nodes_per_elem_dir,Endo_patch_ID,Epi_patch_ID,Base_patch_ID,Apex_hole_patch_ID);
            Lagrange_mesh_Boundaries=find_Boundaries_for_Bfile_from_Lagrange_mesh_cmgui_order_struct(Lagrange_mesh_Nodes_Ref,Lagrange_mesh_Elements_Cmgui,Lagrange_mesh_Neighbours,Nodes_per_elem_dir,Endo_patch_ID,Epi_patch_ID,Base_patch_ID,Apex_hole_patch_ID);
%             %for debugging find_Boundaries_for_Bfile_from_Lagrange_mesh_cmgui_ordering_new
%             Neighbours=Lagrange_mesh_Neighbours; Lagrange_Elements_CMGUI_ordering=Lagrange_mesh_Elements_Cmgui; Lagrange_Nodes=Lagrange_mesh_Nodes_Ref;
           
            % find endocardial epicardial basal and hole elements and
            % nodes:
            endo_elem_indices_in_Bound_mat=find(Lagrange_mesh_Boundaries(:,Nodes_per_elem_dir^2+2)==Endo_patch_ID);
            Endocardial_Elements_IDs=Lagrange_mesh_Boundaries(endo_elem_indices_in_Bound_mat,1);
            Endocardial_Nodes_IDs=sort(unique(reshape(Lagrange_mesh_Boundaries(endo_elem_indices_in_Bound_mat,2:Nodes_per_elem_dir^2+1),1,[]))); % turn endocardial boundary nodes into a row and then make it unique and sort it
            
            epi_elem_indices_in_Bound_mat=find(Lagrange_mesh_Boundaries(:,Nodes_per_elem_dir^2+2)==Epi_patch_ID);
            Epicardial_Elements_IDs=Lagrange_mesh_Boundaries(epi_elem_indices_in_Bound_mat,1);
            Epicardial_Nodes_IDs=sort(unique(reshape(Lagrange_mesh_Boundaries(epi_elem_indices_in_Bound_mat,2:Nodes_per_elem_dir^2+1),1,[]))); % turn epicardial boundary nodes into a row and then make it unique and sort it
            
            base_elem_indices_in_Bound_mat=find(Lagrange_mesh_Boundaries(:,Nodes_per_elem_dir^2+2)==Base_patch_ID);
            Basal_Elements_IDs=Lagrange_mesh_Boundaries(base_elem_indices_in_Bound_mat,1);
            Basal_Nodes_IDs=sort(unique(reshape(Lagrange_mesh_Boundaries(base_elem_indices_in_Bound_mat,2:Nodes_per_elem_dir^2+1),1,[]))); % turn basal boundary nodes into a row and then make it unique and sort it
                   
            if isempty(Apex_hole_patch_ID)==0
                apex_hole_elem_indices_in_Bound_mat=find(Lagrange_mesh_Boundaries(:,Nodes_per_elem_dir^2+2)==Apex_hole_patch_ID);
                Apical_Hole_Elements_IDs=Lagrange_mesh_Boundaries(apex_hole_elem_indices_in_Bound_mat,1);
                Apical_Hole_Nodes_IDs=sort(unique(reshape(Lagrange_mesh_Boundaries(apex_hole_elem_indices_in_Bound_mat,2:Nodes_per_elem_dir^2+1),1,[]))); % turn basal boundary nodes into a row and then make it unique and sort it            
            end
        
        end                
                
        % Create Linear Mesh corresponding to Cubic mesh at Reference frame:
        [Linear_mesh_Nodes,Linear_mesh_Elements,Linear_mesh_Boundaries]=create_Lin_mesh_from_Cubic_mesh_cheart_ordering(Lagrange_mesh_Nodes_Ref,Lagrange_mesh_elements_ala_CHeart_Ref,Lagrange_mesh_Boundaries,Nodes_per_elem_dir);
%         [Linear_mesh_Nodes,Linear_mesh_Elements,Linear_mesh_Boundaries]=create_Lin_mesh_from_Cubic_mesh_cheart_ordering(Cubic_mesh_nodes,Cubic_mesh_elements_ala_CHeart,Cubic_mesh_Boundaries,Nodes_per_elem_dir)


        if flag_correct_pressure_data==1
            [diast_pressures_Pa,Diastolic_Frames_Indices,Diastolic_Frames_Volumes_ml,flag_reference_frame_negative_press]=read_corrected_pressures_from_catheter_data_new_cases(diast_frames_info_folder,chosen_offset,reference_frame);
            if flag_reference_frame_negative_press==1
                disp('problemo!!!! reference frame has p<0 and has been removed from data -  change things..')
            end
            if flag_VF_sim_only__fit_fibers_to_def_configuration==1
                Diastolic_Frames_Indices_TU=Diastolic_Frames_Indices; % do this for all frames of interest (in previous frames did it for all frames but no point since I won't use them.
            else
                Diastolic_Frames_Indices_TU=[];
                ref_fr_ind=find(Diastolic_Frames_Indices==reference_frame);
                Diastolic_Frames_Indices_TU=Diastolic_Frames_Indices(ref_fr_ind+1:end); % removing frames previous to (and including) the reference frame
                diast_pressures_Pa_TU=diast_pressures_Pa(ref_fr_ind+1:end);
                Diastolic_Frames_Volumes_ml_TU=Diastolic_Frames_Volumes_ml(ref_fr_ind+1:end);  
            end
        else
            disp('no alternative yet!--think this through and remember I used this correction for BMMB!')
        end
            
        if flag_VF_sim_only__fit_fibers_to_def_configuration==1
%             Diastolic_Frames_Indices_TU=Diastolic_Frames_Indices; % do this for all frames of interest (in previous frames did it for all frames but no point since I won't use them.
        else %This is a function of the reference frame so you need to re-evaluate this per reference frame.
            
            
            if flag_prep_for_parameter_estimation==1                
                if flag_correct_pressure_data==1
                    Pressure_data_folder_path=[Simulations_mesh_folder_path,'/Ref_Frame_',num2str(reference_frame),'/Corrected_pressures_no_neg'];                    
                end
                mkdir(Pressure_data_folder_path);
                press_data_file=[Pressure_data_folder_path,'/frame_and_press.data'];
                fid=fopen(press_data_file,'w');
                fprintf(fid,'%d\n',length(Diastolic_Frames_Indices_TU)); % number of frames in total
                fprintf(fid,'%d     %21.18f\n',[Diastolic_Frames_Indices_TU.',diast_pressures_Pa_TU.'].'); % Diastolic_Frames_Indices_TU, diast_pressures_Pa_TU are row vectors
                fclose(fid);
                % print the above data..
                %think how it's best to output these (especially
                %diast_pressures_Pa_TU,and Diastolic_Frames_Indices_TU) -
                %the volumes were calculated for original coarser and with
                %collapsed element meshes.
                
            end
        end
        
        if flag_my_fibers==1 && flag_VF_sim_only__fit_fibers_to_def_configuration==0 % if flag_VF_sim_only__fit_fibers_to_def_configuration==1 then I need to fit fibers to every deformed config
            Endo_patch_ID=Endo_patch_ID; Epi_patch_ID=Epi_patch_ID; Base_patch_ID=Base_patch_ID; Apex_hole_patch_ID=Apex_hole_patch_ID; fig_ind=[];
            if flag_disc_Fib==0
                [Fibers_nodes,~,~,~]=fit_fibers_using_Laplace_Dirichlet_solver_no_imbric_angle(fiber_epi_degrees,fiber_endo_degrees,flag_disc_Fib,Lagrange_mesh_Elements_Cmgui,Lagrange_mesh_Nodes_Ref,Lagrange_mesh_Boundaries,Nodes_per_elem_dir,GP_per_elem_dir,Endo_patch_ID,Epi_patch_ID,Base_patch_ID,Apex_hole_patch_ID,fig_ind);
            end
        end
        if flag_fake_fibers==1 % then the frame to which the fibers are fitted is no of no importance because they are constantly parallel to x,y,z
            if flag_disc_Fib==0
                Fibers_nodes=[ones(size(Lagrange_mesh_Nodes_Ref,1),1),zeros(size(Lagrange_mesh_Nodes_Ref,1),1),zeros(size(Lagrange_mesh_Nodes_Ref,1),1),...
                              zeros(size(Lagrange_mesh_Nodes_Ref,1),1),ones(size(Lagrange_mesh_Nodes_Ref,1),1),zeros(size(Lagrange_mesh_Nodes_Ref,1),1),...
                              zeros(size(Lagrange_mesh_Nodes_Ref,1),1),zeros(size(Lagrange_mesh_Nodes_Ref,1),1),ones(size(Lagrange_mesh_Nodes_Ref,1),1)];
            end
        end
        
        if flag_prep_for_parameter_estimation==1
            % create a reference folder with the Reference frame mesh and fibers and Def_nodes at the deformed frames 
            % you don't need anything else and avoid extra folders
            Sim_prep_folder=[Simulations_mesh_folder_path,'/Ref_Frame_',num2str(reference_frame)];
            mkdir(Sim_prep_folder);
                    
            copyfile([Lag_mesh_folder_path,'/Frame_',num2str(reference_frame),'/',Lagrange_meshes_interp_order,'_mesh_FE.X'],[Sim_prep_folder,'/Ref_mesh.X']);
            % Cubic_mesh_FE.T
            copyfile([Lag_mesh_folder_path,'/Frame_',num2str(reference_frame),'/',Lagrange_meshes_interp_order,'_mesh_FE.T'],[Sim_prep_folder,'/Ref_mesh.T']);
            % Cubic_mesh_FE.B
            B_file_name='Ref_mesh.B';
            create_Cheart_Bfile_from_Boundaries_mat(Lagrange_mesh_Boundaries,Nodes_per_elem_dir,Sim_prep_folder,B_file_name);
            % Cubic_mesh_FE.NBT            
            NBT_file_name='Ref_mesh.NBT';
            create_Cheart_NBTfile_from_Neighbours_mat(Lagrange_mesh_Neighbours,Sim_prep_folder,NBT_file_name);
            
            if flag_disc_Fib==0
                fid=fopen([Sim_prep_folder,'/FIBERS.X'],'w');
                fprintf(fid,'%d %d\n',[size(Fibers_nodes,1),size(Fibers_nodes,2)].');
                fprintf(fid,'%21.15f    %21.15f    %21.15f    %21.15f    %21.15f    %21.15f    %21.15f    %21.15f    %21.15f\n',Fibers_nodes.');
                fclose(fid);                    
            end
            
        end
            
        for def_fr_ind=1:length(Diastolic_Frames_Indices_TU)
            
            frame_ID=Diastolic_Frames_Indices_TU(def_fr_ind);
            if flag_VF_sim_only__fit_fibers_to_def_configuration==1 % reference frame is of no consequence, each deformed frame serves as reference frame
                Sim_prep_folder=[Simulations_mesh_folder_path,'/Frame_',num2str(frame_ID)]; 
            else % here the reference frame matters
                if flag_prep_for_parameter_estimation==1 % in this case you don't need a Def_Frame subfolder store everything you need under the 'Ref_Frame' folder.                    
                    % Sim_prep_folder doesn't change per frame (only a
                    % folder per reference frame is needed here and this
                    % was already created.                    
                else 
                    Sim_prep_folder=[Simulations_mesh_folder_path,'/Ref_Frame_',num2str(reference_frame),'/Def_Frame_',num2str(frame_ID)]; 
                end
            end
            if exist(Sim_prep_folder,'dir')==0 % in the case of flag_prep_for_parameter_estimation=1 then only a folder for reference frame will be created.
                mkdir(Sim_prep_folder);                    
            end
            
            
            %% copy and print the required Cubic mesh and Lin mesh files
            if flag_prep_for_parameter_estimation==0 
                if flag_VF_sim_only__fit_fibers_to_def_configuration==1 % if the reference frame is the deformed configuration (at current frame)
                    
%                   % Cubic_mesh_FE.X
                    copyfile([Lag_mesh_folder_path,'/Frame_',num2str(frame_ID),'/',Lagrange_meshes_interp_order,'_mesh_FE.X'],[Sim_prep_folder,'/',Lagrange_meshes_interp_order,'_mesh_FE.X']);
                    % Cubic_mesh_FE.T
                    copyfile([Lag_mesh_folder_path,'/Frame_',num2str(frame_ID),'/',Lagrange_meshes_interp_order,'_mesh_FE.T'],[Sim_prep_folder,'/',Lagrange_meshes_interp_order,'_mesh_FE.T']);
                    % read nodes at deformed configuration:
                    fid=fopen([Lag_mesh_folder_path,'/Frame_',num2str(frame_ID),'/',Lagrange_meshes_interp_order,'_mesh_FE.X']);
                    Nodes_intro=fscanf(fid,'%d',[2,1]).';
                    Lagrange_mesh_Nodes_DEF=fscanf(fid,'%f',[3,Inf]).';
                    fclose(fid); 
                    % create linear mesh out of cubic mesh at frame
                    [Linear_mesh_Nodes,Linear_mesh_Elements,Linear_mesh_Boundaries]=create_Lin_mesh_from_Cubic_mesh_cheart_ordering(Lagrange_mesh_Nodes_DEF,Lagrange_mesh_elements_ala_CHeart_Ref,Lagrange_mesh_Boundaries,Nodes_per_elem_dir);
                else % copy reference frame mesh
                    % Cubic_mesh_FE.X
                    copyfile([Lag_mesh_folder_path,'/Frame_',num2str(reference_frame),'/',Lagrange_meshes_interp_order,'_mesh_FE.X'],[Sim_prep_folder,'/',Lagrange_meshes_interp_order,'_mesh_FE.X']);
                    % Cubic_mesh_FE.T
                    copyfile([Lag_mesh_folder_path,'/Frame_',num2str(reference_frame),'/',Lagrange_meshes_interp_order,'_mesh_FE.T'],[Sim_prep_folder,'/',Lagrange_meshes_interp_order,'_mesh_FE.T']);
                end
                % Cubic_mesh_FE.B
                B_file_name=[Lagrange_meshes_interp_order,'_mesh_FE.B'];
                create_Cheart_Bfile_from_Boundaries_mat(Lagrange_mesh_Boundaries,Nodes_per_elem_dir,Sim_prep_folder,B_file_name);
                % Cubic_mesh_FE.NBT            
                NBT_file_name=[Lagrange_meshes_interp_order,'_mesh_FE.NBT'];
                create_Cheart_NBTfile_from_Boundaries_mat(Lagrange_mesh_Neighbours,Sim_prep_folder,NBT_file_name);

                % Lin_mesh_FE.X
                X_file_name='Lin_mesh_FE.X';
                create_Cheart_Xfile_from_Nodes_mat(Linear_mesh_Nodes,Sim_prep_folder,X_file_name);

                % Lin_mesh_FE.T
                T_file_name='Lin_mesh_FE.T';
                create_Cheart_Tfile_from_Elements_mat(Linear_mesh_Elements,Linear_mesh_Nodes,2,Sim_prep_folder,T_file_name);

                % Lin_mesh_FE.B
                B_file_name='Lin_mesh_FE.B';
                create_Cheart_Bfile_from_Boundaries_mat(Linear_mesh_Boundaries,2,Sim_prep_folder,B_file_name);

                % Lin_mesh_FE.NBT just copy this, it's the same in Cubic and linear meshes      
                copyfile([Sim_prep_folder,'/',Lagrange_meshes_interp_order,'_mesh_FE.NBT'],[Sim_prep_folder,'/Lin_mesh_FE.NBT']);
                if flag_Daves_fibers==1
                    create_P_files_for_Daves_fiber_fitting_tool(Cheart_Fibers_code_folder,Sim_prep_folder,Lagrange_mesh_Nodes_Ref,Basal_Nodes_IDs,Endocardial_Nodes_IDs);

                elseif flag_my_fibers==1
                    % print fibers fitted to ref frame or def frame..:
                    if flag_VF_sim_only__fit_fibers_to_def_configuration==1 % only if this is =1 otherwise fibers need only be fit to reference frame (and you just print those one to the Def Frame where the prescribed end displacement is updated.
                        % if flag_VF_sim_only__fit_fibers_to_def_configuration==1 then I need to fit fibers to every deformed config
                        % read nodes at deformed configuration:
                        fid=fopen([Lag_mesh_folder_path,'/Frame_',num2str(frame_ID),'/',Lagrange_meshes_interp_order,'_mesh_FE.X']);
                        Nodes_intro=fscanf(fid,'%d',[2,1]).';
                        Lagrange_mesh_Nodes_DEF=fscanf(fid,'%f',[3,Inf]).';
                        fclose(fid); 
                        
                        if flag_disc_Fib==0
                            [Fibers_nodes,~,~,~]=fit_fibers_using_Laplace_Dirichlet_solver_no_imbric_angle(fiber_epi_degrees,fiber_endo_degrees,flag_disc_Fib,Lagrange_mesh_Elements_Cmgui,Lagrange_mesh_Nodes_DEF,Lagrange_mesh_Boundaries,Nodes_per_elem_dir,GP_per_elem_dir,Endo_patch_ID,Epi_patch_ID,Base_patch_ID,Apex_hole_patch_ID,[]);
                        end
                    end

                    if flag_disc_Fib==0
                        fid=fopen([Sim_prep_folder,'/FIBERS.X'],'w');
                        fprintf(fid,'%d %d\n',[size(Fibers_nodes,1),size(Fibers_nodes,2)].');
                        fprintf(fid,'%21.15f    %21.15f    %21.15f    %21.15f    %21.15f    %21.15f    %21.15f    %21.15f    %21.15f\n',Fibers_nodes.');
                        fclose(fid);                    
                    end

                elseif flag_fake_fibers==1
                    if flag_disc_Fib==0
                        fid=fopen([Sim_prep_folder,'/FIBERS.X'],'w');
                        fprintf(fid,'%d %d\n',[size(Fibers_nodes,1),size(Fibers_nodes,2)].');
                        fprintf(fid,'%21.15f    %21.15f    %21.15f    %21.15f    %21.15f    %21.15f    %21.15f    %21.15f    %21.15f\n',Fibers_nodes.');
                        fclose(fid);                    
                    end
                end

                if flag_create_VF_epi_NeoHookean_sin_mu3==1 % prep Prescr disps Epi:
                    
                    disp(['== estimating prescr disps Epi for frame: ',num2str(frame_ID),' ==']);
                    disp('');
                    
                    Prescr_Epi_Nodes=build_VF_field_epi_parallel_to_g_theta_prescr_deltaU_EPI_mu_cub(Lagrange_mesh_Nodes_DEF,Lagrange_mesh_Elements_Cmgui,Lagrange_mesh_Boundaries,Nodes_per_elem_dir,GP_per_elem_dir,Endo_patch_ID,Epi_patch_ID,Base_patch_ID,Apex_hole_patch_ID,fig_ind,max_disp);
                    
                    fid=fopen([Sim_prep_folder,'/Tags_Epi_Disp-',num2str(0),'.D'],'w');
                    fprintf(fid,'%d %d\n',[size(Prescr_Epi_Nodes,1),size(Prescr_Epi_Nodes,2)].');
                    fprintf(fid,'%21.15f    %21.15f    %21.15f\n    ',zeros(size(Prescr_Epi_Nodes)));
                    fclose(fid);
                    for n_sim=1:simulation_increments                                                 
                        % Prescr_Epi_Nodes: is zero everywhere except at epi. Prescr_Epi_Nodes is just the displacement vectors to be applied 
                        if isempty(fig_ind)==0   % if fig_ind is empty leave empty                         
                            fig_ind=fig_ind+1;
                        end
                        
                        fid=fopen([Sim_prep_folder,'/Tags_Epi_Disp-',num2str(n_sim),'.D'],'w');
                        fprintf(fid,'%d %d\n',[size(Prescr_Epi_Nodes,1),size(Prescr_Epi_Nodes,2)].');
                        fprintf(fid,'%21.15f    %21.15f    %21.15f\n    ',(Prescr_Epi_Nodes*n_sim/simulation_increments).');% you express the displacement as deformed - reference config
                        fclose(fid);
                    end
                    disp('=========================');
                end
                %% just also store basic data background info
                if flag_VF_sim_only__fit_fibers_to_def_configuration==0 && flag_prep_for_parameter_estimation==0 % || any other VF related flag ==0 and this is not run for parameter estimation prep either (in which case no simulations will be run in the folder)
                    fid=fopen([Sim_prep_folder,'/background_diastolic_data_used.txt'],'w');
                    fprintf(fid,'%s\n','% Frame Pressure Volume');
                    fprintf(fid,'%f   %f   %f\n', [diast_pressures_Pa_TU,Diastolic_Frames_Indices_TU,Diastolic_Frames_Volumes_ml_TU].');
                    fclose(fid);

                     %% infl.data % not needed any more
                    inflation_pressure=diast_pressures_Pa_TU(def_fr_ind);
                    dp_dt=-inflation_pressure/simulation_increments;
    %                 CreateInfldataFile_function(Sim_prep_folder,simulation_increments,-inflation_pressure);

                    %% .P file --this is removed from here to be more generic and avoid confusion
    % % % % %     %             make_cheart_Pfile_for_ventricle(sim_folder_dir,incr_num,C1,C2,C3,C4,patch_Endo,patch_Epi,patch_Base,patch_Apex,flag_base,flag_epi,flag_endo,flag_apex,flag_prescr_APEX_NODE,flag_disc_Fib,flag_fix_base)
    % % % % % %                 make_cheart_Pfile_for_ventricle(Sim_prep_folder,simulation_increments,C1,C2,C3,C4,Endo_patch_ID,Epi_patch_ID,Base_patch_ID,Apex_hole_patch_ID,1,0,0,0,0,flag_disc_Fib,0);

                    %% Prescribed displacements at base (Tags_Base_Disp-...D)
                    if flag_prescr_basal_disps==1
                        fid=fopen([Lag_mesh_folder_path,'/Frame_',num2str(frame_ID),'/',Lagrange_meshes_interp_order,'_mesh_FE.T']);            
                        Element_intro=fscanf(fid,'%d',[2,1]).';
                        Cubic_mesh_elements_ala_CHeart_Def=fscanf(fid,'%f',[64,Inf]).';
                        fclose(fid);

                        fid=fopen([Lag_mesh_folder_path,'/Frame_',num2str(frame_ID),'/',Lagrange_meshes_interp_order,'_mesh_FE.X']);
                        Nodes_intro=fscanf(fid,'%d',[2,1]).';
                        Cubic_mesh_nodes_Def=fscanf(fid,'%f',[3,Inf]).';
                        fclose(fid);

                        Disp=Cubic_mesh_nodes_Def-Lagrange_mesh_Nodes_Ref;


                        % print basal disps specifying files: Tags_Base_Disp-...D
                        Basal_Disps=zeros(size(Disp));
                        Basal_Disps(Basal_Nodes_IDs,:)=Disp(Basal_Nodes_IDs,:); % assign disp values only to basal nodes

                        % print file for Tags_Base_Disp-0.D
                        fid=fopen([Sim_prep_folder,'/Tags_Base_Disp-',num2str(0),'.D'],'w');
                        fprintf(fid,'%d %d\n',[size(Disp,1),size(Disp,2)].');
                        fprintf(fid,'%21.15f    %21.15f    %21.15f\n    ',zeros(size(Basal_Disps)));
                        fclose(fid);
                        % print file for Tags_Base_Disp-n.D for each n=1:simulation_increments
                        for n_sim=1:simulation_increments
                            fid=fopen([Sim_prep_folder,'/Tags_Base_Disp-',num2str(n_sim),'.D'],'w');
                            fprintf(fid,'%d %d\n',[size(Disp,1),size(Disp,2)].');
                            fprintf(fid,'%21.15f    %21.15f    %21.15f\n    ',Basal_Disps.'*(n_sim/simulation_increments));
                            fclose(fid);
                        end
                    end


                end
                if flag_prep_for_parameter_estimation==0 % no sims run if this is just for prepping parameter estimation, so no Pfile needed
                    %% now make Pfile                                                                                      
                    if flag_press_epi==1
                        if dp_dt>0
                            warning('pressure is positive - are you sure you want to apply deflation?')
                        end
                        dp_dt_epi=dp_dt/3;                        
                        press_expr_epi=[num2str(dp_dt_epi),'*t'];
                        press_expr_endo=[num2str(dp_dt),'*t'];
                        make_cheart_Pfile_for_LV_PI_endo_epi_expression_sel_constit(Sim_prep_folder,simulation_increments,press_expr_endo,press_expr_epi,press_expr_sept,material_law,...
                            material_parameters,Endo_patch_ID,Epi_patch_ID,Base_patch_ID,Apex_hole_patch_ID,Septum_patch_ID,flag_base_P_file,...
                            flag_epi_P_file,flag_endo_P_file,flag_apex_P_file,flag_prescr_APEX_NODE_P_file,flag_disc_Fib,flag_fix_base_P_file,...
                            flag_fix_apex_hole_P_file,Nodes_per_elem_dir);
                    else
                        make_cheart_Pfile_for_ventricle_PI_sel_constitutive(Sim_prep_folder,simulation_increments,dp_dt,material_law,material_parameters,...
                            Endo_patch_ID,Epi_patch_ID,Base_patch_ID,Apex_hole_patch_ID,flag_base_P_file,flag_epi_P_file,flag_endo_P_file,flag_apex_P_file,flag_prescr_APEX_NODE_P_file,...
                            flag_disc_Fib,flag_fix_base_P_file,flag_fix_apex_hole_P_file,Nodes_per_elem_dir);                                                                                       
                    end
                end
            else
            
                copyfile([Lag_mesh_folder_path,'/Frame_',num2str(frame_ID),'/',Lagrange_meshes_interp_order,'_mesh_FE.X'],[Sim_prep_folder,'/Def_mesh-',num2str(frame_ID),'.X']);
            end
    
        end
        
    end                    
    
    
else
    disp('modify script for use with other meshes except for new cases')
end

















