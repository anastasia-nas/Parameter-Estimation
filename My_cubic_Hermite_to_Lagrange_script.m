close all; clear all; clc;

% modified 10-October-2019 from My_cubic_Hermite_to_QUAD_Lagrange_script.m to accound for
% different Lagrange interpolations (cubic and quadratic) and also check if the Cubic hermite meshes have issues in their node interpolation order. 
% created on 14/8/2019 by modifying My_cubic_Hermite_to_cubic_Lagrange_script_home_1.m 
%----------------------------------------------------------------------------------------------
% purpose: cub lagrange meshes crash in cheart throwing memory error -
% trying if quad does the job


% created by Anastasia on 16-05-2019
%-------------------------------------------------------------------------------------------------------------------------
% Purpose: 
%-------------------------------------------------------------------------------------------------------------------------
% create a script to get Cubic Lagrange meshes out of the cubic Hermites
% that Pablo gave me (in cmgui .exnode/ .exelem files) from his mesh
% fitting pipeline
%-------------------------------------------------------------------------------------------------------------------------
% Notes on methodology: 
%-------------------------------------------------------------------------------------------------------------------------
% They way to build the CL meshes: I will choose equally spaced locations
% in the CL element natural coordinates space and for each of these points
% with known coordinates in xi1,xi2,xi3 space I will interpolate the
% x,y,z values using CH interpolation. I will assume the nodal values order is
% described correctly in the cmgui .exnode file:
%           x.  Value index= 1, #Derivatives= 7 (d/ds1,d/ds2,d2/ds1ds2,d/ds3,d2/ds1ds3,d2/ds2ds3,d3/ds1ds2ds3), #Versions= 1
%           y.  Value index= 9, #Derivatives= 7 (d/ds1,d/ds2,d2/ds1ds2,d/ds3,d2/ds1ds3,d2/ds2ds3,d3/ds1ds2ds3), #Versions= 1
%           z.  Value index= 17, #Derivatives= 7 (d/ds1,d/ds2,d2/ds1ds2,d/ds3,d2/ds1ds3,d2/ds2ds3,d3/ds1ds2ds3), #Versions= 1
% Remember: in cubic Hermite we are only using data values at the corner
% nodes of the element and the first derivatives, which  for the 3D element
% these are the 7 derivatives shown above (not sure why we take all of them
% other than that it fits the 64 constant requirement for a cubic
% polynomial in 3D (similarly to the 64-node Cub Lag. polynomial)
%-------------------------------------------------------------------------------------------------------------------------
% Input:
%-------------------------------------------------------------------------------------------------------------------------
% cmgui .exnode and .exelem files describing the mesh 
% or other?
% use a function that takes care of this (imports whatever files and exports a CH_Elements- CH_Nodes mesh description) and that can be built on later if
% for example I get another choice of input format


%-------------------------------------------------------------------------------------------------------------------------
% Output:
%-------------------------------------------------------------------------------------------------------------------------
% Cubic Lagrange meshes in Cheart format option to print.
%-------------------------------------------------------------------------------------------------------------------------
% SOS points:
%-------------------------------------------------------------------------------------------------------------------------
% Super-SOS point: Pablo is not doing this right --- He is skipping the
% scaling factors I think (the dx/dxi required to transform the df/dxi to
% df/dx (gradient) ). 
% The correct interpolation for each field coordinate y_i (here the field are the actual coordinates y=X, i=1,2,3) is:
% y_i= y^a_i H^(0,a)+dy_i/dx_j
% There might be an issue with the node numbering in Pablo's script.. Or in
% this one? Also correct it in the the Cubic Lagrange to Cubic Hermite
% version..

% Also the meshes are mm! Turn them into meters!- I chose to turn meshes
% dimensions into meters in order to be consistent and avoid having issues
% from meshes being in different scales across folders. Problem is of
% course I am loosing accuracy this way and increasing round-offs.

%=============================================================================================================================================================

% FLAGS:
res_str='341'; % choose from '341' or '9124'
i_case='OSU';
mesh_name=['MNaNMesh',res_str];
%mesh_name='MNaNMesh341';% Replace with 9124
Lagrange_interpolation_order='Cubic' ;  % choose between 'Cubic' and 'Quadratic'
folder_format_spec='use_new_cases_dir_format'; % this is for allowing some flexibility with folder paths.
CH_mesh_format='cmgui_ex_files';
CL_mesh_output='print';
plot_CL_mesh_at_frame=0; % plot generated CL in matlab for sanity check
flag_test_conversion_from_cmgui_to_cheart=0; %=1 if you want to  plot cmgui mesh after conversion from cmgui to cheart and back to cmgui ordering
flag_also_print_cmgui_ex_files_of_meshes=1; % =1 if you also want to print the meshes as .exnode and .exelem files to help with visualisation and debugging of next steps
% SPECS:
% Nodes_per_elem_dir=3; % QUADS!!!
fig_ind=0; % for numbering plots

% FOLDERS:
% gen_folder=['C:\Users\fitzp\Documents\Josh_Documents\KCL\iBSc\iBSc_4_Stiffness_Project\Manav_Cases_Meshed_',res_str];
gen_folder=['C:\Users\an11\Desktop\Work\Teaching\Stiffness_Project_Joshua_Kader\Manav_Cases_Meshed_',res_str];
% gen_folder='C:\Users\an11\Desktop\Work\new_cases_meshed'; % use from laptop
% gen_folder='/staff/an11/new_cases_meshed'; % use from work desktop
CubHerm_mesh_general_folder=[gen_folder,'/Cubic_Hermite_Meshes_From_Pablos_Toolkit'];

flag_rescale_mesh_to_meters=1; % = 1 if you want to turn mesh dimensions into meeters

% temp_folder_for_home='C:\Users\an11\Desktop\Work\Various_cheart_sims_and_visualisations\anisotropic_Fib_90deg_Daves_fibers_C1_100_alpha_15_Fixed_Base_and_Epi_DEFLATION\Simulation';
%=============================================================================================================================================================
switch Lagrange_interpolation_order
    case 'Cubic'
        Nodes_per_elem_dir=4; % CUBICS!!!
        Lagrange_mesh_general_folder=[gen_folder,'/Cubic_Lagrange_Meshes'];
        if flag_rescale_mesh_to_meters==1
            Lagrange_mesh_general_folder=[Lagrange_mesh_general_folder,'_in_m'];
        end
    case 'Quadratic'
        Nodes_per_elem_dir=3; % QUADS!!!
        Lagrange_mesh_general_folder=[gen_folder,'/Quadratic_Lagrange_Meshes'];
        if flag_rescale_mesh_to_meters==1
            Lagrange_mesh_general_folder=[Lagrange_mesh_general_folder,'_in_m'];
        end
end

CubHerm_mesh_folder_path=[CubHerm_mesh_general_folder,'/case_',i_case];
Lagrange_mesh_folder_path=[Lagrange_mesh_general_folder,'/case_',i_case];
mkdir(Lagrange_mesh_folder_path);

%% read the exnode files in folder to find number of frames:
mesh_exfile_format_str=[mesh_name,'_']; % assumes meshes are in formatl mplampla_10.exnode (unlike test-10.exnode)
% mesh_exfile_format_str=[mesh_name,'-'];
listings=dir([CubHerm_mesh_folder_path,'/',mesh_exfile_format_str,'*.exnode']);
% for each frame read meshes and save them - check Elements matrices are
% the same in each frame..
for n_frame=1:length(listings)
    
    prox=textscan(listings(n_frame).name,[mesh_exfile_format_str,'%d','.exnode']);
    exnode_string=listings(n_frame).name;
    start_ind=regexp(exnode_string,'.exnode','start');
    frame_mesh_name{n_frame}=exnode_string(1:(start_ind-1));
    frame_ID(n_frame)=prox{1};
end
% [frame_ID,frame_mesh_name_ind]=sort(frame_ID); % I don't really need these ordered as long as frame_ID and frame_mesh_name are correctly mapped (i.e. to frame_ID 0 frame_mesh_name MNanMesh9124_00)
% frame_mesh_name=frame_mesh_name(frame_mesh_name_ind);
for n_frame=1:length(frame_ID)
    CH_mesh_exnode_file=[CubHerm_mesh_folder_path,'/',frame_mesh_name{n_frame},'.exnode']; % this assumes the files are in format: Mesh9123_00.exnode 
    CH_mesh_exelem_file=[CubHerm_mesh_folder_path,'/',frame_mesh_name{n_frame},'.exelem']; % this assumes the files are in format: Mesh9123_00.exelem
    %%-----------CALL A FUNCTION TO READ THE CH MESH ELEMENTS AND NODES---------
    CH_Nodes_of_Fields_in_file=read_CH_Nodes_from_cmgui_exnode_file(CH_mesh_exnode_file);
    [CH_Elements,CH_Elements_scale_factors]=read_CH_Elements_from_cmgui_exelem_file(CH_mesh_exelem_file);
    %%------------FIND THE xi COORDINATES FOR THE INTERNAL LAGRANGE NODES-------
    % you only have nodes at xi=0 and xi=1 and you need the intermediate ones
    % (xi/3,2xi/3)
    for n1=1:Nodes_per_elem_dir
        for n2=1:Nodes_per_elem_dir
            for n3=1:Nodes_per_elem_dir
                ksi1=(n1-1)/(Nodes_per_elem_dir-1);
                ksi2=(n2-1)/(Nodes_per_elem_dir-1);
                ksi3=(n3-1)/(Nodes_per_elem_dir-1);
                Internal_node_locations_for_CL_mesh((n3-1)*Nodes_per_elem_dir^2+(n2-1)*Nodes_per_elem_dir+n1,:)=[ksi1,ksi2,ksi3];
            end
        end
    end

    %% -----------INTERPOLATE THE FIELD VALUES TO THE xi COORDINATES WITH CH INTEPROLATION-------------------
    field_to_interpolate='coordinate'; % this depends on specific cmgui mesh (could be 'Space' or 'coordinate' for X, 'Disp', 'Fib' whatever.
    fields_under_nodes=fieldnames(CH_Nodes_of_Fields_in_file); % read which fields were read from the exnode file
    for n_f=1:length(fields_under_nodes)
        curr_field_str=fields_under_nodes{n_f};
        coord_ID_vec(n_f)=textscan(curr_field_str,[field_to_interpolate,'_%d']); % vector with coordinate IDs under a certain 'field_to_interpolate'            
    end
    % do interpolation per coordinate
    
    for n_coord=1:length(coord_ID_vec)
        coun_new_node=0; % keep it inside the coordinate loop
        coord=coord_ID_vec{n_coord};
        CH_Nodes_field_coord=CH_Nodes_of_Fields_in_file.(fields_under_nodes{n_coord});
        Lagrange_Nodes(1:size(CH_Nodes_field_coord,1),coord)=CH_Nodes_field_coord(:,1); % apothikeyeis prwta ta nodes pou yparxoun sto CH mesh -  the first column is the nodal values you want the other 7 cols are derivatives.
        for m_el=1:size(CH_Elements,1)
            Elem_Nodes=CH_Nodes_field_coord(CH_Elements(m_el,:),:); % 8 nodes per element and 8 coordinates per node (1 value 7 derivatives following cmgui format)
            Elem_Weights_row=CH_Elements_scale_factors(m_el,:);%ordered according to cmgui ordering 1 value and 7derivatives per node
            Elem_nodes_row=[Elem_Nodes(1,:),Elem_Nodes(2,:),Elem_Nodes(3,:),Elem_Nodes(4,:),Elem_Nodes(5,:),Elem_Nodes(6,:),Elem_Nodes(7,:),Elem_Nodes(8,:)];
            Elem_weighted_nodes_row=Elem_nodes_row.*Elem_Weights_row; % size: 8 rows (one row per node ) 8 columns (one column per value, node or derivative
            

            for int_nod=1:size(Internal_node_locations_for_CL_mesh,1)
                ksi_location=Internal_node_locations_for_CL_mesh(int_nod,:);
                %make sure you don't reassign corner nodes -no point in that
                if sum(abs(ksi_location-[0,0,0]))~=0 && sum(abs(ksi_location-[1,0,0]))~=0 && sum(abs(ksi_location-[0,1,0]))~=0 && sum(abs(ksi_location-[1,1,0]))~=0 && sum(abs(ksi_location-[0,0,1]))~=0 && sum(abs(ksi_location-[1,0,1]))~=0 && sum(abs(ksi_location-[0,1,1]))~=0 && sum(abs(ksi_location)-abs([1,1,1]))~=0
                    coun_new_node=coun_new_node+1;
                    new_node_index_in_mesh=size(CH_Nodes_field_coord,1)+coun_new_node;
                    Shape_function_vector_cmgui_ordering=Hermite_basis_functions_in_3D_VECTOR_cmgui_ordering(ksi_location(1),ksi_location(2),ksi_location(3)); % row vector size (1,64)
                    Lagrange_Nodes(new_node_index_in_mesh,coord)=sum(Elem_weighted_nodes_row.*Shape_function_vector_cmgui_ordering);
                    Lagrange_Elements(m_el,int_nod)=new_node_index_in_mesh;
                else
                    if sum(abs(ksi_location-[0,0,0]))==0
                        Lagrange_Elements(m_el,int_nod)=CH_Elements(m_el,1);
                    elseif sum(abs(ksi_location-[1,0,0]))==0
                        Lagrange_Elements(m_el,int_nod)=CH_Elements(m_el,2);
                    elseif sum(abs(ksi_location-[0,1,0]))==0
                        Lagrange_Elements(m_el,int_nod)=CH_Elements(m_el,3);
                    elseif sum(abs(ksi_location-[1,1,0]))==0
                        Lagrange_Elements(m_el,int_nod)=CH_Elements(m_el,4);
                    elseif sum(abs(ksi_location-[0,0,1]))==0
                        Lagrange_Elements(m_el,int_nod)=CH_Elements(m_el,5);
                    elseif sum(abs(ksi_location-[1,0,1]))==0
                        Lagrange_Elements(m_el,int_nod)=CH_Elements(m_el,6);
                    elseif sum(abs(ksi_location-[0,1,1]))==0
                        Lagrange_Elements(m_el,int_nod)=CH_Elements(m_el,7);
                    elseif sum(abs(ksi_location-[1,1,1]))==0
                        Lagrange_Elements(m_el,int_nod)=CH_Elements(m_el,8);
                    end
                end
                        
                
            end
        end
    end
    
   
    %% remove multiple nodes and make node ids in elements unique      
    if n_frame==1
        % use the same unique sorting for all frames- so just do it for
        % one.
%         [CL_Nodes_UNIQUE_matlab1,i_mult1,i_uniq1]=unique(CL_Nodes,'rows'); %to idio mou vgazei me to uniquetol, alla gia asfaleia xrisimopoiise to uniquetol dioti mporei sto test case (BAL) na min yparxei round off error alla isws se alla cases na yparxei.
%         %remember that [C,i_a,i_c]=unique(A,'rows') means: C(i_c,:)=A(i_a,:) 
        [Lagrange_Nodes_UNIQUE_matlab,i_mult,i_uniq]=uniquetol(Lagrange_Nodes,10^(-5),'ByRows',true);
        for n=1:length(i_uniq)
            for m_el=1:size(Lagrange_Elements,1)
                Lagrange_Elements_UNIQUE_matlab(m_el,:)=i_uniq(Lagrange_Elements(m_el,:));
            end
        end
    else
        Lagrange_Nodes_UNIQUE_matlab=Lagrange_Nodes(i_mult,:);
    end
    
    disp(['unique sorting finished for frame ',num2str(frame_ID(n_frame))]);  
    if plot_CL_mesh_at_frame==1
        fig_ind=fig_ind+1;
        figure(fig_ind)
        plot_cmgui_ordered_mesh(Lagrange_Elements_UNIQUE_matlab,Lagrange_Nodes_UNIQUE_matlab,Nodes_per_elem_dir,3); % this works fine, so cmgui ordering is fine
    end
    
    
    %% Check ksi1-ksi2-ksi3 form a right hand system (this was a problem with some of the meshes - all except for WIL as far as I could tell). -- this is important for doing calculations in the elements and not getting negative J's (thetaX/thetaksi)
    % not sure where the original problem comes from - Pablo checked the
    % template mesh he was using and the xi direction was correct there (I
    % also saw that)
    Lagrange_Elements_UNIQUE_matlab_pre=Lagrange_Elements_UNIQUE_matlab;
    Lagrange_Elements_UNIQUE_matlab_fin=Lagrange_Elements_UNIQUE_matlab;
    for m_el=1:size(Lagrange_Elements_UNIQUE_matlab,1) % note elements are in cmgui ordering (no corner node priority)
        ksi1_dir=(Lagrange_Nodes_UNIQUE_matlab(Lagrange_Elements_UNIQUE_matlab(m_el,Nodes_per_elem_dir),:)-Lagrange_Nodes_UNIQUE_matlab(Lagrange_Elements_UNIQUE_matlab(m_el,1),:))/norm(Lagrange_Nodes_UNIQUE_matlab(Lagrange_Elements_UNIQUE_matlab(m_el,Nodes_per_elem_dir),:)-Lagrange_Nodes_UNIQUE_matlab(Lagrange_Elements_UNIQUE_matlab(m_el,1),:));
        ksi2_dir=(Lagrange_Nodes_UNIQUE_matlab(Lagrange_Elements_UNIQUE_matlab(m_el,Nodes_per_elem_dir*(Nodes_per_elem_dir-1)+1),:)-Lagrange_Nodes_UNIQUE_matlab(Lagrange_Elements_UNIQUE_matlab(m_el,1),:))/norm(Lagrange_Nodes_UNIQUE_matlab(Lagrange_Elements_UNIQUE_matlab(m_el,Nodes_per_elem_dir*(Nodes_per_elem_dir-1)+1),:)-Lagrange_Nodes_UNIQUE_matlab(Lagrange_Elements_UNIQUE_matlab(m_el,1),:));
        ksi3_dir=(Lagrange_Nodes_UNIQUE_matlab(Lagrange_Elements_UNIQUE_matlab(m_el,Nodes_per_elem_dir^2*(Nodes_per_elem_dir-1)+1),:)-Lagrange_Nodes_UNIQUE_matlab(Lagrange_Elements_UNIQUE_matlab(m_el,1),:))/norm(Lagrange_Nodes_UNIQUE_matlab(Lagrange_Elements_UNIQUE_matlab(m_el,Nodes_per_elem_dir^2*(Nodes_per_elem_dir-1)+1),:)-Lagrange_Nodes_UNIQUE_matlab(Lagrange_Elements_UNIQUE_matlab(m_el,1),:));
        
        if dot(cross(ksi1_dir,ksi2_dir),ksi3_dir)<0 % if cross product of ksi1, ksi2 has opposite direction to ksi3 - then change ksi1 order (I keep ksi2, ksi3 because these are the way I like them by convention)
            % reorder along ksi1: 
            % I will change ksi1 because I want to preserve ksi2, ksi3 order as these work well in Pablo's meshes. (also because just changing one direction in a 3 vector system that is not right handed is enough to make it right handed)
            New_node_order=find_node_order_for_change_in_ksi_direction_in_Lagrange_meshes(Nodes_per_elem_dir,1);
            Lagrange_Elements_UNIQUE_matlab_fin(m_el,:)=Lagrange_Elements_UNIQUE_matlab(m_el,New_node_order);
            disp(['reordering element ': num2str(m_el)]);
        end
        
    end
    Lagrange_Elements_UNIQUE_matlab=Lagrange_Elements_UNIQUE_matlab_fin; % assign original name so there is no change in script if no change in ordering takes place
    
    %% check the scaling is right and save in structure: (mesh BAL is in mm - not sure this is the case for the remaining meshes-hence do checks and turn everything to meters to be compatible with Cheart simulations
    if flag_rescale_mesh_to_meters==1
        dim_X=abs(max(Lagrange_Nodes_UNIQUE_matlab(:,1))-min(Lagrange_Nodes_UNIQUE_matlab(:,1)));
        dim_Y=abs(max(Lagrange_Nodes_UNIQUE_matlab(:,2))-min(Lagrange_Nodes_UNIQUE_matlab(:,2)));
        dim_Z=abs(max(Lagrange_Nodes_UNIQUE_matlab(:,3))-min(Lagrange_Nodes_UNIQUE_matlab(:,3)));

        % taking into advantage that adult human heart dimensions will surely
        % be between 2 and 20 cm:
        if max([dim_X,dim_Y,dim_Z])<=0.2 && max([dim_X,dim_Y,dim_Z])>=0.02 % mesh is in m - all good!


        elseif max([dim_X,dim_Y,dim_Z])<=2 && max([dim_X,dim_Y,dim_Z])>=0.2 % mesh is in dm - divide by 10
            scale_TU=0.1;
            Lagrange_Nodes_UNIQUE_matlab=scale_TU*Lagrange_Nodes_UNIQUE_matlab;

        elseif max([dim_X,dim_Y,dim_Z])<=20 && max([dim_X,dim_Y,dim_Z])>=2 % mesh is in cm - divide by 100
            scale_TU=0.01;
            Lagrange_Nodes_UNIQUE_matlab=scale_TU*Lagrange_Nodes_UNIQUE_matlab;

        elseif max([dim_X,dim_Y,dim_Z])<=200 && max([dim_X,dim_Y,dim_Z])>=20 % mesh is in mm - divide by 1000
            scale_TU=0.001;
            Lagrange_Nodes_UNIQUE_matlab=scale_TU*Lagrange_Nodes_UNIQUE_matlab;
        end
    end
        
    Lagrange_Nodes_per_frame.(['Frame_',num2str(frame_ID(n_frame))])=Lagrange_Nodes_UNIQUE_matlab;
        
    %% -----------CALL FUNCTION TO PRINT TO CL FORMAT IN NEW FOLDER--------------
    %use flag folder_format_spec input to specify the format of directory where
    %you will print this.
%     1) turn cmgui ordering into cheart ordering

%     2) print each cheart mesh
    mesh_folder=[Lagrange_mesh_folder_path,'/Frame_',num2str(frame_ID(n_frame))];
    mkdir(mesh_folder);
    [Lagrange_Elements_UNIQUE_matlab_ala_cheart,error]=turn_Cmgui_node_ordering_to_Cheart(Lagrange_Elements_UNIQUE_matlab,Nodes_per_elem_dir);
    % test your conversions back and forth work:
    if flag_test_conversion_from_cmgui_to_cheart==1
        [Lagrange_Elements_UNIQUE_matlab_ala_CMGUI_From_cheart,error]=turn_Cheart_node_ordering_to_Cmgui(Lagrange_Elements_UNIQUE_matlab_ala_cheart,Nodes_per_elem_dir);
        if plot_CL_mesh_at_frame==1
            fig_ind=fig_ind+1;
            figure(fig_ind)
            plot_cmgui_ordered_mesh(Lagrange_Elements_UNIQUE_matlab_ala_CMGUI_From_cheart,Lagrange_Nodes_UNIQUE_matlab,Nodes_per_elem_dir,3);
        end
    end
    elements_str=[];
    for n_el=1:Nodes_per_elem_dir^3
        elements_str=[elements_str,'    %d'];
    end
    elements_str=[elements_str,'\n'];
    fid=fopen([mesh_folder,'/',Lagrange_interpolation_order,'_mesh_FE.T'],'w');
    fprintf(fid,'    %d    %d\n',[size(Lagrange_Elements_UNIQUE_matlab_ala_cheart,1),size(Lagrange_Nodes_UNIQUE_matlab,1)].');
%     fprintf(fid,'    %d    %d    %d    %d    %d    %d    %d    %d    %d    %d    %d    %d    %d    %d    %d    %d    %d    %d    %d    %d    %d    %d    %d    %d    %d    %d    %d    %d    %d    %d    %d    %d    %d    %d    %d    %d    %d    %d    %d    %d    %d    %d    %d    %d    %d    %d    %d    %d    %d    %d    %d    %d    %d    %d    %d    %d    %d    %d    %d    %d    %d    %d    %d    %d\n',CL_Elements_UNIQUE_matlab_ala_cheart.');
    fprintf(fid,elements_str,Lagrange_Elements_UNIQUE_matlab_ala_cheart.');
    fclose(fid);        
    
    fid=fopen([mesh_folder,'/',Lagrange_interpolation_order,'_mesh_FE.X'],'w');
    fprintf(fid,'    %d    %d\n',[size(Lagrange_Nodes_UNIQUE_matlab,1),3].');
    fprintf(fid,'    %21.16f    %21.16f    %21.16f\n',Lagrange_Nodes_UNIQUE_matlab.');
    fclose(fid);
    
    
    if flag_also_print_cmgui_ex_files_of_meshes==1
        view_all_mesh_folder=[Lagrange_mesh_folder_path,'/Visualise_all_frames'];
        mkdir(view_all_mesh_folder);
        GroupName=['vis_frame_',num2str(frame_ID(n_frame))];
        exelem_filename=[GroupName,'.exelem'];
        exnode_filename=[GroupName,'.exnode'];
        create_exnode_exelem_files_for_Space_var_and_Lagrange_interp(GroupName,exelem_filename,exnode_filename,view_all_mesh_folder,Lagrange_Nodes_UNIQUE_matlab,Lagrange_Elements_UNIQUE_matlab,Nodes_per_elem_dir)
        
        if n_frame==1 && strcmp(i_case,'BAL')==0 % no need to copy this for each frame.. 
%            copyfile([gen_folder,'/new_cases_hole/Cubic_Lagrange_Meshes','/case_BAL/Visualise_all_frames/view_beat.com'],[Lagrange_mesh_folder_path,'/Visualise_all_frames/view_beat.com']);
        end
    end
    
end

mal=date;
% save(['workspace_My_cubic_Hermite_to_quad_Lagrange_script_case_',i_case,'_',mal]);

