function [Wint_wrt_alpha_per_frame, virtual_Wint_wrt_alpha_per_frame]=calc_Wint_and_virtWint_per_frame_and_a_for_CF_Guccione(alpha_PS_vector,r_ff,r_fs,r_sn,diastolic_frames_list,Ref_Nodes,Elements_cmgui,Fibers_discon_nodes,Fibers_discon_elements_cmgui,Def_Nodes_per_Frame_struct,virt_disp_Nodes_per_Frame_struct,GP_per_elem_dir,Nodes_per_elem_dir)
% % % %  %% for debug then remove:        
% % % % diastolic_frames_list=Diastolic_Frames_Indices_TU;
% % % % Ref_Nodes=Ref_mesh_Nodes;
% % % % Elements_cmgui=Ref_mesh_Elements_ala_Cmgui;
% % % % Fibers_discon_nodes=Fib_mesh_Nodes;
% % % % Fibers_discon_elements_cmgui=Ref_mesh_Elements_ala_Cmgui;

% created by Anastasia on 23-02-2020
%-------------------------------------------------------------------------------------------------------------------------
% Purpose: 
%-------------------------------------------------------------------------------------------------------------------------
% put the estimation of Wint and delta_Wint into a function so this can be
% recycled after debugging
%-------------------------------------------------------------------------------------------------------------------------
% Notes on methodology:
%-------------------------------------------------------------------------------------------------------------------------
% maximum vectorisation, preallocation and economy
%-------------------------------------------------------------------------------------------------------------------------
% Input:
%-------------------------------------------------------------------------------------------------------------------------
% alpha_PS_vector: vector of alpha values in taken in parameter sweep
% diastolic_frames_list: row vector with the IDs of the diastolic frames (e.g. [25,26,17,1,2,3])
% Ref_Nodes: mesh nodes at reference frame
% Elements_cmgui:elements with with cmgui ordering (no priority at corner nodes)
% Fibers_discon_nodes: discontinuous fibers mesh nodes
% Fibers_discon_elements: element connectivity of the discontinuous fiber mesh (if mesh is continuous then use "Elements_cmgui"
% Def_Nodes_per_Frame_struct: a structure containing the deformed mesh at each frame (e.g. Def_Nodes_per_Frame_struct.Frame_25,Def_Nodes_per_Frame_struct.Frame_27, Def_Nodes_per_Frame_struct.Frame_1 etc..
% virt_disp_Nodes_per_Frame_struct: virtual field used for each deformed state, again in a structure: virt_disp_Nodes_per_Frame_struct.Frame_25, virt_disp_Nodes_per_Frame_struct.Frame_27, virt_disp_Nodes_per_Frame_struct.Frame_1 etc
%                                   if =[] then no virtual works calculations will take place and all respective VFB variables will be left empty
% GP_per_elem_dir: number of GPs per element dimension
% Nodes_per_elem_dir: Lagrange mesh nodes number per element dimension indicating the interpolation order (=4 for Cubic).


%-------------------------------------------------------------------------------------------------------------------------
% Output:
%-------------------------------------------------------------------------------------------------------------------------
% Wint_wrt_alpha: matrix of Wint per alpha in alpha_PS_vector per frame (size: [length(alpha_PS_vector),length(diastolic_frames_list)])
% virtual_Wint_wrt_alpha: matrix of \deltaWint per alpha in alpha_PS_vector per frame (size: [length(alpha_PS_vector),length(diastolic_frames_list)])
%-------------------------------------------------------------------------------------------------------------------------
% SOS points:
%-------------------------------------------------------------------------------------------------------------------------
% Throughout this process C1 is taken =1.

%% ==============================================START==============================================================
%% anisotropy ratios matrix:
r_coeff_row=[r_ff,r_fs,r_fs,r_fs,r_sn,r_sn,r_fs,r_sn,r_sn];
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

 
%% first do work on quantities calculated on the reference frame: thetaX/thetaksi, inv(thetaX_thetaksi) dV=det(thetaX/thetaksi), Rotation tensor
%preallocations:
Rotation_matrix_per_GP_mat=zeros(size(Elements_cmgui,1)*size(GPcoords,1),9);
Total_weights_per_GP_mat=zeros(size(Elements_cmgui,1)*size(GPcoords,1),1);
dV_per_GP_mat=zeros(size(Elements_cmgui,1)*size(GPcoords,1),1);
inv_thetaX_thetaksi_per_GP_mat=zeros(size(Elements_cmgui,1)*size(GPcoords,1),9);
thetaX_thetaksi_per_GP_mat=zeros(size(Elements_cmgui,1)*size(GPcoords,1),9);
for m_el=1:size(Elements_cmgui,1)
    X_at_Element_Nodes=Ref_Nodes(Elements_cmgui(m_el,:),:);
    Fibers_at_Element_Nodes=Fibers_discon_nodes(Fibers_discon_elements_cmgui(m_el,:),:);
    for n_GP=1:size(GPcoords,1)
        GP_ID_total=(m_el-1)*size(GPcoords,1)+n_GP;
        thetaN_thetaksi=squeeze(thetaN_thetaksi_per_GP_mat(:,:,n_GP)); % [Nodes_per_elem_dir^3,3]
        N_atksi=N_at_ksi_per_GP_mat(:,n_GP); % col [Nodes_per_elem_dir^3,1]
        Fiber_vectors_all=N_atksi.'*Fibers_at_Element_Nodes; % [1,9] row vector
        e_f=Fiber_vectors_all(1:3); e_s=Fiber_vectors_all(4:6); e_n=Fiber_vectors_all(7:9);
        Rotation_matrix_per_GP_mat(GP_ID_total,:)=[dot(e_f,e_x),dot(e_f,e_y),dot(e_f,e_z),dot(e_s,e_x),dot(e_s,e_y),dot(e_s,e_z),dot(e_n,e_x),dot(e_n,e_y),dot(e_n,e_z)]; % when this is turned back into a 3x3 matrix it will serve for turning vector from cartesian to fiber coors: v_f=R*v_c
        thetaX_thetaksi=X_at_Element_Nodes.'*thetaN_thetaksi;
        inv_thetaX_thetaksi=inv(thetaX_thetaksi);
        thetaX_thetaksi_per_GP_mat(GP_ID_total,:)=[thetaX_thetaksi(1,:),thetaX_thetaksi(2,:),thetaX_thetaksi(3,:)]; %thetaX/thetaksi1,thetaXthetaksi2,thetaXthetaksi3,thetaY/thetaksi1,thetaYthetaksi2,thetaYthetaksi3,thetaZ/thetaksi1,thetaZthetaksi2,thetaZthetaksi3]
        inv_thetaX_thetaksi_per_GP_mat(GP_ID_total,:)=[inv_thetaX_thetaksi(1,:),inv_thetaX_thetaksi(2,:),inv_thetaX_thetaksi(3,:)]; 
        dV_per_GP_mat(GP_ID_total,:)=det(thetaX_thetaksi);
        Total_weights_per_GP_mat(GP_ID_total,:)=GPweights_mat(n_GP,1)*GPweights_mat(n_GP,2)*GPweights_mat(n_GP,3);
    end
    
end

%% now read each frame:
% preallocate:
Wint_wrt_alpha_per_frame=zeros(length(alpha_PS_vector),length(diastolic_frames_list));
virtual_Wint_wrt_alpha_per_frame=zeros(length(alpha_PS_vector),length(diastolic_frames_list));
Q_over_a_per_GP_mat=zeros(size(Elements_cmgui,1)*size(GPcoords,1),1);
r_Efib_DEdeltau_per_GP_mat=zeros(size(Elements_cmgui,1)*size(GPcoords,1),1);
for n_DF=1:length(diastolic_frames_list)
    Def_Nodes=Def_Nodes_per_Frame_struct.(['Frame_',num2str(diastolic_frames_list(n_DF))]);
    if isempty(virt_disp_Nodes_per_Frame_struct)==0
        Virtual_Disp_Nodes=virt_disp_Nodes_per_Frame_struct.(['Frame_',num2str(diastolic_frames_list(n_DF))]);
    else
        Virtual_Disp_Nodes=[];
    end
    
    % calculate Wint per alpha per frame:
    for m_el=1:size(Elements_cmgui,1)
        x_at_Element_Nodes=Def_Nodes(Elements_cmgui(m_el,:),:);
        if isempty(Virtual_Disp_Nodes)==0
            virtu_at_Element_Nodes=Virtual_Disp_Nodes(Elements_cmgui(m_el,:),:);
        end
        for n_GP=1:size(GPcoords,1)
            GP_ID_total=(m_el-1)*size(GPcoords,1)+n_GP;
            thetaN_thetaksi=squeeze(thetaN_thetaksi_per_GP_mat(:,:,n_GP)); % [Nodes_per_elem_dir^3,3]
            thetax_thetaksi=x_at_Element_Nodes.'*thetaN_thetaksi;
            thetaksi_thetaX=[inv_thetaX_thetaksi_per_GP_mat(GP_ID_total,1:3);inv_thetaX_thetaksi_per_GP_mat(GP_ID_total,4:6);inv_thetaX_thetaksi_per_GP_mat(GP_ID_total,7:9)];
            F=thetax_thetaksi*thetaksi_thetaX;
            E=0.5*(F.'*F-eye(3));
            % turn everything to fiber coords:
            R=[Rotation_matrix_per_GP_mat(GP_ID_total,1:3);Rotation_matrix_per_GP_mat(GP_ID_total,4:6);Rotation_matrix_per_GP_mat(GP_ID_total,7:9)];
            E_fib=R*E*R.';
            E_fib_row=[E_fib(1,:),E_fib(2,:),E_fib(3,:)];
            Q_over_a_per_GP_mat(GP_ID_total,:)=sum(r_coeff_row.*(E_fib_row.^2));
            if isempty(Virtual_Disp_Nodes)==0
                thetavirtu_thetaksi=virtu_at_Element_Nodes.'*thetaN_thetaksi;
                thetavirtu_thetaX=thetavirtu_thetaksi*thetaksi_thetaX;
                DE_virtu=0.5*(thetavirtu_thetaX.'*F+F.'*thetavirtu_thetaX);
                % turn everything to fiber coords:
                DE_virtu_fib=R*DE_virtu*R.';
                DE_virtu_fib_row=[DE_virtu_fib(1,:),DE_virtu_fib(2,:),DE_virtu_fib(3,:)];
                r_Efib_DEdeltau_per_GP_mat(GP_ID_total,:)=sum(r_coeff_row.*E_fib_row.*DE_virtu_fib_row); % (r\hadamard E_f) : DE_f[du] - multiply this by alpha*e^Q dV w_GP and sum over GPs for virt Wint
            end                                    
        end
    end
    
    % now do the parameter sweep.
    [Q_over_a_mat,Alpha_sweep_mat]=meshgrid(Q_over_a_per_GP_mat,alpha_PS_vector);
    Q_sweep_mat=Alpha_sweep_mat.*Q_over_a_mat; % col: GP, row: alpha
    e_Q_sweep_mat=exp(Q_sweep_mat);% col: GP, row: alpha
    Wint_temp=0.5*(e_Q_sweep_mat-ones(size(e_Q_sweep_mat)))*(dV_per_GP_mat.*Total_weights_per_GP_mat); % col vector      
    Wint_wrt_alpha_per_frame(:,n_DF)=Wint_temp;
    if isempty(Virtual_Disp_Nodes)==0
        virt_Wint_temp=(e_Q_sweep_mat.*Alpha_sweep_mat)*(r_Efib_DEdeltau_per_GP_mat.*dV_per_GP_mat.*Total_weights_per_GP_mat); % col vector
        virtual_Wint_wrt_alpha_per_frame(:,n_DF)=virt_Wint_temp;
    else
        virtual_Wint_wrt_alpha_per_frame=[];
    end
end