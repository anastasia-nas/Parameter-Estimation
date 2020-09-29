function make_cheart_Pfile_for_LV_PI_endo_epi_expression_sel_constit(sim_folder_dir,incr_num,press_expr_endo,press_expr_epi,press_expr_sept,material_law,material_parameters,patchID_Endo,patchID_Epi,patchID_Base,patchID_Apex,patchID_Sept,flag_prescr_bound_base,flag_prescr_bound_epi,flag_prescr_bound_endo,flag_prescr_bound_apex_hole,flag_prescr_APEX_NODE,flag_disc_Fib,flag_fix_base,flag_fix_apex_hole,Nodes_per_elem_dir)

% modified 24/02/2020 from make_cheart_Pfile_for_ventricle_new.m to account
% for different material laws as well as parameters. Also it allows you to
% prescribe the passive inflation with an expression without the need for a
% infl.data file.

% newer version of script created 15 Jan 2020 to account for Quadratic
% interpolation for disp/Space/fibers.

% function make_cheart_Pfile_for_ventricle(sim_folder_dir,incr_num,Disp_interpl,Pres_intrpl,Fib_intrpl,mat_law,mat_params,patch_Endo,patch_Epi,patch_Base,patch_Apex,flag_base,flag_epi,flag_endo,flag_apex,flag_disc_Fib)

% press_expr_endo: (=[] for no pressure at endo) - otherwise an expression with dp/dt*t where dp/dt must be <0 for inflation 
%                  (e.g. '-10*t' or '-10*t*Space.3/Z_scale' where Z_scale the maximum dimension along Z (to scale things)
% press_expr_epi: (=[] for no pressure at epi) - same as above but usually no gradient here..
% press_expr_sept: (=[]) if no septal ID is defined or no pressure is
% assigned to it...
% dp_dt: (=[] for no inflation) the time derivative of pressure (the pressure increment per time increment e.g. for dp_dt=-10 the infl.data file would be equivalent to 0,-10,-20,-30,....
% material_law: choose from 'Guccione', 'NeoHookean', ... add more Laws..
% material_parameters: [C1,C2,C3,C4] (row vector) in case of Guccione, C (scalar) in case of neohookean
% sim_folder_dir: full path to simulation folder
% incr_num: total number of simulation increments
% Disp_interpl: the interpolation used for the displacement field : 'Cub' or
%   'Quad' -> this script only runs for Cub currently..
% Pres_intrpl: interpolation of the pressure field: can be 'Lin' or 'Quad',
%    works only for Lin currently
% Fib_intrpl: Interpolation of fiber field is assumed to be the same as
%   displacement field, unless fiber field is discontinuous..
% mat_law: material  law currently only 'Guccione'
% mat_params=[C1,C2,C3,C4]; for Guccione material law
% how many boundary conditions i apply
% patch_Endo: the number ID of endocardial patch, same for patch_Epi,
% patch_Base (usually it's :patch_Endo=1, patch_Epi=2, patch_Base=3 for
% Roberts meshes)
% patchID_Sept: leave empty (=[]) if no septum is defined, otherwise 
%           assign a boundary patch if you have a septal part in the epi (I chose a group of
%           elements that are loaded with increased epicardial pressure in some
%           simulations to simulate the RV pressure on the septum)

% flag_base: =1 if you prescribe displacements at the entire base and thus
% need to define an extra variable the Base_Disp whose values at each
% increment will be updated by reading the file Tags_Base_Disp-incr_num.D
% or =0 if you don't prescribe disps that way at the base

% flag_epi, flag_endo: =1/0 same as flag_base

% flag_disc_Fib: =1 if the fiber field is discontinuous

% Nodes_per_elem_dir: Lagrange Interpolation (2:linear, 3: quad, 4:cubic)
% gia na fix base : gia in silico sims pou thela na dokimasw : an
% flag_fix_base=1 fixed base alliws vaze to 0..
%ftiakse kai to Ellipse
    curr_frame_new_sim_folder=sim_folder_dir;
    total_incr_curr_frame=incr_num;
%     if flag_apex==1
%         %         ftiakse to kainourgio .B file:
%         Cubic_Boundaries_New=[Patch_1; Patch_2; Patch_3; Patch_4];
%         
%         fid=fopen([curr_frame_new_sim_folder,'/Cubic_mesh_FE.B'],'w');
%         fprintf(fid,'%d\n',size(Cubic_Boundaries_New,1));
%         for mboun=1:size(Cubic_Boundaries_New,1)
%             for nboun=1:size(Cubic_Boundaries_New,2)
%                 fprintf(fid,'%d       ',Cubic_Boundaries_New(mboun,nboun));
%             end
%             fprintf(fid,'\n');
%         end
%         fclose(fid);
%     end
    %ftiakse kai to Ellipse
% %     if isempty(dp_dt)==0
% %         if dp_dt>0
% %             warning('!pressure time derivative is positive - this is equivalent to !DEFLATION in "make_cheart_Pfile_for_ventricle_PI_sel_constitutive" ')
% %         end
% %     end
    if strcmp(material_law,'Guccione')==1
        if length(material_parameters)~=4
            warning('error in "make_cheart_Pfile_for_ventricle_PI_sel_constitutive" size of material_paramemeters vector should be 4 and it is not!!!');
            disp('stopping script.....')
            disp('');
            return
            
        end
    elseif strcmp(material_law,'NeoHookean')==1
        if length(material_parameters)~=1
            warning('error in "make_cheart_Pfile_for_ventricle_PI_sel_constitutive" size of material_paramemeters vector should be 1 and it is not!!!');
            disp('stopping script.....')
            disp('');
            return
            
        end
    else % either a typo in material_law specification or material law not created here- correct this:
        warning('unspecified material_law in "make_cheart_Pfile_for_ventricle_PI_sel_constitutive" - correct this');
        disp('stopping script.....')
        return
    end
    fid=fopen([curr_frame_new_sim_folder,'/Ellipse.P'],'w');
    fprintf(fid,'%s\n','%  ---  Solver Groups & Matrices  --------------------------------------------');
    fprintf(fid,'%s\n','!DefSolverGroup={SolidTest|timestepping|SolidMatrix}');
    fprintf(fid,'%s\n','!DefSolverSubGroup={SolidTest|SOLVER_SEQUENTIAL_FIXED_POINT|USE_LINESEARCH|SolidMatrix}');
    fprintf(fid,'%s\n','!DefSolverMatrix={SolidMatrix|SOLVER_SuperLU|solid}');
    fprintf(fid,'%s\n','!SetSolverGroup={SolidTest|L2TOL|5e-9}');
    if flag_prescr_bound_base==1
        fprintf(fid,'%s\n','!SetSolverGroup={SolidTest|AddVariables|BaseDisp}');
    end
    if flag_prescr_bound_epi==1
        fprintf(fid,'%s\n','!SetSolverGroup={SolidTest|AddVariables|EpiDisp}');
    end
    if flag_prescr_bound_endo==1
        fprintf(fid,'%s\n','!SetSolverGroup={SolidTest|AddVariables|EndoDisp}');
    end
    if flag_prescr_bound_apex_hole==1
        fprintf(fid,'%s\n','!SetSolverGroup={SolidTest|AddVariables|ApexDisp}');
    end
    
    fprintf(fid,'%s\n','!DefTimeStepScheme={TimeStepping}');
    fprintf(fid,'   %d   %d   %d\n\n',[1,total_incr_curr_frame,1].');
    fprintf(fid,'%s\n','%  ---  Basis function Definitions  --------------------------------------------');
    fprintf(fid,'%s\n','!UseBasis={LinBasis|HEXAHEDRAL_ELEMENT|NODAL_LAGRANGE1|GAUSS_LEGENDRE8}');
    if Nodes_per_elem_dir==4
        fprintf(fid,'%s\n\n','!UseBasis={CubicBasis|HEXAHEDRAL_ELEMENT|NODAL_LAGRANGE3|GAUSS_LEGENDRE8}');
    elseif Nodes_per_elem_dir==3
        fprintf(fid,'%s\n\n','!UseBasis={QuadBasis|HEXAHEDRAL_ELEMENT|NODAL_LAGRANGE2|GAUSS_LEGENDRE8}');
    end
    fprintf(fid,'%s\n','%  ---  Topology Definitions  --------------------------------------------');
    if Nodes_per_elem_dir==4
        fprintf(fid,'%s\n','!DefTopology={TP3|Cubic_mesh|CubicBasis|*.B}');
    elseif Nodes_per_elem_dir==3
        fprintf(fid,'%s\n','!DefTopology={TP2|Quadratic_mesh|QuadBasis|*.B}');
    end
    fprintf(fid,'%s\n','!DefTopology={TP1|Lin_mesh|LinBasis}');
    if flag_disc_Fib==1
        if Nodes_per_elem_dir==4
            fprintf(fid,'%s\n','!DefTopology={TP3Discont| DisconFibMesh |CubicBasis}');
        elseif Nodes_per_elem_dir==3
            fprintf(fid,'%s\n','!DefTopology={TP2Discont| DisconFibMesh |QuadBasis}');
        end
    end
    fprintf(fid,'%s\n','%  ---  Interface Definitions  --------------------------------------------');
    if Nodes_per_elem_dir==4
        fprintf(fid,'%s\n','!DefInterface={OneToOne|TP3|TP1}');
    elseif Nodes_per_elem_dir==3
        fprintf(fid,'%s\n','!DefInterface={OneToOne|TP2|TP1}');
    end
    if flag_disc_Fib==1
        if Nodes_per_elem_dir==4
            fprintf(fid,'%s\n','!DefInterface={OneToOne|TP3Discont|TP3}');
            fprintf(fid,'%s\n\n','!DefInterface={OneToOne|TP3Discont|TP1}');
        elseif Nodes_per_elem_dir==3
            fprintf(fid,'%s\n','!DefInterface={OneToOne|TP2Discont|TP2}');
            fprintf(fid,'%s\n\n','!DefInterface={OneToOne|TP2Discont|TP1}');
        end
    end
    fprintf(fid,'%s\n','%  ---  Variable Pointer Definitions  --------------------------------------------');
    if Nodes_per_elem_dir==4
        fprintf(fid,'%s\n','!DefVariablePointer={Space|TP3|Cubic_mesh}');
    elseif Nodes_per_elem_dir==3
        fprintf(fid,'%s\n','!DefVariablePointer={Space|TP2|Quadratic_mesh}');
    end
    fprintf(fid,'   %d\n',3);
    if Nodes_per_elem_dir==4
        fprintf(fid,'%s\n','!DefVariablePointer={Disp|TP3|3}');
    elseif Nodes_per_elem_dir==3
        fprintf(fid,'%s\n','!DefVariablePointer={Disp|TP2|3}');
    end
    fprintf(fid,'%s\n','!DefVariablePointer={Pres|TP1|1}');
    if strcmp(material_law,'Guccione')==1 % || add any other laws that use fibers (holzapfel- ogden,...)
        if flag_disc_Fib==1
            if Nodes_per_elem_dir==4
                fprintf(fid,'%s\n','!DefVariablePointer={Fib|TP3Discont| DisconFibMesh}');
            elseif Nodes_per_elem_dir==3
                fprintf(fid,'%s\n','!DefVariablePointer={Fib|TP2Discont| DisconFibMesh}');
            end
        else
            if Nodes_per_elem_dir==4
                fprintf(fid,'%s\n','!DefVariablePointer={Fib|TP3| FIBERS.X}');
            elseif Nodes_per_elem_dir==3
                fprintf(fid,'%s\n','!DefVariablePointer={Fib|TP2| FIBERS.X}');
            end
        end
        fprintf(fid,'   %d\n',9);
    end
    coun_patches_TU=0; % initialisation
    if isempty(press_expr_endo)==0
        coun_patches_TU=coun_patches_TU+1; % afou exw LV pressure     
    end
    if isempty(press_expr_epi)==0
        coun_patches_TU=coun_patches_TU+1; % afou exw LV pressure     
    end
    if isempty(press_expr_sept)==0
        coun_patches_TU=coun_patches_TU+1; % afou exw LV pressure     
    end
    if flag_prescr_bound_base==1
        if Nodes_per_elem_dir==4
            fprintf(fid,'%s\n','!DefVariablePointer={BaseDisp|TP3|3}');
        elseif Nodes_per_elem_dir==3
            fprintf(fid,'%s\n','!DefVariablePointer={BaseDisp|TP2|3}');
        end
        fprintf(fid,'%s\n','!SetVariablePointer={BaseDisp|TEMPORAL_UPDATE_FILE|Tags_Base_Disp*.D}');
        coun_patches_TU=coun_patches_TU+1;
    end
    
    if flag_prescr_bound_epi==1
        if Nodes_per_elem_dir==4
            fprintf(fid,'%s\n','!DefVariablePointer={EpiDisp|TP3|3}');
        elseif Nodes_per_elem_dir==3
            fprintf(fid,'%s\n','!DefVariablePointer={EpiDisp|TP2|3}');
        end
        fprintf(fid,'%s\n','!SetVariablePointer={EpiDisp|TEMPORAL_UPDATE_FILE|Tags_Epi_Disp*.D}');
        coun_patches_TU=coun_patches_TU+1;
    end
    
    if flag_prescr_bound_endo==1
        if Nodes_per_elem_dir==4
            fprintf(fid,'%s\n','!DefVariablePointer={EndoDisp|TP3|3}');
        elseif Nodes_per_elem_dir==3
            fprintf(fid,'%s\n','!DefVariablePointer={EndoDisp|TP2|3}');
        end
        fprintf(fid,'%s\n','!SetVariablePointer={EndoDisp|TEMPORAL_UPDATE_FILE|Tags_Endo_Disp*.D}');
         coun_patches_TU=coun_patches_TU+1;
    end
    if flag_prescr_bound_apex_hole==1
        if Nodes_per_elem_dir==4
            fprintf(fid,'%s\n','!DefVariablePointer={ApexDisp|TP3|3}');
        elseif Nodes_per_elem_dir==3
            fprintf(fid,'%s\n','!DefVariablePointer={ApexDisp|TP2|3}');
        end
        fprintf(fid,'%s\n','!SetVariablePointer={ApexDisp|TEMPORAL_UPDATE_FILE|Tags_Apex_Disp*.D}');
         coun_patches_TU=coun_patches_TU+1;
    end
    if flag_prescr_APEX_NODE==1
        coun_patches_TU=coun_patches_TU+1;
    end
    
    if flag_fix_base==1
        coun_patches_TU=coun_patches_TU+1;
    end
    if flag_fix_apex_hole==1
        coun_patches_TU=coun_patches_TU+1;
    end
        
    fprintf(fid,'\n%s\n','%  ----------   Problem Definitions -----------------------------------------------');
    fprintf(fid,'%s\n','!DefProblem={solid|STEADY_QUASI_STATIC_ELASTICITY}');
    fprintf(fid,'%s\n','!UseVariablePointer={Space|Space}');
    fprintf(fid,'%s\n','!UseVariablePointer={Displacement|Disp}');
    fprintf(fid,'%s\n','!UseVariablePointer={Pressure|Pres}');
    if strcmp(material_law,'Guccione')==1
        fprintf(fid,'%s\n\n','!UseVariablePointer={Fibers|Fib}');
        fprintf(fid,'%s\n','!ConstitutiveLaw={orthofiber-classic}');
        C1=material_parameters(1); % C_1
        C2=material_parameters(2);  % b_ff
        C3=material_parameters(3); % b_sn- this must be it as it multiplies Ess^2, Enn^2 and Esn^2+Ens^2 whereas b_fs multiplies only off-diagonal terms that can be grouped into two: (Efs^2+Esf^2) and (Efn^2+Enf^2)
        C4=material_parameters(4); % bfs
        fprintf(fid,'%f  %f  %f  %f  %f  %f  %f\n',[C1,C2,C4,C4,C3,C3,C3].');
    elseif strcmp(material_law,'NeoHookean')==1
        fprintf(fid,'%s\n','!ConstitutiveLaw={NeoHookean}');
        fprintf(fid,'%f \n',material_parameters.');
    end
    fprintf(fid,'%s\n','!SetProblemMatrixCalculation={EVALUATE_EVERY_BUILD}');
    fprintf(fid,'%s\n','!Number-of-boundary-patches');
    fprintf(fid,'   %d\n\n',coun_patches_TU);
    fprintf(fid,'%s\n','');
    if isempty(press_expr_endo)==0
        fprintf(fid,'%s\n','!DefExpression={Pressure_endo}');
        fprintf(fid,'%s\n',['    ',press_expr_endo]);
    end
    if isempty(press_expr_epi)==0
        fprintf(fid,'%s\n','!DefExpression={Pressure_epi}');
        fprintf(fid,'%s\n',['    ',press_expr_epi]);
    end
    if isempty(press_expr_sept)==0
        fprintf(fid,'%s\n','!DefExpression={Pressure_sept}');
        fprintf(fid,'%s\n',['    ',press_expr_sept]);
    end
    fprintf(fid,'%s\n','!Boundary-patch-definitions');
    if isempty(press_expr_endo)==0
        fprintf(fid,'     %d',patchID_Endo);
        fprintf(fid,'%s\n','  Disp  SCALED_NORMAL Pressure_endo');  
    end
    if isempty(press_expr_epi)==0
        fprintf(fid,'     %d',patchID_Epi);
        fprintf(fid,'%s\n','  Disp  SCALED_NORMAL Pressure_epi');  
    end
    if isempty(press_expr_sept)==0
        fprintf(fid,'     %d',patchID_Sept);
        fprintf(fid,'%s\n','  Disp  SCALED_NORMAL Pressure_sept');  
    end
%     fprintf(fid,'%s\n','  Disp  SCALED_NORMAL infl.data');  % replaced this with expression. 
    if flag_prescr_bound_base==1
        fprintf(fid,'     %d',patchID_Base);
        fprintf(fid,'%s\n','  Disp  DIRICHLET  BaseDisp');
    end
    if flag_fix_base==1
        fprintf(fid,'     %d',patchID_Base);
        fprintf(fid,'%s\n','  Disp  DIRICHLET  0 0 0');
    end
    if flag_fix_apex_hole==1
        fprintf(fid,'     %d',patchID_Apex);
        fprintf(fid,'%s\n','  Disp  DIRICHLET  0 0 0');
    end
    if flag_prescr_bound_epi==1
        fprintf(fid,'     %d',patchID_Epi);
        fprintf(fid,'%s\n','  Disp  DIRICHLET  EpiDisp');
    end
    if flag_prescr_bound_endo==1
        fprintf(fid,'     %d',patchID_Endo);
        fprintf(fid,'%s\n','  Disp  DIRICHLET  EndoDisp');
    end
    if flag_prescr_bound_apex_hole==1
        fprintf(fid,'     %d',patchID_Apex);
        fprintf(fid,'%s\n','  Disp  DIRICHLET  ApexDisp');
    end
    if flag_prescr_APEX_NODE==1
        fprintf(fid, '%s\n','specific Disp.1 specific Prescr_Disp_1.FE');
        fprintf(fid, '%s\n','specific Disp.2 specific Prescr_Disp_2.FE');
        fprintf(fid, '%s\n','specific Disp.3 specific Prescr_Disp_3.FE');
    end
    if strcmp(material_law,'Guccione')==1
        fprintf(fid,'!VisualizeVariables={Space|Pres|Disp|Fib|TimeSequence[%d-%d-%d]}',[1,total_incr_curr_frame,1].');
    elseif strcmp(material_law,'NeoHookean')==1 % no Fib field to visualise here
        fprintf(fid,'!VisualizeVariables={Space|Pres|Disp|TimeSequence[%d-%d-%d]}',[1,total_incr_curr_frame,1].');
    end
    fclose(fid);