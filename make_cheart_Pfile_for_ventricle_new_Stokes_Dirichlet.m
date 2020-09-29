function make_cheart_Pfile_for_ventricle_new_Stokes_Dirichlet(sim_folder_dir,incr_num,endo_elements_new_patch_ID,endo_patch_ID,epi_patch_ID,base_patch_ID,apex_patch_ID,flag_Dirichlet_Endo,flag_Dirichlet_Epi,Nodes_per_elem_dir)

% newer version of script created 15 Jan 2020 to account for Quadratic
% interpolation for disp/Space/fibers.

% from function make_cheart_Pfile_for_ventricle(sim_folder_dir,incr_num,Disp_interpl,Pres_intrpl,Fib_intrpl,mat_law,mat_params,patch_Endo,patch_Epi,patch_Base,patch_Apex,flag_base,flag_epi,flag_endo,flag_apex,flag_disc_Fib)

% sim_folder_dir: full path to simulation folder
% incr_num: total number of simulation increments

% patch_Endo: the number ID of endocardial patch, same for patch_Epi,
% patch_Base (usually it's :patch_Endo=1, patch_Epi=2, patch_Base=3
% patch_Apex=4)

% flag_Dirichlet_Endo: =1 if you create prescr disps file for endo - else=0
% flag_Dirichlet_Epi: =1 if you create prescr disps for epi - else=0
    
if Nodes_per_elem_dir==4
       interp_order_str='Cubic';
elseif Nodes_per_elem_dir==3
       interp_order_str='Quad';
end
coun_boundary_patches=2; % because you fix apex and base
%%
%     fprintf(fid,'     %d',patch_Apex);
fid=fopen([sim_folder_dir,'/Ellipse_hole_Stokes_Dirichlet.P'],'w');
fprintf(fid,'%s\n','!DefSolverGroup={SGroup|TimeStepping|SMatrix}');
fprintf(fid,'%s\n','!DefSolverSubGroup={SGroup|SOLVER_SEQUENTIAL|SMatrix}');
if flag_Dirichlet_Epi==1
    fprintf(fid,'%s\n','!SetSolverGroup={SGroup|AddVariables|EpiVel}');
end
if flag_Dirichlet_Endo==1
    fprintf(fid,'%s\n','!SetSolverGroup={SGroup|AddVariables|EndoVel}');
end
fprintf(fid,'%s\n','!DefTimeStepScheme={TimeStepping}');
fprintf(fid,'%d  %d  %d\n',[1,incr_num,1]);
fprintf(fid,'%s\n','');
fprintf(fid,'%s\n','!DefSolverMatrix={SMatrix|SOLVER_MUMPS|STEADY_STOKES_PROBLEM}');
fprintf(fid,'%s\n','');
fprintf(fid,'%s\n','%  ---  Topologies & Basis sets  --------------------------------------------');
fprintf(fid,'%s\n','!UseBasis={LinBasis|HEXAHEDRAL_ELEMENT|NODAL_LAGRANGE1|GAUSS_LEGENDRE8}');
str=['!UseBasis={',interp_order_str,'Basis|HEXAHEDRAL_ELEMENT|NODAL_LAGRANGE',num2str(Nodes_per_elem_dir-1),'|GAUSS_LEGENDRE8}'];
fprintf(fid,'%s\n',str);
str=['!DefTopology={TP',num2str(Nodes_per_elem_dir-1),'|',interp_order_str,'_mesh|',interp_order_str,'Basis|*.B}']; % !DefTopology={TP2|Quad_mesh|QuadBasis|*.B}
fprintf(fid,'%s\n',str);  
fprintf(fid,'%s\n','!DefTopology={TP1|Lin_mesh|LinBasis}');
str=['!DefInterface={OneToOne|TP',num2str(Nodes_per_elem_dir-1),'|TP1}'];
fprintf(fid,'%s\n',str); %!DefInterface={OneToOne|TP2|TP1}
fprintf(fid,'%s\n','% --- Variables ------------------------------------------------');
str=['!DefVariablePointer={Space|TP',num2str(Nodes_per_elem_dir-1),'|',interp_order_str,'_mesh}']; % !DefVariablePointer={Space|TP2|Quad_mesh}
fprintf(fid,'%s\n',str);
fprintf(fid,'%d\n',3);
str=['!DefVariablePointer={Vel|TP',num2str(Nodes_per_elem_dir-1),'|3}'];
fprintf(fid,'%s\n',str);
fprintf(fid,'%s\n','!DefVariablePointer={Pres|TP1|1}');
if flag_Dirichlet_Epi==1
    str=['!DefVariablePointer={EpiVel|TP',num2str(Nodes_per_elem_dir-1),'|3}'];
    fprintf(fid,'%s\n',str);
    fprintf(fid,'%s\n','!SetVariablePointer={EpiVel|TEMPORAL_UPDATE_FILE|Prescr_Vel_Epi*.D}');
    coun_boundary_patches=coun_boundary_patches+1;
end
if flag_Dirichlet_Endo==1
    str=['!DefVariablePointer={EndoVel|TP',num2str(Nodes_per_elem_dir-1),'|3}'];
    fprintf(fid,'%s\n',str);
    fprintf(fid,'%s\n','!SetVariablePointer={EndoVel|TEMPORAL_UPDATE_FILE|Prescr_Vel_Endo*.D}');
    coun_boundary_patches=coun_boundary_patches+1;
end
fprintf(fid,'%s\n','');
fprintf(fid,'%s\n','% %!SetProblemMatrixCalculation={EVALUATE_EVERY_BUILD} % do not think I need this here');
fprintf(fid,'%s\n','');
fprintf(fid,'%s\n','!Number-of-boundary-patches');
fprintf(fid,'%d\n',coun_boundary_patches);
fprintf(fid,'%s\n','');
fprintf(fid,'%s\n','!DefProblem={STEADY_STOKES_PROBLEM|STEADY_STOKES_FLOW}');
fprintf(fid,'%s\n','!UseVariablePointer={Space   |Space}');
fprintf(fid,'%s\n','!UseVariablePointer={Velocity|Vel}');
fprintf(fid,'%s\n','!UseVariablePointer={Pressure|Pres}');
fprintf(fid,'%s\n','!True-Navier-Poisson');
fprintf(fid,'%s\n','!Viscosity={200.0}'); % leave it be for now
fprintf(fid,'%s\n','!Density={0.1}');
fprintf(fid,'%s\n','!boundary-patch-definitions');
str=['    ',num2str(base_patch_ID),'  Vel    DIRICHLET 0.0 0.0 0.0 % zero at base'];
fprintf(fid,'%s\n',str);
str=['    ',num2str(apex_patch_ID),'  Vel    DIRICHLET 0.0 0.0 0.0 % zero at apex hole'];
fprintf(fid,'%s\n',str);
if flag_Dirichlet_Epi==1
    str=['    ',num2str(epi_patch_ID),'  Vel    DIRICHLET  EpiVel'];
    fprintf(fid,'%s\n',str);
end
if flag_Dirichlet_Endo==1
    str=['    ',num2str(endo_elements_new_patch_ID),'  Vel    DIRICHLET  EndoVel'];
    fprintf(fid,'%s\n',str);
end
fprintf(fid,'%s\n','');

str=['!VisualizeVariables={Space|Pres|Vel|TimeSequence[1-',num2str(incr_num),'-1]}']; % !VisualizeVariables={Space|Pres|Vel|TimeSequence[1-30-1]}
fprintf(fid,'%s\n',str);
fclose(fid);

