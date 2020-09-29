function create_Cheart_Bfile_from_Boundaries_mat(Lagrange_mesh_Boundaries,Nodes_per_elem_dir,Sim_prep_folder,B_file_name)

% created by Anastasia on 26-02-2020
%-------------------------------------------------------------------------------------------------------------------------
% Purpose: 
%-------------------------------------------------------------------------------------------------------------------------
% a function for printing a Boundaries matrix into a .B file for use in
% simulations or data processing following Cheart style.

%-------------------------------------------------------------------------------------------------------------------------
% Notes on methodology:
%-------------------------------------------------------------------------------------------------------------------------

%-------------------------------------------------------------------------------------------------------------------------
% Input:
%-------------------------------------------------------------------------------------------------------------------------
% Lagrange_mesh_Boundaries: boundaries mat in format: [element_ID,nodes on boundary (Nodes_per_elem_dir^2 in size), patch_ID]
% Nodes_per_elem_dir: indicates Lagrange interpolation order (for Cubic mesh =4)
% Sim_prep_folder: folder where the B.file will be printed
% B_file_name: the full filename you want to use (e.g. B_file_name='Cubic_mesh_FE.B') --  the file ending in ".B" eg Cubic_mesh_FE.B (That's the default in Cheart)

%-------------------------------------------------------------------------------------------------------------------------
% Output:
%-------------------------------------------------------------------------------------------------------------------------

%-------------------------------------------------------------------------------------------------------------------------
% SOS points:
%-------------------------------------------------------------------------------------------------------------------------

% % % B_file_name=[Lagrange_meshes_interp_order,'_mesh_FE.B']

% Cubic_mesh_FE.B
boundary_format_spec=[];
for n_n=1:(Nodes_per_elem_dir^2+2)
    boundary_format_spec=[boundary_format_spec,'%d     '];
end
boundary_format_spec=[boundary_format_spec,'\n'];

fid=fopen([Sim_prep_folder,'/',B_file_name],'w');
fprintf(fid,'%d\n',size(Lagrange_mesh_Boundaries,1));
fprintf(fid,boundary_format_spec,Lagrange_mesh_Boundaries.');
fclose(fid);
