function create_Cheart_NBTfile_from_Neighbours_mat(Lagrange_mesh_Neighbours,Sim_prep_folder,NBT_file_name)

% created by Anastasia on 26-02-2020
%-------------------------------------------------------------------------------------------------------------------------
% Purpose: 
%-------------------------------------------------------------------------------------------------------------------------
% a function for printing a Neighbours matrix into a .NBT file for use in
% simulations or data processing following Cheart style.

%-------------------------------------------------------------------------------------------------------------------------
% Notes on methodology:
%-------------------------------------------------------------------------------------------------------------------------

%-------------------------------------------------------------------------------------------------------------------------
% Input:
%-------------------------------------------------------------------------------------------------------------------------
% Lagrange_mesh_Neighbours: Neighoburs mat in format: [-ksi3,-ksi2,-ksi1,ksi,ksi2,ksi3] where each row corresponds to an element
% Nodes_per_elem_dir: indicates Lagrange interpolation order (for Cubic mesh =4)
% Sim_prep_folder: folder where the B.file will be printed
% NBT_file_name: the full filename you want to use (e.g. NBT_file_name='Cubic_mesh_FE.NBT') --  the file ending in ".NBT" eg Cubic_mesh_FE.NBT (That's the default in Cheart)

%-------------------------------------------------------------------------------------------------------------------------
% Output:
%-------------------------------------------------------------------------------------------------------------------------

%-------------------------------------------------------------------------------------------------------------------------
% SOS points:
%-------------------------------------------------------------------------------------------------------------------------

% NBT_file_name=[Lagrange_meshes_interp_order,'_mesh_FE.NBT']

% Cubic_mesh_FE.NBT
fid=fopen([Sim_prep_folder,'/',NBT_file_name],'w');
fprintf(fid,'%d\n',size(Lagrange_mesh_Neighbours,1));
fprintf(fid,'%d   %d   %d   %d   %d   %d\n',Lagrange_mesh_Neighbours.');
fclose(fid);