function create_Cheart_Tfile_from_Elements_mat(Lagrange_mesh_Elements,Lagrange_mesh_Nodes,Nodes_per_elem_dir,Sim_prep_folder,T_file_name)

% created by Anastasia on 26-02-2020
%-------------------------------------------------------------------------------------------------------------------------
% Purpose: 
%-------------------------------------------------------------------------------------------------------------------------
% a function for printing a Element connectivity matrix into a .T file for use in
% simulations or data processing following Cheart style.

%-------------------------------------------------------------------------------------------------------------------------
% Notes on methodology:
%-------------------------------------------------------------------------------------------------------------------------

%-------------------------------------------------------------------------------------------------------------------------
% Input:
%-------------------------------------------------------------------------------------------------------------------------
% Lagrange_mesh_Elements: elements connectivity matrix with size(#number of elements, Nodes_per_elem_dir^3)
% Lagrange_mesh_Nodes : Nodes matrix each row: node ID, each col: coordinate
% Nodes_per_elem_dir: indicates Lagrange interpolation order (for Cubic mesh =4)
% Sim_prep_folder: folder where the B.file will be printed
% T_file_name: the full filename you want to use (e.g. T_file_name='Cubic_mesh_FE.T') --  the file ending in ".T" eg Cubic_mesh_FE.T (That's the default in Cheart)

%-------------------------------------------------------------------------------------------------------------------------
% Output:
%-------------------------------------------------------------------------------------------------------------------------

%-------------------------------------------------------------------------------------------------------------------------
% SOS points:
%-------------------------------------------------------------------------------------------------------------------------
% % % Lagrange_mesh_Elements=Linear_mesh_Elements;
% % % Lagrange_mesh_Nodes=Linear_mesh_Nodes;
% % % T_file_name='Lin_mesh_FE.T';
element_format_spec=[];
for n_n=1:Nodes_per_elem_dir^3
    element_format_spec=[element_format_spec,'%d     '];
end
element_format_spec=[element_format_spec,'\n'];
fid=fopen([Sim_prep_folder,'/',T_file_name],'w');
fprintf(fid,'%d %d\n',[size(Lagrange_mesh_Elements,1),size(Lagrange_mesh_Nodes,1)].');
fprintf(fid,element_format_spec,Lagrange_mesh_Elements.');
fclose(fid);