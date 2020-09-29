function create_exnode_exelem_files_for_Space_var_and_Lagrange_interp(GroupName,exelem_filename,exnode_filename,save_to_folder,Nodes,Elements_in_cmgui_ordering,Nodes_per_elem_dir)

% created by Anastasia --
% modified from create_exnode_exelem_files_for_Space_var_and_CL.m on 1 Oct
% 2019 to also print quadratic Lagrange meshes.
%exnode_filename: string in filename.exnode format
% exelem_filename: string in filename.exelem format
% create exnode exelem files of group name  'GroupName'
% Nodes: 3 col matrix,where no of row is Node ID and in each col is the
% node coordinate in x,y,z
% Elements_in_cmgui_ordering: matrix where in ith row is element i and 
% jth col is the id of the jth internal Node numbered according to  cmgui (no coorner
% priority like in Cheart)

if Nodes_per_elem_dir==3
    interp_str='q';
elseif Nodes_per_elem_dir==4
    interp_str='c';
elseif Nodes_per_elem_dir==2
    interp_str='l';
else
    disp('error! this only works up to cubic order of Lagrange interpolation -- check your interpolation scheme--quiting!')
    return;
end

fid=fopen([save_to_folder,'/',exnode_filename],'w');
fprintf(fid,' Group name: %s\n',GroupName);
fprintf(fid,'%s\n',' #Fields=1');
fprintf(fid,'%s\n',' 1) Space, coordinate, rectangular cartesian, #Components=3');
fprintf(fid,'%s\n','   1.  Value index=1, #Derivatives=0');
fprintf(fid,'%s\n','   2.  Value index=2, #Derivatives=0');
fprintf(fid,'%s\n','   3.  Value index=3, #Derivatives=0');
for n=1:size(Nodes,1)
    fprintf(fid,' Node:           %d\n',n);
    for m=1:3
        fprintf(fid,'  %25.20g\n',Nodes(n,m));
    end
end
fclose(fid);

fid=fopen([save_to_folder,'/',exelem_filename],'w');
fprintf(fid,' Group name: %s\n',GroupName);
fprintf(fid,'%s\n',' Shape.  Dimension=3');
fprintf(fid,'%s\n',' #Scale factor sets= 0');
fprintf(fid,' #Nodes=%d\n',Nodes_per_elem_dir^3); 
fprintf(fid,'%s\n',' #Fields=1');
fprintf(fid,'%s\n',' 1) Space, coordinate, rectangular cartesian, #Components=3');
for n_dim=1:3
%     fprintf(fid,'  %d.  c.Lagrange*c.Lagrange*c.Lagrange, no modify, standard node based.\n',n_dim);
    fprintf(fid,'%s\n',['  ',num2str(n_dim),'.  ',interp_str,'.Lagrange*',interp_str,'.Lagrange*',interp_str,'.Lagrange, no modify, standard node based.']);
    fprintf(fid,'     #Nodes= %d\n', Nodes_per_elem_dir^3);
    for n_nod=1:Nodes_per_elem_dir^3
        fprintf(fid,'     %d.  #Values=1\n',n_nod);
        fprintf(fid,'%s\n','       Value indices:     1');
        fprintf(fid,'%s\n','       Scale factor indices:    0');
    end
end
for n_el=1:size(Elements_in_cmgui_ordering,1)
    fprintf(fid,'Element:           %d  0  0\n',n_el);
    fprintf(fid,'%s\n','    Nodes:');
    fprintf(fid,'            %d',Elements_in_cmgui_ordering(n_el,:).');
    fprintf(fid,'\n');
end
fclose(fid);