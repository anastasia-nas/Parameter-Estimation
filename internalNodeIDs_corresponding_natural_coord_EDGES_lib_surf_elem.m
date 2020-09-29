function [indicesEDGE_ksi1eq0_surf, indicesEDGE_ksi1eq1_surf, indicesEDGE_ksi2eq0_surf, indicesEDGE_ksi2eq1_surf]=internalNodeIDs_corresponding_natural_coord_EDGES_lib_surf_elem(Nodes_per_elem_dir)

% deixnei poia nodal ids anikoun se kathe akmi(EDGE) se natural coordinates (ta ids
% einai arithmimena me proteraiotita kata ksi1, meta ksi2 
% (diladi symfwna me cmgui kai oxi Cheart!!!) mesa se ena surface element
% (me synolika nodes : Nodes_per_elem_dir^2)

indicesEDGE_ksi1eq0_surf=([1:Nodes_per_elem_dir]-1)*Nodes_per_elem_dir+1;
indicesEDGE_ksi1eq1_surf=[1:Nodes_per_elem_dir]*Nodes_per_elem_dir;
indicesEDGE_ksi2eq0_surf=1:Nodes_per_elem_dir;
indicesEDGE_ksi2eq1_surf=Nodes_per_elem_dir*(Nodes_per_elem_dir-1)+[1:Nodes_per_elem_dir];