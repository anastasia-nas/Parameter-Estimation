function [indicesEDGE_ksi1ksi2eq0, indicesEDGE_ksi1ksi2eq1, indicesEDGE_ksi1eq1ksi2eq0, indicesEDGE_ksi1eq0ksi2eq1,...
    indicesEDGE_ksi2ksi3eq0, indicesEDGE_ksi2ksi3eq1, indicesEDGE_ksi2eq1ksi3eq0, indicesEDGE_ksi2eq0ksi3eq1,...
    indicesEDGE_ksi3ksi1eq0, indicesEDGE_ksi3ksi1eq1, indicesEDGE_ksi3eq0ksi1eq1, indicesEDGE_ksi3eq1ksi1eq0]=internalNodeIDs_corresponding_natural_coord_EDGES_library(Nodes_per_elem_dir)

% deixnei poia nodal ids anikoun se kathe akmi(EDGE) se natural coordinates (ta ids
% einai arithmimena me proteraiotita kata ksi1, meta ksi2 kai meta ksi 3
% (diladi symfwna me cmgui kai oxi Cheart!!!)

indicesEDGE_ksi1ksi2eq0=Nodes_per_elem_dir^2*([1:Nodes_per_elem_dir]-1)+1;
indicesEDGE_ksi1ksi2eq1=Nodes_per_elem_dir^2*[1:Nodes_per_elem_dir];
indicesEDGE_ksi1eq1ksi2eq0=Nodes_per_elem_dir^2*([1:Nodes_per_elem_dir]-1)+Nodes_per_elem_dir;
indicesEDGE_ksi1eq0ksi2eq1=Nodes_per_elem_dir^2*[1:Nodes_per_elem_dir]-(Nodes_per_elem_dir-1);
indicesEDGE_ksi2ksi3eq0=1:Nodes_per_elem_dir;
indicesEDGE_ksi2ksi3eq1=Nodes_per_elem_dir^3-(Nodes_per_elem_dir-[1:Nodes_per_elem_dir]);
indicesEDGE_ksi2eq1ksi3eq0=Nodes_per_elem_dir^2-(Nodes_per_elem_dir-[1:Nodes_per_elem_dir]);
indicesEDGE_ksi2eq0ksi3eq1=(Nodes_per_elem_dir-1)*Nodes_per_elem_dir^2+(1:Nodes_per_elem_dir);
indicesEDGE_ksi3ksi1eq0=Nodes_per_elem_dir*((1:Nodes_per_elem_dir)-1)+1;
indicesEDGE_ksi3ksi1eq1=Nodes_per_elem_dir^2*(Nodes_per_elem_dir-1)+Nodes_per_elem_dir*(1:Nodes_per_elem_dir);
indicesEDGE_ksi3eq0ksi1eq1=Nodes_per_elem_dir*(1:Nodes_per_elem_dir);
indicesEDGE_ksi3eq1ksi1eq0=Nodes_per_elem_dir^2*(Nodes_per_elem_dir-1)+Nodes_per_elem_dir*((1:Nodes_per_elem_dir)-1)+1;