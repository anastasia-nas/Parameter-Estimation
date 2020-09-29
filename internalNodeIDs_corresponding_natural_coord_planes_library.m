function [indices_plane_ksi1eq0,indices_plane_ksi1eq1,indices_plane_ksi2eq0,indices_plane_ksi2eq1,indices_plane_ksi3eq0,indices_plane_ksi3eq1,...
    indices_plane_diag_EDGE1_ksi1ksi2eq0_EDGE2_ksi1ksi2eq1,indices_plane_diag_EDGE1_ksi1eq1ksi2eq0_EDGE2_ksi1eq0ksi2eq1,...
    indices_plane_diag_EDGE1_ksi2ksi3eq0_EDGE2_ksi2ksi3eq1,indices_plane_diag_EDGE1_ksi2eq1ksi3eq0_EDGE2_ksi2eq0ksi3eq1,...
    indices_plane_diag_EDGE1_ksi3ksi1eq0_EDGE2_ksi3ksi1eq1,indices_plane_diag_EDGE1_ksi3eq0ksi1eq1_EDGE2_ksi3eq1ksi1eq0]=internalNodeIDs_corresponding_natural_coord_planes_library(Nodes_per_elem_dir)

% deixnei poia nodal ids anikoun sto kathe natural coordinate plane (ta ids
% einai arithmimena me proteraiotita kata ksi1, meta ksi2 kai meta ksi 3
% (diladi symfwna me cmgui kai oxi Cheart!!!)

indices_plane_ksi3eq0=1:Nodes_per_elem_dir^2;
indices_plane_ksi3eq1=(Nodes_per_elem_dir-1)*Nodes_per_elem_dir^2+[1:Nodes_per_elem_dir^2];
[meshX,meshY]=meshgrid([1:Nodes_per_elem_dir],[1:Nodes_per_elem_dir]);
indices_plane_ksi1eq0=reshape(Nodes_per_elem_dir^2*(meshX-1)+(meshY-1)*Nodes_per_elem_dir,[1,size(meshX,1)*size(meshX,2)])+1;
indices_plane_ksi1eq1=reshape(Nodes_per_elem_dir^2*(meshX-1)+(meshY-1)*Nodes_per_elem_dir,[1,size(meshX,1)*size(meshX,2)])+Nodes_per_elem_dir;
indices_plane_ksi2eq0=reshape(Nodes_per_elem_dir^2*(meshX-1)+meshY,[1,numel(meshX)]);
indices_plane_ksi2eq1=reshape(Nodes_per_elem_dir^2*(meshX-1)+(Nodes_per_elem_dir*(Nodes_per_elem_dir-1))+meshY,[1,numel(meshX)]);
% episis gia to sygkekrimeno in silico mesh to XZ plane
% mporei na to petyxainei mesa sto psaxno alla to vriskei
% diagwnia me 2 tropous : eite sto epipedo apo akmi ksi1=0, ksi3=0 pros
% akmi ksi=1 ksi3=1 (indices_diag_ksi3ksi1eq0_ksi3ksi1eq1) (zwgrafise 1 kyvo me to pws paei to
% ksi1, ksi2, ksi3 gia na katalaveis); eite sto epipedo apo
% akmi ksi1=0 ksi3=1 ews akmi ksi1=1 ksi3=0
% (indices_diag_ksi3eq0ksi1eq1_ksi3eq1ksi1eq0) klp gia ta
% 6 diagwnia pou enwnoun antidiametrikes akmes

% don't touch this: einai swsto!
indices_plane_diag_EDGE1_ksi3ksi1eq0_EDGE2_ksi3ksi1eq1=reshape(Nodes_per_elem_dir^2*(meshX-1)+Nodes_per_elem_dir*(meshY-1)+meshX,[1,numel(meshX)]);
indices_plane_diag_EDGE1_ksi3eq0ksi1eq1_EDGE2_ksi3eq1ksi1eq0=reshape(Nodes_per_elem_dir^2*(meshX-1)+Nodes_per_elem_dir*(meshY-1)+(sortrows(meshY,-1)).',[1,numel(meshX)]); % to sort (sortrows(meshY,-1)-1).' einai to meshX pou apla exei antistrafei i seira (paw apo to 4 sto 1 se kathe col, anti gia apo to 1 sto 4)
indices_plane_diag_EDGE1_ksi2ksi3eq0_EDGE2_ksi2ksi3eq1=reshape(Nodes_per_elem_dir^2*(meshX-1)+Nodes_per_elem_dir*(meshX-1)+meshY,[1,numel(meshX)]);
indices_plane_diag_EDGE1_ksi2eq1ksi3eq0_EDGE2_ksi2eq0ksi3eq1=reshape(Nodes_per_elem_dir^2*(meshX-1)+Nodes_per_elem_dir*(sortrows(meshY,-1).'-1)+meshY,[1,numel(meshX)]); 
indices_plane_diag_EDGE1_ksi1ksi2eq0_EDGE2_ksi1ksi2eq1=reshape(Nodes_per_elem_dir^2*(meshX-1)+Nodes_per_elem_dir*(meshY-1)+meshY,[1,numel(meshX)]);
indices_plane_diag_EDGE1_ksi1eq1ksi2eq0_EDGE2_ksi1eq0ksi2eq1=reshape(Nodes_per_elem_dir^2*(meshX-1)+Nodes_per_elem_dir*(meshY-1)+sortrows(meshY,-1),[1,numel(meshX)]);