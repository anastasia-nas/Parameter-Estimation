function ksi_location_3D=find_ksi_location_of_local_node_ID_in_3D_element_cmgui_ordering(local_node_ID,Nodes_per_elem_dir)


% created on 2 Dec 2019 
% Input: 
% local_node_ID: a local node ID in surface element (from 1 to Nodes_per_elem_dir^3) 
% Output: 
% ksi_location_3D=[ksi1,ksi2,ksi3]; % row vector with the ksi_location of the
% node in a row vector where 0<=ksi(i)<=1,for i=1,2,3



local_node_ID_TU=local_node_ID-1; % turn ids into indices from 0 to Nodes_per_elem_dir^3-1 so that I can quickly use modular arithmetic
ksi_location_3D(1)=mod(local_node_ID_TU,Nodes_per_elem_dir)/(Nodes_per_elem_dir-1);
ksi_location_3D(2)=(mod(local_node_ID_TU,Nodes_per_elem_dir^2)-mod(local_node_ID_TU,Nodes_per_elem_dir))/Nodes_per_elem_dir/(Nodes_per_elem_dir-1);

ksi_location_3D(3)=((local_node_ID_TU-mod(local_node_ID_TU,Nodes_per_elem_dir^2))/Nodes_per_elem_dir^2)/(Nodes_per_elem_dir-1); 


