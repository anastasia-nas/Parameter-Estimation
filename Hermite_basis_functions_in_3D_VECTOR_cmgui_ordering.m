function Shape_function_vector_cmgui_ordering=Hermite_basis_functions_in_3D_VECTOR_cmgui_ordering(ksi1,ksi2,ksi3)

% created by Anastasia on 23-05-2019
%-------------------------------------------------------------------------------------------------------------------------
% Purpose: 
%-------------------------------------------------------------------------------------------------------------------------
% cubic hermite basis functions in 3D you just need
% products of these basic basis function in 1D --see paper by Smith et al
% 2004 where the CH interpolation is explained explicity
%-------------------------------------------------------------------------------------------------------------------------
% Notes on methodology: 
%-------------------------------------------------------------------------------------------------------------------------
%  follows Smith 2004, Multiscale computational modelling of the heart
% in 1D the cubic hermite interpolation of variable u:
% u = u(0)*H_0_0 + u(1)*H_0_1 + du/dksi(0)*H_1_0 + du/dksi(1)*H_1_1
% hence H_0_0:multiplies variable value u at node with ksi=0
%       H_0_1:multiplies variable value u at node with ksi=1
%       H_1_0:multiplies variable gradient du/dksi at node with ksi=0
%       H_1_1:multiplies variable gradient du/dksi at node with ksi=1
%-------------------------------------------------------------------------------------------------------------------------
% Input:
%-------------------------------------------------------------------------------------------------------------------------
% local coordinate value in 1D. For use in 3D you will use the local
% coordinate ksi1,ksi2 or ksi3 depending on the local coordinate the basis
% function corresponds to 


%-------------------------------------------------------------------------------------------------------------------------
% Output:
%-------------------------------------------------------------------------------------------------------------------------
% Shape_function_vector_cmgui_ordering : row vector size (1,64) with interpolation values at each node of element following cmgui
% ordering for values and derivatives (first priority the values inside the node according to cmgui:
% (u,du/dksi1,du/dksi2,du^2/dksi1dksi2,du/dksi3,du^2/dksi1dksi3,du^2/dksi2dksi3,du3/dksi1dksi2dksi3)
% and then from node to node following the ksi1,ksi2,ksi3 ordering.
%-------------------------------------------------------------------------------------------------------------------------
% SOS points:
%-------------------------------------------------------------------------------------------------------------------------
% ksi belongs to [0,1]
% nodal values follow cmgui ordering
%=============================================================================================================================================================

[H_0_0_ksi1,H_0_1_ksi1,H_1_0_ksi1,H_1_1_ksi1]=Hermite_basis_functions_1D(ksi1);

[H_0_0_ksi2,H_0_1_ksi2,H_1_0_ksi2,H_1_1_ksi2]=Hermite_basis_functions_1D(ksi2);

[H_0_0_ksi3,H_0_1_ksi3,H_1_0_ksi3,H_1_1_ksi3]=Hermite_basis_functions_1D(ksi3);

H_mat=[H_0_0_ksi1,H_0_1_ksi1,H_1_0_ksi1,H_1_1_ksi1;...
       H_0_0_ksi2,H_0_1_ksi2,H_1_0_ksi2,H_1_1_ksi2;...
       H_0_0_ksi3,H_0_1_ksi3,H_1_0_ksi3,H_1_1_ksi3]; % row: ksi1,ksi2,ksi3 direction, col:1-2 value at node 1 and 2 and 3-4 derivative at node 1 and 2 

Shape_function_vector_cmgui_ordering=zeros(1,64);
% Internal node ordering first along ksi1, then along ksi2, and then ksi3
for N1=1:2
    for N2=1:2
        for N3=1:2
            internal_node_ID=(N3-1)*2*2 + (N2-1)*2 + N1;
            
            nodal_ksi_Values=[N1-1,N2-1,N3-1]; % ksi1, ksi2, ksi3 values at each corner node of Hermite Element:
%          ordering followed:   u,du/dksi1,du/dksi2,du^2/dksi1dksi2,du/dksi3,du^2/dksi1dksi3,du^2/dksi2dksi3,du3/dksi1dksi2dksi3

            Shape_function_vector_cmgui_ordering((internal_node_ID-1)*8+1:internal_node_ID*8)=[H_mat(1,N1)*H_mat(2,N2)*H_mat(3,N3),... %  u
                                                                                               H_mat(1,2+N1)*H_mat(2,N2)*H_mat(3,N3),... % du/dksi1
                                                                                               H_mat(1,N1)*H_mat(2,2+N2)*H_mat(3,N3),... % du/dksi2
                                                                                               H_mat(1,2+N1)*H_mat(2,2+N2)*H_mat(3,N3),... % du^2/dksi1dksi2
                                                                                               H_mat(1,N1)*H_mat(2,N2)*H_mat(3,2+N3),...  % du/dksi3
                                                                                               H_mat(1,2+N1)*H_mat(2,N2)*H_mat(3,2+N3),... % du^2/dksi1dksi3
                                                                                               H_mat(1,N1)*H_mat(2,2+N2)*H_mat(3,2+N3),... % du^2/dksi2dksi3
                                                                                               H_mat(1,2+N1)*H_mat(2,2+N2)*H_mat(3,2+N3)]; % du3/dksi1dksi2dksi3                                                                                           
        end
    end
end
                                                                                               




