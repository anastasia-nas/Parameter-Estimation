function [H_0_0,H_0_1,H_1_0,H_1_1]=Hermite_basis_functions_1D(ksi)

% created by Anastasia on 23-05-2019
%-------------------------------------------------------------------------------------------------------------------------
% Purpose: 
%-------------------------------------------------------------------------------------------------------------------------
% just writes down hermite basis functions in 1D - for 3D you just need
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
% Cubic Lagrange meshes in Cheart format option to print.
%-------------------------------------------------------------------------------------------------------------------------
% SOS points:
%-------------------------------------------------------------------------------------------------------------------------
% ksi belongs to [0,1]
%=============================================================================================================================================================

H_0_0=1-3*ksi^2+2*ksi^3;

H_0_1=ksi^2*(3-2*ksi);

H_1_0=ksi*(ksi-1)^2;

H_1_1=ksi^2*(ksi-1);
