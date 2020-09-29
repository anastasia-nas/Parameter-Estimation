function [GPcoords GPweights]=Gauss_Points_In_2D_Element(n1,n2) % n1: #GP along ksi1, n2: #GP along ksi2

% akrivws idio me GaussPointsInElementNewSurf3D.m apla rename ekana gia
% clarity dioti den allazei kati an oloklirwnw se surface 3D element i an
% oloklirwnw ena 2D element... 

[xGP1 wGP1 xGP2 wGP2]=GaussPointsCoordsScaled(n1,n2);
for mn2=1:n2
    for mn1=1:n1
        GPcoords(mn1+n1*(mn2-1),:)=[xGP1(mn1),xGP2(mn2)];  % GPcoords: oisyntetagmenes kata ksi1, ksi2 tou kathe GP: row:#GP, 1st col:ksi1 coord, 2nd col:jksi2 coord
        GPweights(mn1+n1*(mn2-1),:)=[wGP1(mn1),wGP2(mn2)];  % GPweights: row: #GP, 1stcol: wight of this GP along ksi1, 2nd col: weight of GP along ksi2
    end
end



% morfi GPcoords:
% GPcoords=[ksi1GP1 ksi2GP1;
%           ksi1GP2 ksi2GP2;
%           ksi1GP3 ksi2GP3;
%           ksi1GP4 ksi2GP4]
%       
%       where: e.g. for n1=3, n2=2:
%       
%       ksi2
%       ^
%       |
%       |
%       |-------------------|
%       |                   |
%       |  4      5       6 |
%       |                   |
%       |                   |
%       |                   |
%       |                   |
%       |  1      2       3 |
%       |-------------------|------->ksi1
