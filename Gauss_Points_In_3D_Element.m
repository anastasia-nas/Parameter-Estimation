function [GPcoords GPweights]=Gauss_Points_In_3D_Element(n1,n2,n3) % n1: #GP along ksi1, n2: #GP along ksi2

% einai oloidio me to GaussPointsInElementNew3D.m apla to kanw rename gia
% na tairiazei me ta Gauss_Points_In_2D_Element.m kai
% Gauss_Points_In_1D_Element.m
[xGP1 wGP1 xGP2 wGP2 xGP3 wGP3]=GaussPointsCoordsScaled3D(n1,n2,n3); % ta XGP ta lew X alla einai ksi sti pragmatikotita
for mn3=1:n3
    for mn2=1:n2
        for mn1=1:n1
            GPcoords(mn1+(mn2-1)*n1+(mn3-1)*n2*n1,:)=[xGP1(mn1),xGP2(mn2),xGP3(mn3)];  % GPcoords: oisyntetagmenes kata ksi1, ksi2 tou kathe GP: row:#GP, 1st col:ksi1 coord, 2nd col:jksi2 coord
            GPweights(mn1+(mn2-1)*n1+(mn3-1)*n2*n1,:)=[wGP1(mn1),wGP2(mn2),wGP3(mn3)];  % GPweights: row: #GP, 1stcol: wight of this GP along ksi1, 2nd col: weight of GP along ksi2
        end
    end
end

%  sto 3D ta GPpoints arithmountai ws : prwta exantleis ti ksi1, meta to
%  exantleis kata ksi2 kai telos kata ksi3



% morfi GPcoords:
% GPcoords=[ksi1GP1 ksi2GP1 ksi3GP1;
%           ksi1GP2 ksi2GP2 ksi3GP2;
%           ksi1GP3 ksi2GP3 ksi3GP3;
%           ksi1GP4 ksi2GP4 ksi3GP4;
%           ........              ];
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
