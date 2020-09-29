function [xGP1, wGP1, xGP2, wGP2, xGP3, wGP3]=GaussPointsCoordsScaled3D(n1,n2,n3) % n1: #GP along ksi1, n2: #GP along ksi2
[yx1, yw1]=GaussPointLibrary(n1); % oi theseis twn GP kata ksi1 an to ksi1 ekfrazotan sto [-1,1];
[yx2, yw2]=GaussPointLibrary(n2);
[yx3, yw3]=GaussPointLibrary(n3);
% assume a  ksi in [0,1] 
%use ksi1max, ksi1min instead of 1, 0 here if you ksi is not in [0,1], same for ksi2 
ksiAver=(0+1)/2;
scale=(1-0)/(1-(-1)); % by what the integral has  to be multiplied to account for the difference between ksi [-1,1] and anything else-> dioti to mikos sto [-1,1] einai diplo apo to [0,1]
xGP1=yx1*scale+ksiAver;  % xGP1, wGP1: positions and weights of the gauss points in the ksi1 direction same for ksi2 : SOS! xGP1 contains only 1dimensional info: only positions along ksi1
xGP2=yx2*scale+ksiAver; % i diataxi twn simeiwn einai apo to pio negative sta pio positive
xGP3=yx3*scale+ksiAver;
wGP1=yw1*scale;  % ta vari mpainoun epi scale dioti skepsou tio gewmwtriki representation twn weights: einai se poso diastima aristera dexia antiproswpeyei to kathe GP-> ayta exoun ypologistei gia ena diastima 1-(-1)=2 enw twra mpainou se ena diastima 1-0=0 ara tha prepei na mpoun misa!!!!
wGP2=yw2*scale;
wGP3=yw3*scale;
%xGP=[xGP1;xGP2]; % 1st row: along ksi1 axis, 2nd row: ksi2 axis size: 2(dimensions)*numberOfGaussPoints
%wGP=[wGP1;wGP2];


% o olos logos pou vazw ta xGP1,xGP2,xGP3 is to allow for diffrent number
% of GPs along each direction, but if n1=n2=n3 then xGP1,xGP2,xGP3 are
% exactly the same