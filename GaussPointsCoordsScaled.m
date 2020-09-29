function [xGP1 wGP1 xGP2 wGP2]=GaussPointsCoordsScaled(n1,n2) % n1: #GP along ksi1, n2: #GP along ksi2
[yx1, yw1]=GaussPointLibrary(n1); % oi theseis twn GP kata ksi1 an to ksi1 ekfrazotan sto [-1,1];
[yx2, yw2]=GaussPointLibrary(n2);
% assume a  ksi in [0,1] 
%use ksi1max, ksi1min instead of 1, 0 here if you ksi is not in [0,1], same for ksi2 
ksiAver=(0+1)/2;
scale=(1-0)/(1-(-1)); % by what the integral has  to be multiplied to account for the difference between ksi [-1,1] and anything else
xGP1=yx1*scale+ksiAver;  % xGP1, wGP1: positions and weights of the gauss points in the ksi1 direction same for ksi2 : SOS! xGP1 contains only 1dimensional info: only positions along ksi1
xGP2=yx2*scale+ksiAver; % i diataxi twn simeiwn einai apo to pio negative sta pio positive
wGP1=yw1*scale;  % ta vari mpainoun epi scale dioti skepsou tio gewmwtriki representation twn weights: einai se poso diastima aristera dexia antiproswpeyei to kathe GP-> ayta exoun ypologistei gia ena diastima 1-(-1)=2 enw twra mpainou se ena diastima 1-0=0 ara tha prepei na mpoun misa!!!!
wGP2=yw2*scale;
%xGP=[xGP1;xGP2]; % 1st row: along ksi1 axis, 2nd row: ksi2 axis size: 2(dimensions)*numberOfGaussPoints
%wGP=[wGP1;wGP2];
