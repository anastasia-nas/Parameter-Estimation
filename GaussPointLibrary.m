function [x,w]=GaussPointLibrary(n)  % yx: thesi GP sto [-1,1],  yw: weight GP sti thesi yx, n: total number of GP n=[n1,n2], n1-> #GP alomg ksi1, n2# GP along ksi2
    if  n==1
        x=0;
        w=2;
    end
    if  n==2
        x=[-1/sqrt(3),1/sqrt(3)];
        w=[1,1];
    end
    if  n==3
        x=[-sqrt(0.6),0,sqrt(0.6)];
        w=[5/9,8/9,5/9];
    end
    if  n==4
        x=[-0.861136311594053, -0.339981043584856, 0.339981043584856, 0.861136311594053];
        w=[0.347854845137454, 0.652145154862546, 0.652145154862546, 0.347854845137454];
    end
    if  n==5
        x=[-0.906179845938664,-0.538469310105683,0,0.538469310105683,0.906179845938664];
        w=[0.236926885056189, 0.478628670499366, 0.568888888888889, 0.478628670499366, 0.236926885056189];
    end
    if  n==6
        x=[-0.932469514203152, -0.661209386466265, -0.238619186083197, 0.238619186083197, 0.661209386466265, 0.932469514203152];
        w=[0.171324492379170, 0.360761573048139, 0.467913934572691, 0.467913934572691, 0.360761573048139, 0.171324492379170];
    end   
    
    % here the weights include the scaling  of the [-1,1] length segment
    % with respect to a unit segment [0,1], so the weights are all basically
    % multiplied by 2...
    % the GP positions are also spread over the space [-1,1] and ned to be
    % adjusted for the [0,1] segment (ksi e [0,1])
    
    
    
