function thetaN_thetaksi=Shapefunction_derivative(ksi_location,Nodes_per_Elem_dir)

%-----------------------------------------------------------------------------------
% Created ~ April 2014
% INPUT: Element nodes in matrix: row: internal node order number
%                                 col: i coordinate  of node's position
% OUTPUT:  thetaN_thetaksi: (Nodes_per_Elem_dir^3,3): row: i internal node's i shapefunction derivative
%                                                     col: j: shapefunction derivative in j direction
%---------------------------------------------------------------------------------------------


% vriskei ti paragwgo tou N me to ksi
% thetaN_thetaksi: (Nodes_per_Elem_dir^3,3) -> anagjkastika tha exei komvous sta rows kai ksi directions sta cols 
% ksi_location=[ksi1, ksi2, ksi3];

% internal_node_order : KATA cmgui oxi cheart:

for n_1=1:Nodes_per_Elem_dir % kata ksi1
    for n_2=1:Nodes_per_Elem_dir %ksi2
        for n_3=1:Nodes_per_Elem_dir %ksi3
            internal_node_order((n_3-1)*Nodes_per_Elem_dir^2+(n_2-1)*Nodes_per_Elem_dir+n_1,:)=[n_1,n_2,n_3];
        end
    end
end

for n_nod=1:Nodes_per_Elem_dir
    ksi_along_elem_dir(n_nod)=(n_nod-1)/(Nodes_per_Elem_dir-1);
end

thetaN_thetaksi=ones(Nodes_per_Elem_dir^3,3);
for nod=1:Nodes_per_Elem_dir^3 %total nodes per element
    
    for th_ksi=1:3
        for n_ksi=1:3
            
            num(n_ksi)=1;
            denom(n_ksi)=1;
            % along ksi n_ksi
            for n_dir=1:Nodes_per_Elem_dir % kata ksi1
                if n_dir~=internal_node_order(nod,n_ksi)
                    if n_ksi~=th_ksi
                        num(n_ksi)=num(n_ksi)*(ksi_location(n_ksi)-ksi_along_elem_dir(n_dir));
                    end
                    denom(n_ksi)=denom(n_ksi)*(ksi_along_elem_dir(internal_node_order(nod,n_ksi))-ksi_along_elem_dir(n_dir));
                end
            end
%             if n_ksi~=th_ksi
%                 thetaN_thetaksi(th_ksi,nod)=thetaN_thetaksi(th_ksi,nod)*num/denom;
%             end
            if n_ksi==th_ksi
                num(n_ksi)=0; % tha ksanaoristei 
                % calculate the term that is f'gh+fg'h+fgh'
                for n_dir=1:Nodes_per_Elem_dir
                    factor(n_dir)=1;
                    if n_dir~=internal_node_order(nod,n_ksi)
                        for m_dir=1:Nodes_per_Elem_dir
                            if m_dir~=internal_node_order(nod,n_ksi) && m_dir ~=n_dir
                                factor(n_dir)=factor(n_dir)*(ksi_location(n_ksi)-ksi_along_elem_dir(m_dir));
                            end
                        end
                        num(n_ksi)=num(n_ksi)+factor(n_dir);
                    end
                    
                end
            end
        end
            
        for n_ksi=1:3
            thetaN_thetaksi(nod,th_ksi)=thetaN_thetaksi(nod,th_ksi)*num(n_ksi)/denom(n_ksi);
        end
   
    end
end