function thetaNsurf_thetaksi=Surface_Shapefunction_derivative(ksi_surf_location,Nodes_per_Elem_dir)

% vriskei ti paragwgo tou N me to ksi
% thetaN_thetaksi: (Nodes_per_Elem_dir^2,3) -> anagjkastika tha exei komvous sta rows kai ksi directions sta cols 
% ksi_location=[ksi, ni]; the natural coordinates of the surface

% internal_surf_node_order : KATA cmgui oxi cheart:
for n_1=1:Nodes_per_Elem_dir % kata ksi1
    for n_2=1:Nodes_per_Elem_dir %ksi2
       
        internal_surf_node_order((n_2-1)*Nodes_per_Elem_dir+n_1,:)=[n_1,n_2];
        
    end
end

for n_nod=1:Nodes_per_Elem_dir
    ksi_along_elem_dir(n_nod)=(n_nod-1)/(Nodes_per_Elem_dir-1);
end

thetaNsurf_thetaksi=ones(Nodes_per_Elem_dir^2,2);
for nod=1:Nodes_per_Elem_dir^2 %total nodes per surface element
    
    for th_ksi=1:2
        for n_ksi=1:2
            
            num(n_ksi)=1;
            denom(n_ksi)=1;
            % along ksi n_ksi
            for n_dir=1:Nodes_per_Elem_dir % kata ksi1
                if n_dir~=internal_surf_node_order(nod,n_ksi)
                    if n_ksi~=th_ksi
                        num(n_ksi)=num(n_ksi)*(ksi_surf_location(n_ksi)-ksi_along_elem_dir(n_dir));
                    end
                    denom(n_ksi)=denom(n_ksi)*(ksi_along_elem_dir(internal_surf_node_order(nod,n_ksi))-ksi_along_elem_dir(n_dir));
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
                    if n_dir~=internal_surf_node_order(nod,n_ksi)
                        for m_dir=1:Nodes_per_Elem_dir
                            if m_dir~=internal_surf_node_order(nod,n_ksi) && m_dir ~=n_dir
                                factor(n_dir)=factor(n_dir)*(ksi_surf_location(n_ksi)-ksi_along_elem_dir(m_dir));
                            end
                        end
                        num(n_ksi)=num(n_ksi)+factor(n_dir);
                    end
                    
                end
            end
        end
            
        for n_ksi=1:2
            thetaNsurf_thetaksi(nod,th_ksi)=thetaNsurf_thetaksi(nod,th_ksi)*num(n_ksi)/denom(n_ksi);
        end
   
    end
end