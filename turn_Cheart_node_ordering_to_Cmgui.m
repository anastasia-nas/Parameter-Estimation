function [Reordered_Elements_Cubic,error]=turn_Cheart_node_ordering_to_Cmgui(Elements_Cubic,Nodes_per_elem_dir)


% modified  on October 1st 2019 to account for different interpolation
% orders (quadratic and cubic).. 


% the function reorders cubic lagrange elemnts so that they have the same
% node ordering as in cmgui (no priority to corner nodes )

%Elements_Cubic: ola ta Elements tou mesh sou: rows: element no, columns:
%nodes as ordered in element

% Nodes_per_Elem_dir: nodes per ksi directions in each element: e.g. 4 for
% cubic interpolation..

% if Nodes_per_Elem_dir~=4
%     error=1;
% else
%     error=0;
% end
% Reordered_Elements_Cubic=[Elements_Cubic(:,1),Elements_Cubic(:,9:10),Elements_Cubic(:,2),Elements_Cubic(:,11:18),Elements_Cubic(:,3),Elements_Cubic(:,19:20),Elements_Cubic(:,4),Elements_Cubic(:,21:52),Elements_Cubic(:,5),Elements_Cubic(:,53:54),Elements_Cubic(:,6),Elements_Cubic(:,55:62),Elements_Cubic(:,7),Elements_Cubic(:,63:64),Elements_Cubic(:,8)];

Cmgui_order_old=[1,9:10,2,11:18,3,19:20,4,21:52,5,53:54,6,55:62,7,63:64,8];
%% first build a cmgui order to cheart order map 


Cheart_order=zeros(1,Nodes_per_elem_dir^3); % preallocate

for m_el=1:Nodes_per_elem_dir^3 
    ksi_3=floor((m_el-1)/Nodes_per_elem_dir^2)+1; % ranging from 1 to Nodes_per_elem_dir
    ksi_2=floor(((m_el-(ksi_3-1)*Nodes_per_elem_dir^2)-1)/Nodes_per_elem_dir)+1; % ranging from 1 to Nodes_per_elem_dir
    ksi_1=(m_el-((ksi_3-1)*Nodes_per_elem_dir^2+(ksi_2-1)*Nodes_per_elem_dir)); %  ranging from 1 to Nodes_per_elem_dir
    flag_corner_node=0;
    % if it's a corner node give priority
    if ksi_1==1 || ksi_1==Nodes_per_elem_dir
        if ksi_2==1 || ksi_2==Nodes_per_elem_dir
            if ksi_3==1 || ksi_3==Nodes_per_elem_dir
                % using round(ksi_3/Nodes_per_elem_dir)-- for corner nodes this will be =0 or=1-- to convert the ksi
                % layer ID in the higher interpolation (Nodes_per_elem_dir)
                % order to linear in order to find corner node node
                current_corner_node_ID=round(ksi_3/Nodes_per_elem_dir)*4+round(ksi_2/Nodes_per_elem_dir)*2+round(ksi_1/Nodes_per_elem_dir)+1; %  the last corner node that is filled in the element 
                Cheart_order(m_el)=current_corner_node_ID;
                flag_corner_node=1;
            end 
        end
    end
    
   if flag_corner_node==0 % an to node den einai corner node
       Cheart_order(m_el)=m_el+8-current_corner_node_ID; % to ID tou komvou m_el ala cheart einai oso alla cmgui -(8+x) pou einai ta corner nodes pou mpikan prwta meion osa corner nodes (x) proigountai tou komvou m_el sto cmgui ordering (opote den epetai o ekastote komvos 8 corner nodes alla (8-x) corner nodes)
   end
end

       
 Reordered_Elements_Cubic=Elements_Cubic(:,Cheart_order);
    
    
    
    %% edw eftiaxna to map kateytheian - alla poly complique kai variemai
% %     if m_el<=8 % for the first 8 corner nodes of a hexadral element there is prioritisation in cheart
% %         ksi_3=floor((m_el-1)/4); % this will be 0 (for nodes 1:4) or 1 (for nodes 7:8)
% %         ksi_2=floor(((m_el-ksi_3*4)-1)/2); % again 0 (for nodes 1,2,5,6) or 1 (for nodes 3,4,7,8)
% %         ksi_1=(m_el-(ksi_3*4+ksi_2*2))-1; % 0 (for nodes 1,3,5,7) or 1 for nodes 2,4,6,8)
% %         Cheart_order(m_el)=ksi_3*(Nodes_per_elem_dir^3-Nodes_per_elem_dir^2)+ksi_2*(Nodes_per_elem_dir^2-Nodes_per_elem_dir)+ksi1*Nodes_per_elem_dir;
% %     else 
        
    

    
error=0; % no use for this anymore - just leave it for compatibility with earlier scripts that call this function