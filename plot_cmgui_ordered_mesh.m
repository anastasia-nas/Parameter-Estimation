function plot_cmgui_ordered_mesh(Elements_ala_cmgui,Nodes,Nodes_per_elem_dir,dimension)

%% created 2 July 2017
% purpose: function for plotting meshes (surface meshes or 3D meshes) in 3
% dimensions. 

%----under construction---!
% Elements : the elements that you want to plot
% Nodes : the full Nodes matrix that corresponds to elements
% dimension: 
% %         =2   ; if it is a surface mesh (Nodes_per_elem_dir^2 nodes per element)
% %         =3   ; if it is a surface mesh (Nodes_per_elem_dir^3 nodes per element)

% % prerequisite: name figure and hold it on: e.g.
% %                                              figure(2000); hold on;

also_print_nodes_IDs=1; % =1 to print also nodes IDs at nodes
only_plot_nodes_sans_connectivity=0; % don't plot element connectivity only scatter nodes.

[indicesEDGE_ksi1ksi2eq0, indicesEDGE_ksi1ksi2eq1, indicesEDGE_ksi1eq1ksi2eq0, indicesEDGE_ksi1eq0ksi2eq1,...
        indicesEDGE_ksi2ksi3eq0, indicesEDGE_ksi2ksi3eq1, indicesEDGE_ksi2eq1ksi3eq0, indicesEDGE_ksi2eq0ksi3eq1,...
        indicesEDGE_ksi3ksi1eq0, indicesEDGE_ksi3ksi1eq1, indicesEDGE_ksi3eq0ksi1eq1, indicesEDGE_ksi3eq1ksi1eq0]=internalNodeIDs_corresponding_natural_coord_EDGES_library(Nodes_per_elem_dir); 
if dimension==3; % an to mesh einai 3D

    if only_plot_nodes_sans_connectivity==1
        scatter3(Nodes(:,1),Nodes(:,2),Nodes(:,3),'.')
        
    else
        for m_el=1:size(Elements_ala_cmgui,1)
            % gia kathe mia apo tis 12 edges tou 3D elem:
            seg1=indicesEDGE_ksi1ksi2eq0;
    %         seg=[seg1,seg1(1)]; % giati to kanw kleisto kyklwma?? --modified
    %         april 18 2018 :vazw seg=seg1; an den doulevei kane undo kai
    %         ksanavale seg=[seg1,seg1(1)];
            seg=seg1;
            line(Nodes(Elements_ala_cmgui(m_el,seg),1),Nodes(Elements_ala_cmgui(m_el,seg),2),Nodes(Elements_ala_cmgui(m_el,seg),3),'Color','g');
            seg1=indicesEDGE_ksi1ksi2eq1;
    %         seg=[seg1,seg1(1)];
            seg=seg1;
            line(Nodes(Elements_ala_cmgui(m_el,seg),1),Nodes(Elements_ala_cmgui(m_el,seg),2),Nodes(Elements_ala_cmgui(m_el,seg),3),'Color','g');
            seg1=indicesEDGE_ksi1eq1ksi2eq0;
    %         seg=[seg1,seg1(1)];
            seg=seg1;
            line(Nodes(Elements_ala_cmgui(m_el,seg),1),Nodes(Elements_ala_cmgui(m_el,seg),2),Nodes(Elements_ala_cmgui(m_el,seg),3),'Color','g');
            seg1=indicesEDGE_ksi1eq0ksi2eq1;
    %         seg=[seg1,seg1(1)];
            seg=seg1;
            line(Nodes(Elements_ala_cmgui(m_el,seg),1),Nodes(Elements_ala_cmgui(m_el,seg),2),Nodes(Elements_ala_cmgui(m_el,seg),3),'Color','g');
            seg1=indicesEDGE_ksi2ksi3eq0;
    %         seg=[seg1,seg1(1)];
            seg=seg1;
            line(Nodes(Elements_ala_cmgui(m_el,seg),1),Nodes(Elements_ala_cmgui(m_el,seg),2),Nodes(Elements_ala_cmgui(m_el,seg),3),'Color','g');
            seg1=indicesEDGE_ksi2ksi3eq1;
    %         seg=[seg1,seg1(1)];
            seg=seg1;
            line(Nodes(Elements_ala_cmgui(m_el,seg),1),Nodes(Elements_ala_cmgui(m_el,seg),2),Nodes(Elements_ala_cmgui(m_el,seg),3),'Color','g');
            seg1=indicesEDGE_ksi2eq1ksi3eq0;
    %         seg=[seg1,seg1(1)];
            seg=seg1;
            line(Nodes(Elements_ala_cmgui(m_el,seg),1),Nodes(Elements_ala_cmgui(m_el,seg),2),Nodes(Elements_ala_cmgui(m_el,seg),3),'Color','g');
            seg1=indicesEDGE_ksi2eq0ksi3eq1;
    %         seg=[seg1,seg1(1)];
            seg=seg1;
            line(Nodes(Elements_ala_cmgui(m_el,seg),1),Nodes(Elements_ala_cmgui(m_el,seg),2),Nodes(Elements_ala_cmgui(m_el,seg),3),'Color','g');
            seg1=indicesEDGE_ksi3ksi1eq0;
    %         seg=[seg1,seg1(1)];
            seg=seg1;
            line(Nodes(Elements_ala_cmgui(m_el,seg),1),Nodes(Elements_ala_cmgui(m_el,seg),2),Nodes(Elements_ala_cmgui(m_el,seg),3),'Color','g');
            seg1=indicesEDGE_ksi3ksi1eq1;
    %         seg=[seg1,seg1(1)];
            seg=seg1;
            line(Nodes(Elements_ala_cmgui(m_el,seg),1),Nodes(Elements_ala_cmgui(m_el,seg),2),Nodes(Elements_ala_cmgui(m_el,seg),3),'Color','g');
            seg1=indicesEDGE_ksi3eq0ksi1eq1;
    %         seg=[seg1,seg1(1)];
            seg=seg1;
            line(Nodes(Elements_ala_cmgui(m_el,seg),1),Nodes(Elements_ala_cmgui(m_el,seg),2),Nodes(Elements_ala_cmgui(m_el,seg),3),'Color','g');
            seg1=indicesEDGE_ksi3eq1ksi1eq0;
    %         seg=[seg1,seg1(1)];
            seg=seg1;
            line(Nodes(Elements_ala_cmgui(m_el,seg),1),Nodes(Elements_ala_cmgui(m_el,seg),2),Nodes(Elements_ala_cmgui(m_el,seg),3),'Color','g');
            text(mean(Nodes(Elements_ala_cmgui(m_el,:),1)),mean(Nodes(Elements_ala_cmgui(m_el,:),2)),mean(Nodes(Elements_ala_cmgui(m_el,:),3)),num2cell(m_el))
            if also_print_nodes_IDs==1
                text(Nodes(Elements_ala_cmgui(m_el,:),1),Nodes(Elements_ala_cmgui(m_el,:),2),Nodes(Elements_ala_cmgui(m_el,:),3),num2cell(Elements_ala_cmgui(m_el,:)))
            end
        end
    end
end

if dimension==2
%     scatter3(Nodes(:,1),Nodes(:,2),Nodes(:,3),'*')

    for m_el=1:size(Elements_ala_cmgui,1)
        seg=1:Nodes_per_elem_dir;
        line(Nodes(Elements_ala_cmgui(m_el,seg),1),Nodes(Elements_ala_cmgui(m_el,seg),2),Nodes(Elements_ala_cmgui(m_el,seg),3),'Color','b');
        seg=(0:(Nodes_per_elem_dir-1))*Nodes_per_elem_dir+1;
        line(Nodes(Elements_ala_cmgui(m_el,seg),1),Nodes(Elements_ala_cmgui(m_el,seg),2),Nodes(Elements_ala_cmgui(m_el,seg),3),'Color','b');
        seg=Nodes_per_elem_dir*(1:Nodes_per_elem_dir);
        line(Nodes(Elements_ala_cmgui(m_el,seg),1),Nodes(Elements_ala_cmgui(m_el,seg),2),Nodes(Elements_ala_cmgui(m_el,seg),3),'Color','b');
        seg=Nodes_per_elem_dir^2-Nodes_per_elem_dir+1:Nodes_per_elem_dir^2;
        line(Nodes(Elements_ala_cmgui(m_el,seg),1),Nodes(Elements_ala_cmgui(m_el,seg),2),Nodes(Elements_ala_cmgui(m_el,seg),3),'Color','b');
        text(mean(Nodes(Elements_ala_cmgui(m_el,:),1)),mean(Nodes(Elements_ala_cmgui(m_el,:),2)),mean(Nodes(Elements_ala_cmgui(m_el,:),3)),num2cell(m_el));
    end
end
