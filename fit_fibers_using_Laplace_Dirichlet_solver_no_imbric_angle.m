function [Fibers_nodes_con_or_discon,flag_non_perpendicular_grads_mu_lambda,Fibers_discon_nodes_position,Fibers_discon_elem]=fit_fibers_using_Laplace_Dirichlet_solver_no_imbric_angle(fiber_epi_degrees,fiber_endo_degrees,flag_build_discont_fiber_mesh,Elements_cmgui,Nodes,Boundaries,Nodes_per_elem_dir,GP_per_elem_dir,endo_patch_ID,epi_patch_ID,base_patch_ID,apex_hole_patch_ID,fig_ind)

% created by Anastasia on 10-01-2020
%-------------------------------------------------------------------------------------------------------------------------
% Purpose: 
%-------------------------------------------------------------------------------------------------------------------------
% create my own fiber fitting tool based on solving the Laplace Dirichlet
% problem. For this version fibers are assumed to be tangential to the the
% circumferential-longitudinal (azimuthial) surfaces in the heart (no
% imbrication angle). General methodology follows Bayer 2012.
%-------------------------------------------------------------------------------------------------------------------------
% Notes on methodology:
%-------------------------------------------------------------------------------------------------------------------------
% find local apicobasal and transmural directions by solving the Laplace
% Dirichlet problem with known boundaries the apex-base and endo-epi
% respectively. The local circumferential will be found by the cross
% product of the two.
% Sheet direction is kept constant pointing from endo to epi. This leads to
% unphysiologically oriented cleavage planes (compare to cleavage plane
% angle in LeGrice &Smaill 1995 varying from -90 to ~90 epi to endo).
% However as I only use a transversely isotropic material description this
% is irrelevant in my case
%-------------------------------------------------------------------------------------------------------------------------
% Input:
%-------------------------------------------------------------------------------------------------------------------------
% fiber_epi_degrees: fiber angle at epi in degrees (fiber at epi is usually negative)
% fiber_endo: fiber angle at endo in degrees
% Nodes: Nodes of FE mesh that you need to create the Dirichlet BCs at epi for
% Elements_cmgui: element connectivity of FE mesh
% Boundaries: boundary matrix accoring to cheart
% Neighbours: neighbours matrix ala Cheart (assuming the neighbouring element order is [-ksi3,-ksi2,-ksi1,ksi1,ksi2,ksi3] )
% Nodes_per_elem_dir: Lagrange interpolation order used (4 for cubic Lagrange)
% GP_per_elem_dir: number of GPs per element direction used for quadrature
% endo_patch_ID: ID of endo boundary in "Boundaries"
% epi_patch_ID: ID of epi boundary in "Boundaries"
% base_patch_ID: ID of base boundary in "Boundaries"
% apex_hole_patch_ID: ID of apical hole boundary in "Boundaries" --use =[] if there is none, although all my scripts for solving Laplace work only with "cylindirical shapes" 
%                       (where there are two holes in the topology) and also I think this is beneficial possibly for the Stokes flow problem (where you're essentially imposing a torsion) 
%                       although this could work just by assigning delta_u=0 at the apical nodes as well) -- just rewrite the solve Laplace script in that case.
% fig_ind: leave empty (=[]) if you don't want to plot this or specify figure window if you want to plot this

% flag_build_discont_fiber_mesh: =1 if you create discontinuous fibers - do it for collapsed element meshes or if you get very high differences in normal gradient directions from the Laplace solution 
%                                           for each neighbouring element at the same node.
%-------------------------------------------------------------------------------------------------------------------------
% Output:
%-------------------------------------------------------------------------------------------------------------------------
% Fibers_nodes_con_or_discon: the Nodes matrix size: nodes_num x 9 (this will either correspond to the continuous input mesh (Nodes) or will 
%                             correspond to the discontinuous mesh if flag_build_discont_fiber_mesh=1 
% flag_non_perpendicular_grads_mu_lambda: =1 if the g_lambda unit vectors found from Laplace (with endo - epi fixed boundaries) is 
%                             not perpendicular to the g_mu from the Laplace (with apex hole and base as fixed boundaries)
% Fibers_discon_nodes_position: the discontinuous mesh nodes positions in space (size: Discon_mesh_Nodes x 3). If discontinuous mesh 
%                             not selected will be left empty.
% Fibers_discon_elem: discontinuous mesh element connectivity (Node Ids correspond to rows in Fibers_nodes_con_or_discon matrix 
%                             (if mesh discontinuous) and Fibers_discon_nodes_position). If flag_build_discont_fiber_mesh=0 and
%                             then Fibers_discon_elem will be empty (=[])
%  
%-------------------------------------------------------------------------------------------------------------------------
% SOS points:
%-------------------------------------------------------------------------------------------------------------------------
% NO imbrication angle!
% Neighbour matrix directions:[-ksi3,-ksi2,-ksi1,ksi1,ksi2,ksi3] 
% Sheet direction is kept constant pointing from endo to epi, which leads to unphysiologically oriented cleavage planes (see Methodology above),
%            but OK for transversely isotropic material analysis

%% START CODE:

fiber_endo_pi=fiber_endo_degrees/180*pi;
fiber_epi_pi=fiber_epi_degrees/180*pi;
if isempty(fig_ind)==1
    flag_plot=0;  % for plotting Laplace problem outputs in terms of scalar field and its gradients
else
    flag_plot=1;
end

boundary_1_ID_Phi_1=epi_patch_ID; boundary_2_ID_Phi_0=endo_patch_ID; % this field has direction endo to epi
[Lambda_Nodes_CG,Grad_Lambda_Nodes_CG,Grad_Lambda_dist_CG,Grad_Lambda_Nodes_CG_discont,Discon_mesh_Nodes1,Discon_mesh_Elements1,fig_ind]=solve_Laplace_between_oppos_boundaries_v1(boundary_1_ID_Phi_1,boundary_2_ID_Phi_0,Nodes,Elements_cmgui,Boundaries,Nodes_per_elem_dir,GP_per_elem_dir,flag_plot,fig_ind);

boundary_1_ID_Phi_1=base_patch_ID; boundary_2_ID_Phi_0=apex_hole_patch_ID; % field has zero at apex and 1 at base
[Mu_Nodes_CG,Grad_Mu_Nodes_CG,Grad_Mu_dist_CG,Grad_Mu_Nodes_CG_discont,Discon_mesh_Nodes2,Discon_mesh_Elements2,fig_ind]=solve_Laplace_between_oppos_boundaries_v1(boundary_1_ID_Phi_1,boundary_2_ID_Phi_0,Nodes,Elements_cmgui,Boundaries,Nodes_per_elem_dir,GP_per_elem_dir,flag_plot,fig_ind);

% check discontinuous meshes coincide

if max(sum(abs(round((Discon_mesh_Nodes1-Discon_mesh_Nodes2)*10^8)),2))>0 || max(sum(abs(Discon_mesh_Elements1-Discon_mesh_Elements2),2))>0
    disp('error discontinuous memshes in Laplace do not coincide -- correct this and do not use script for generating discontinuous fibers')
    flag_build_discont_fiber_mesh=0; % should be by default otherwise won't work
    return
else
    Discon_mesh_Nodes=Discon_mesh_Nodes1;
    Discon_mesh_Elements=Discon_mesh_Elements1;
end
flag_non_perpendicular_grads_mu_lambda=0;

% if you choose to create a discontinuous fiber mesh
if flag_build_discont_fiber_mesh==1
    Fibers_discon_nodes_position=Discon_mesh_Nodes;
    Fibers_discon_elem=Discon_mesh_Elements;
    Fiber_mesh_Nodes_pos_TU=Discon_mesh_Nodes; 
    Fiber_mesh_Elements_TU=Discon_mesh_Elements;
    Grad_Mu_Nodes_TU=Grad_Mu_Nodes_CG_discont;
    Grad_Lambda_Nodes_TU=Grad_Lambda_Nodes_CG_discont;
    
    Grad_Theta_Nodes_CG_discont=zeros(size(Grad_Lambda_Nodes_CG_discont)); % preallocate
    for n_nod=1:size(Grad_Lambda_Nodes_CG_discont,1)
        grad_lambda=Grad_Lambda_Nodes_CG_discont(n_nod,:);
        grad_mu=Grad_Mu_Nodes_CG_discont(n_nod,:);
        % check grad mu and grad lambda are perpendicular within numerical accuracy:
        if round(10^10*dot(grad_mu,grad_lambda))~=0
            flag_non_perpendicular_grads_mu_lambda=1;
            % I remove the component parallet to grad_mu from the grad_lambda field to make sure the resulting (grad_lambda,grad_mu,grad_theta) coordinate basis is orthonormal            
            %  -- assuming that grad_lambda is maybe more prone to mistakes as it is a more constrained problem (although it actually may be the other way round)
            Grad_Lambda_Nodes_TU(n_nod,:)=(grad_lambda-dot(grad_lambda,grad_mu)*grad_mu)/norm(grad_lambda-dot(grad_lambda,grad_mu)*grad_mu); 
        end
        Grad_Theta_Nodes_CG_discont(n_nod,:)=cross(grad_mu,grad_lambda)/norm(cross(grad_mu,grad_lambda)); % normalised --anticlockwise direction in top view so that theta-mu-lambda form a right hand system
    end    
    Grad_Theta_Nodes_TU=Grad_Theta_Nodes_CG_discont;
    disp('calculate Lambda_Nodes_TU for discontin mesh!!!'); % Find repeated nodes and then use Lambda_Nodes_CG to build the discontinuous mesh version of it
    
% if you choose to create a continuous fiber mesh
else  
    Fibers_discon_nodes_position=[];
    Fibers_discon_elem=[];
    Fiber_mesh_Nodes_pos_TU=Nodes; 
    Fiber_mesh_Elements_TU=Elements_cmgui;
    Grad_Mu_Nodes_TU=Grad_Mu_Nodes_CG;
    Grad_Lambda_Nodes_TU=Grad_Lambda_Nodes_CG;
    Lambda_Nodes_TU=Lambda_Nodes_CG;
    Grad_Theta_Nodes_CG=zeros(size(Grad_Lambda_Nodes_CG)); % preallocate
    for n_nod=1:size(Grad_Lambda_Nodes_CG,1)
        grad_lambda=Grad_Lambda_Nodes_CG(n_nod,:);
        grad_mu=Grad_Mu_Nodes_CG(n_nod,:);
        % check grad mu and grad lambda are perpendicular within numerical accuracy:
        if round(10^10*dot(grad_mu,grad_lambda))~=0
            flag_non_perpendicular_grads_mu_lambda=1;
            % I remove the component parallet to grad_mu from the grad_lambda field to make sure the resulting (grad_lambda,grad_mu,grad_theta) coordinate basis is orthonormal            
            %  -- assuming that grad_lambda is maybe more prone to mistakes as it is a more constrained problem (although it actually may be the other way round)
            Grad_Lambda_Nodes_TU(n_nod,:)=(grad_lambda-dot(grad_lambda,grad_mu)*grad_mu)/norm(grad_lambda-dot(grad_lambda,grad_mu)*grad_mu); 
        end
        Grad_Theta_Nodes_CG(n_nod,:)=cross(grad_mu,grad_lambda)/norm(cross(grad_mu,grad_lambda)); % normalised --anticlockwise direction in top view so that theta-mu-lambda form a right hand system
    end
   
    Grad_Theta_Nodes_TU=Grad_Theta_Nodes_CG;
end



% create local fiber, sheet and sheet normal directions. I assume f lies on circumferential-apicobasal plane (no imbrication angle) and that sheet
% direction (s) always points from endo to epi (see Dave's thesis p. 186, varying linearly endo to epi.) --also checked with LeGrice 1995 (e.g. Fig 12) 
% where it shows that cleavage planes run radially so that s is always transmural (since f is always tangential) 
Fibers_nodes_con_or_discon=zeros(size(Grad_Mu_Nodes_TU,1),9);
for n_nod=1:size(Grad_Mu_Nodes_TU,1)
%     lambda: you need to find to weigh this so that it varies linearly from 0 at endo to 1 at epi.
%     as the Laplace problem for bounded endo-epi gives an approximate
%     linearly varying field with transmural direction, I'll just use
%     Lambda for the weight.
    
    lambda=Lambda_Nodes_TU(n_nod,:); %lambda=1 @epi and =0 @endo:
    fiber_angle=(1-lambda)*fiber_endo_pi+lambda*fiber_epi_pi;
    
    % so to sum up now we have grad_theta pointing anticlockwise as you
    % look at base from atria, grad_mu pointing apex to base and
    % grad_lambda endo to epi, they're unit vectors and normal to each
    % other so the vectors (grad_theta,grad_mu,grad_lambda) form a right hand
    % system.
    
    f_vec=cos(fiber_angle)*Grad_Theta_Nodes_TU(n_nod,:)+sin(fiber_angle)*Grad_Mu_Nodes_TU(n_nod,:);
    s_vec=Grad_Lambda_Nodes_TU(n_nod,:);
    n_vec=cross(f_vec,s_vec)/norm(cross(f_vec,s_vec));
    Fibers_nodes_con_or_discon(n_nod,:)=[f_vec,s_vec,n_vec];                
end

if isempty(fig_ind)==0
    fig_ind=fig_ind+1;   
    figure(fig_ind);
    scatter3(Nodes(:,1),Nodes(:,2),Nodes(:,3),20*ones(size(Mu_Nodes_CG)),Mu_Nodes_CG,'f');
    colorbar;
    hold on;
    title('Mu')
    hold off;
    
    Basal_Nodes=sort(unique(reshape(Boundaries(find(Boundaries(:,Nodes_per_elem_dir^2+2)==base_patch_ID),2:Nodes_per_elem_dir^2+1),[],1))); % col vector
    Apical_Hole_Nodes=sort(unique(reshape(Boundaries(find(Boundaries(:,Nodes_per_elem_dir^2+2)==apex_hole_patch_ID),2:Nodes_per_elem_dir^2+1),[],1))); % col vector
    
    vector_scale=norm(Nodes(Basal_Nodes(1),:)-Nodes(Apical_Hole_Nodes(1),:))/20; % you want vector length to be approximately 1/20 of the long axis
    
    fig_ind=fig_ind+1;
    figure(fig_ind);
    plot_cmgui_ordered_mesh(Elements_cmgui,Nodes,Nodes_per_elem_dir,3);
    hold on;
    quiver3(Fiber_mesh_Nodes_pos_TU(:,1),Fiber_mesh_Nodes_pos_TU(:,2),Fiber_mesh_Nodes_pos_TU(:,3),Fibers_nodes_con_or_discon(:,1)*vector_scale,Fibers_nodes_con_or_discon(:,2)*vector_scale,Fibers_nodes_con_or_discon(:,3)*vector_scale);
    title('fiber vectors')
    hold off;
    
    fig_ind=fig_ind+1;
    figure(fig_ind);
    plot_cmgui_ordered_mesh(Elements_cmgui,Nodes,Nodes_per_elem_dir,3);
    hold on;
    quiver3(Fiber_mesh_Nodes_pos_TU(:,1),Fiber_mesh_Nodes_pos_TU(:,2),Fiber_mesh_Nodes_pos_TU(:,3),Fibers_nodes_con_or_discon(:,4)*vector_scale,Fibers_nodes_con_or_discon(:,5)*vector_scale,Fibers_nodes_con_or_discon(:,6)*vector_scale);
    title('sheet vectors')
    hold off;

    fig_ind=fig_ind+1;
    figure(fig_ind);
    plot_cmgui_ordered_mesh(Elements_cmgui,Nodes,Nodes_per_elem_dir,3);
    hold on;
    quiver3(Fiber_mesh_Nodes_pos_TU(:,1),Fiber_mesh_Nodes_pos_TU(:,2),Fiber_mesh_Nodes_pos_TU(:,3),Fibers_nodes_con_or_discon(:,7)*vector_scale,Fibers_nodes_con_or_discon(:,8)*vector_scale,Fibers_nodes_con_or_discon(:,9)*vector_scale);
    title('sheet normal vectors')
    hold off;
end





