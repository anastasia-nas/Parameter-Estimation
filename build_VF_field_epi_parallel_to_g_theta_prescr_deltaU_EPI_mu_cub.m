function Prescr_Epi_Nodes=build_VF_field_epi_parallel_to_g_theta_prescr_deltaU_EPI_mu_cub(Nodes,Elements_cmgui,Boundaries,Nodes_per_elem_dir,GP_per_elem_dir,endo_patch_ID,epi_patch_ID,base_patch_ID,apex_hole_patch_ID,fig_ind,max_disp)


% created by Anastasia on 11-03-2020
%-------------------------------------------------------------------------------------------------------------------------
% Purpose: 
%-------------------------------------------------------------------------------------------------------------------------
% to use for improving my material parameter estimation routine by removing
% the effect of the unknown boundary tractions.
% Here I will create a VF that is tangential to epi (and thus the pressure 
% at epi will be zero) and  zero at the nodes that are at base and
% hole to ensure delta_Wext_base=0 and delta_Wext_hole=0. To do this I will
% run a cheart simulation of a Stokes problem (satisfies incompressibility)
% and with Dirichlet conditions at epi  which satisfy delta_u_epi: tangential
% to epi and zero at apex and base. The Dirichlet BCs at Epi will be created
% based on finding the g_theta field and using a sin(mu) weight (where mu
% is the apicobasal variable which will scaled to be  within [0,pi]
% I have also used this for solving with a NeoHookean material
%-------------------------------------------------------------------------------------------------------------------------
% Notes on methodology:
%-------------------------------------------------------------------------------------------------------------------------
% modified from
% build_VF_field_epi_parallel_to_g_theta_create_prescr_deltaU_EPI.m
% to find a better way to shift mu towards a more even distribution from
% apex to base. Exploiting the fact that the distribution of mu with
% apicobasal distance resembles x^(1/2) or x^(1/3) I will raise mu from
% the solution of the Laplace problem in the power of 2 or 3.. (for *cub
% version I raise to cubic) - also here I specify the maximum displacement
% at the meridian (where mu =0.5) so I have more control of the field I am
% creating.
%---notes on methodology from previous version
%(build_VF_field_epi_parallel_to_g_theta_create_prescr_deltaU_EPI.m):----
% find local apicobasal and transmural directions by solving the Laplace
% Dirichlet problem with known boundaries the apex-base and endo-epi
% respectively. The local circumferential will be found by the cross
% product of the two.
% then this "g_theta" field (local circumferential) will be scaled by a
% weight that follows the sin(mu) distribution where mu is the normalised
% scalar quantity which is known from solving the Laplace problem but now
% will be scaled to belong to [0,pi] instead of [0,1] as was the case for
% the laplace problem. Then the weight applied to the g_theta field will be 
% sin(mu) which should yield a tangential field to epi with vectors at base
% and hole equal to zero.
%-------------------------------------------------------------------------------------------------------------------------
% Input:
%-------------------------------------------------------------------------------------------------------------------------
% Nodes: Nodes of FE mesh that you need to create the Dirichlet BCs at epi for
% Elements_cmgui: element connectivity of FE mesh
% Boundaries: boundary matrix accoring to cheart
% Nodes_per_elem_dir: Lagrange interpolation order used (4 for cubic Lagrange)
% GP_per_elem_dir: number of GPs per element direction used for quadrature
% endo_patch_ID: ID of endo boundary in "Boundaries" -- this will be a list of patch IDs which will be assigned to endocardial elements (to allow for different BCs at endo)
% epi_patch_ID: ID of epi boundary in "Boundaries"
% base_patch_ID: ID of base boundary in "Boundaries"
% apex_hole_patch_ID: ID of apical hole boundary in "Boundaries" --use =[] if there is none, although all my scripts for solving Laplace work only with "cylindirical shapes" 
%                       (where there are two holes in the topology) and also I think this is beneficial possibly for the Stokes flow problem (where you're essentially imposing a torsion) 
%                       although this could work just by assigning delta_u=0 at the apical nodes as well) -- just rewrite the solve Laplace script in that case.
% fig_ind: leave empty (=[]) if you don't want to plot this or specify figure window if you want to plot this
% flag_scale=1 if you want to scale the prescribed nodes to be ~ 1/20 of mesh length
%-------------------------------------------------------------------------------------------------------------------------
% Output:
%-------------------------------------------------------------------------------------------------------------------------
% Prescr_Epi_Nodes: Nodes size matrix with prescribed velocities at the epicardial nodes and zero everywhere else.
%-------------------------------------------------------------------------------------------------------------------------
% SOS points:
%-------------------------------------------------------------------------------------------------------------------------
% use of mesh with 2 holes so far due to current problem definintions in solve_Laplace_between_oppos_boundaries_v1
if isempty(fig_ind)==0
    flag_plot=1;
else
    flag_plot=0;
end
boundary_1_ID_Phi_1=endo_patch_ID; boundary_2_ID_Phi_0=epi_patch_ID;
%  [Phi_Nodes_CG,Grad_Phi_Nodes_CG,Grad_Phi_dist_CG,Grad_Phi_Nodes_CG_discont,Discon_mesh_Nodes,Discon_mesh_Elements]=solve_Laplace_between_oppos_boundaries_v1(boundary_1_ID,boundary_2_ID,Nodes,Elements_cmgui,Boundaries,Nodes_per_elem_dir,GP_per_elem_dir,flag_plot);
[Lambda_Nodes_CG,Grad_Lambda_Nodes_CG,Grad_Lambda_dist_CG,Grad_Lambda_Nodes_CG_discont,Discon_mesh_Nodes,Discon_mesh_Elements,fig_ind]=solve_Laplace_between_oppos_boundaries_v1(boundary_1_ID_Phi_1,boundary_2_ID_Phi_0,Nodes,Elements_cmgui,Boundaries,Nodes_per_elem_dir,GP_per_elem_dir,flag_plot,fig_ind);

boundary_1_ID_Phi_1=base_patch_ID; boundary_2_ID_Phi_0=apex_hole_patch_ID;
[Mu_Nodes_CG,Grad_Mu_Nodes_CG,Grad_Mu_dist_CG,Grad_Mu_Nodes_CG_discont,Discon_mesh_Nodes,Discon_mesh_Elements,fig_ind]=solve_Laplace_between_oppos_boundaries_v1(boundary_1_ID_Phi_1,boundary_2_ID_Phi_0,Nodes,Elements_cmgui,Boundaries,Nodes_per_elem_dir,GP_per_elem_dir,flag_plot,fig_ind);

Grad_Theta_Nodes_CG=zeros(size(Grad_Lambda_Nodes_CG));
for n_nod=1:size(Grad_Lambda_Nodes_CG,1)
    grad_lambda=Grad_Lambda_Nodes_CG(n_nod,:);
    grad_mu=Grad_Mu_Nodes_CG(n_nod,:);
    Grad_Theta_Nodes_CG(n_nod,:)=cross(grad_lambda,grad_mu)/norm(cross(grad_lambda,grad_mu)); % normalised
end


%% the profile of Phi apicobasally ( due to the div(gradPhi) constraint in the Laplace problem
% that ensures incompressibility of Phi ) changes swiftly around the apex
% and very slowly towards the base. 

% A)) So if you were to plot Phi with apicobasal distance the profile would look like that if you prescribe
% Phi=1 at base and Phi=0 at apex:
%    1|       %%%%%%%%%%%%
%     |     %%%  
%  Mu |   %%
%     |  %
%    0| %
%     --------------------------------> apicobasal distance
%     apex                            base

% B)) whereas if you prescribe Phi=1 at apex and Phi=0 at base, the profile of Mu against the apicobasal distance would look like this:
%   1 |    %
%     |     %
% Mu  |      %%
%     |        %%%
%   0 |           %%%%%%%%%%%%%
%     ---------------------------------> apicobasal distance
%    apex                              base

% Observing case A)) looks like the log function plot so I will turn Mu
% into exp(Mu) and normalise so that Mu_mod still is in [0,1].
% Mu_Nodes_CG_mod_pre=exp(Mu_Nodes_CG+100); % I add 5 to make mu lie in the range [5,6] where the effect of exp will be more pronounced so that indeed the values towards the base will be more differentiated 
% Mu_Nodes_CG_mod_norm=(Mu_Nodes_CG_mod_pre-min(Mu_Nodes_CG_mod_pre))/(max(Mu_Nodes_CG_mod_pre)-min(Mu_Nodes_CG_mod_pre)); % normalised between 0 1
% 
% % vale kai mia riza na meiwseis diafora mi
% % case B looks like 1/exp function
% 
% Mu_Nodes_CG_pi=Mu_Nodes_CG*pi;
% Mu_Nodes_CG_mod_pi=Mu_Nodes_CG_mod_norm*pi;
% for a mu profile that resembles x^(1/3) or x^(1/2) raising mu to the
% cubic might solve the problem (and mu will still be in [0,1]) :
Mu_Nodes_CG_mod=Mu_Nodes_CG.^3;
VF_length_weight=sin(Mu_Nodes_CG_mod*pi); % den einai kalo dioti logw tis katanomis tou phi stous komvous pou dinei phi poly mikro konta sto base kai mikro sto apex (dioti to provlima Laplace eksasfalizei div(gradPhi)=0 ara to gradPhi einai katallilo gia incompressiblity (an ki emena den me noiazei to incompressibility aytou tou field afou to xrisimopoiw mono gia na parw to prescribed disp @Epi.)
%skepsou mia synartisi san to sin(mu) pou na dinei 0 sta akra alla pio
%symmetriki katanomi gia ena profil mu typou: (opou to max mu den vrisketai
%sti mesi tis apostasis alla pros to mu 0 (i to mu pi arkei na valw to pi
%sto apex (anastrepse tous orismous boundary1 kai boundary2 sto provlima
%Laplace Dirichlet.)
               %                
            %%  %
         %%%     %
      %%%        %
   %%%            %
%%%               % 

%%
VF_nodes=([VF_length_weight,VF_length_weight,VF_length_weight].*Grad_Theta_Nodes_CG)*max_disp;
basal_Nodes_IDs=sort(unique(reshape(Boundaries(find(Boundaries(:,Nodes_per_elem_dir^2+2)==base_patch_ID),2:Nodes_per_elem_dir^2+1),[],1)));
apical_hole_Nodes_IDs=sort(unique(reshape(Boundaries(find(Boundaries(:,Nodes_per_elem_dir^2+2)==apex_hole_patch_ID),2:Nodes_per_elem_dir^2+1),[],1)));
VF_nodes(basal_Nodes_IDs,:)=zeros(length(basal_Nodes_IDs),3); % make sure nodes at base and hole have VF=0 and there are no round-offs
VF_nodes(apical_hole_Nodes_IDs,:)=zeros(length(apical_hole_Nodes_IDs),3); % make sure nodes at base and hole have VF=0 and there are no round-offs


epicardial_Nodes_IDs=sort(unique(reshape(Boundaries(find(Boundaries(:,Nodes_per_elem_dir^2+2)==epi_patch_ID),2:Nodes_per_elem_dir^2+1),[],1)));
delta_u_Nodes_epi=VF_nodes(epicardial_Nodes_IDs,:);

Prescr_Epi_Nodes=zeros(size(Nodes)); % Prescribed disp/vel file for cheart: values only at epi nodes and zeros for all other nodes..
Prescr_Epi_Nodes(epicardial_Nodes_IDs,:)=delta_u_Nodes_epi;


% present:
if isempty(fig_ind)==0
    fig_ind=fig_ind+1;
    figure(fig_ind);
    scatter3(Nodes(:,1),Nodes(:,2),Nodes(:,3),20*ones(size(Mu_Nodes_CG)),Mu_Nodes_CG,'f');
    colorbar;
    hold on;
    title('Mu')
    hold off;
    
    fig_ind=fig_ind+1;
    figure(fig_ind);
    scatter3(Nodes(:,1),Nodes(:,2),Nodes(:,3),20*ones(size(Mu_Nodes_CG_mod_norm)),Mu_Nodes_CG_mod_norm,'f');
    colorbar;
    hold on;
    title('Mu field mod and VF field');
    plot_cmgui_ordered_mesh(Elements_cmgui,Nodes,Nodes_per_elem_dir,3)
    hold on;
    quiver3(Nodes(:,1),Nodes(:,2),Nodes(:,3),VF_nodes(:,1),VF_nodes(:,2),VF_nodes(:,3));
    hold off;

    VF_nodes_norm=sqrt(sum(VF_nodes.^2,2));
    fig_ind=fig_ind+1;
    scatter3(Nodes(:,1),Nodes(:,2),Nodes(:,3),20*ones(size(VF_nodes_norm)),VF_nodes_norm,'f');
    colorbar;
    title('VF field norm')
    %einai akoma xamila...
    
    fig_ind=fig_ind+1;
    figure(fig_ind);
    plot_cmgui_ordered_mesh(Elements_cmgui,Nodes,Nodes_per_elem_dir,3);
    hold on;
    quiver3(Nodes(:,1),Nodes(:,2),Nodes(:,3),VF_nodes(:,1),VF_nodes(:,2),VF_nodes(:,3));
    title('all VF nodes');
    hold off;
    
    fig_ind=fig_ind+1;
    figure(fig_ind);
    plot_cmgui_ordered_mesh(Elements_cmgui,Nodes,Nodes_per_elem_dir,3);
    hold on;
    quiver3(Nodes(epicardial_Nodes_IDs,1),Nodes(epicardial_Nodes_IDs,2),Nodes(epicardial_Nodes_IDs,3),delta_u_Nodes_epi(:,1),delta_u_Nodes_epi(:,2),delta_u_Nodes_epi(:,3));
    title(' VF nodes just at epi');
    hold off;
end
 
