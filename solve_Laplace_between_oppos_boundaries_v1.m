function [Phi_Nodes_CG,Grad_Phi_Nodes_CG,Grad_Phi_dist_CG,varargout]=solve_Laplace_between_oppos_boundaries_v1(boundary_1_ID_Phi1,boundary_2_ID_Phi0,Nodes,Elements_cmgui,Boundaries,Nodes_per_elem_dir,GP_per_elem_dir,flag_plot,fig_ind_in)

% created by Anastasia on 4-11-2019
%-------------------------------------------------------------------------------------------------------------------------
% Usage: 
%-------------------------------------------------------------------------------------------------------------------------
% for continuous fiber field:
% [Phi_Nodes_CG,Grad_Phi_Nodes_CG,Grad_Phi_dist_CG]=solve_Laplace_between_oppos_boundaries_v1(boundary_1_ID,boundary_2_ID,Nodes,Elements_cmgui,Boundaries,Nodes_per_elem_dir,GP_per_elem_dir,flag_plot,fig_ind_in);
% for discontinuous fiber field:
%  [Phi_Nodes_CG,Grad_Phi_Nodes_CG,Grad_Phi_dist_CG,Grad_Phi_Nodes_CG_discont,Discon_mesh_Nodes,Discon_mesh_Elements,fig_ind_out]=solve_Laplace_between_oppos_boundaries_v1(boundary_1_ID,boundary_2_ID,Nodes,Elements_cmgui,Boundaries,Nodes_per_elem_dir,GP_per_elem_dir,flag_plot,fig_ind_in);
%-------------------------------------------------------------------------------------------------------------------------
% Purpose: 
%-------------------------------------------------------------------------------------------------------------------------
% A solver for the laplace equation where BCs at two opposite boundaries (boundary_1, boundary_2) are prescribed. 
% This is an intermediate step for my fiber fitting tool, where I solve the Laplace problem twice, once setting
% BCs at endo and epi to find the transmural direction and once at apex and
% base to find apicobasal
%-------------------------------------------------------------------------------------------------------------------------
% Notes on methodology:
%-------------------------------------------------------------------------------------------------------------------------
%  I follow the methodology in Bayer 2005 and Bayer 2012. The Laplace equation
%  thetaphi/thetaxi/thetaxi=0 (with BCs phi=1 at boundary_1 and phi=0 at boundary_2)
% is turned into a weak formulation, from where a symmetric "stiffness" matrix 
% can emerge (eq. 7 in Bayer 2005) and taking in mind that the surface integral at the boundary patches where
% no BCs are prescribed is 0 because the gradient of the scalar quantity phi runs tangential to the surfaces 
% and is thus perpendicular to the surface normal.
% Works with Lagrange meshes so far.
%-------------------------------------------------------------------------------------------------------------------------
% Input:
%-------------------------------------------------------------------------------------------------------------------------
% a FE mesh and boundary ID's definitions.
% Nodes: standard nodes matrix
% Elements_cmgui: element connectivity matrix with no corner node priority
% Boundaries: boundary file ala cheart
% Neighbours: neighbour file ala cheart 
% Nodes_per_elem_dir: Lagrange interpolation order
% flag_plot=1 if you want to plot gradients and Phi at nodes, =0 if you want no plots, =2 if you want to plot outputs for both congugate gradient and inv(K) --generally no need to print results for inv(K) since it's rarely correct
% boundary_1_ID= patch ID where Phi = 1 (this can be a list of patch IDs which will be assigned to endocardial elements for example (to allow for different BCs at endo)
% boundary_2_ID= patch ID where Phi = 0
% fig_ind_in= index of previous figure window that is currently open
%-------------------------------------------------------------------------------------------------------------------------
% Output:
%-------------------------------------------------------------------------------------------------------------------------
% Phi_Nodes_CG: the scalar quantity values at the nodes. Obtained from solving for K with the pcg (preconditioned conjugate gradient) matlab function with TOL +10^-7 according to Bayer 2012. The gradient of the quantity denotes the direction between the boundaries
% Grad_Phi_Nodes_CG: gradient vector at each node, size: [total_nodes,3]. The gradient is obtained at each node by averaging over the calculated gradients for the same node at all elements containing the node. (Continuous mesh)
% Grad_Phi_dist_CG: maximum distance between gradient vectors corresponding to the same node calculated at the different elements containing the node size: [total_nodes,3]
% Grad_Phi_Nodes_CG_discont: grad phi at nodes for discontinuous mesh, size: [Nodes_per_elem_dir^3*total_elements,3]
% Discon_mesh_Nodes: coordinates of the discontinuous mesh nodes (size: [Nodes_per_elem_dir^3*total_elements,3])
% Discon_mesh_Elements: element connectivity of the discontinuous mesh (size: [total_elements, Nodes_per_elem_dir^3)
% fig_ind_out: the index of the last figure that was output in this script -- so you can easily use it as input for the next call of the script
%-------------------------------------------------------------------------------------------------------------------------
% SOS points:
%-------------------------------------------------------------------------------------------------------------------------
% BCs are prescribed by patch IDs so far. For finding transmural direction
% that is fine, but for apicobasal will be problem if meshes don't have a
% hole (where you can use as input patch IDs at base and apex hole). In
% that case you need to allow the function to have as input Node IDs that
% will serve as Dirchlet BCs for the Laplace problem.

% Solving Phi_U=(KFF)^-1*(K_FB*Phi_B)::
% inverting worked for cube-linear Lagrange and case JUL_new with quadratic Lagrange
% for half sphere mesh with cubic Lagrange K_FF was near singular (too
% sparse) and result was cccrazzzy. 
% I currently use the matlab pcg function for preconditioned conjugate
% gradient with TOL=10^(-7) as recommended by Bayer 2012.
% older script versions concerns- now solved or no longer pertinent
% (introduced in previous versions and still adopted
% % % I am having issues with the solution. It seems to go from 0 to negative
% % % (~-0.8) and then abruptly to 1 at the other boundary. This might mean the
% % % surface normal is somewhere violated. I have made changes in the shift
% % % from ksi_2D coordinates (at surface GPs) to corresponding ksi_3D location
% % % because in the previous script i used () I was assigning the remaining
% % % ksi_3D locations to ksi_surf(1) and ksi_surf(2) without proper checks.. I
% % % created function "turn_ksi_2D_on_element_boundary_into_ksi_3D.m" for that
% % % purpose..
% % % I also replaced the use of the function
% % % "find_ksi_location_from_XYZ_coords_with_Newton_solve_in_3D_elem.m" at the
% % % part where I estimate the gradient, with new function 
% % % "find_ksi_location_of_local_node_ID_in_3D_element_cmgui_ordering"in order
% % % to make it faster. 



%-------------------------------------------------------------------------------------------------------------------------
% NOTES on version COMPARED TO PREVIOUS SCRIPTS:
%-------------------------------------------------------------------------------------------------------------------------
% Dec 9th 2019 - modified from
% solve_Laplace_between_oppos_boundaries_new_vec_no_parfor_pre.m where I
% am no longer using parfor and this actually helps speed - as parfor was
% slowing this version of the script down. No major vectorisation was done
% either.
% modified from solve_Laplace_between_oppos_boundaries_new.m aiming to
% speed up the code a bit by further vectorisation. I have also erased all
% the commented bits where I used to estimate the surface integrals (and
% which are no longer needed as N_A=0 there). 
% For speeding up what I need to do is store all GPs values required in the
% volume integral for each element and then assign them to nodes. That way
% you don't make the estimations several times for the same element (as all
% the code is "element-o-centric" due to the piecewise interpolation used).
% In the current (as was in solve_Laplace_between_oppos_boundaries_new.m
% version) version of things I estimate values per nodes A and B and
% summing over the elements that include these- but that includes heavy
% operations that include the whole element and numerical integrations and
% are only assigned to one entry (A,B) in K_FF when I should be sparing my
% self the time and go over each element finish every computation and then
% assign what is needed to each (A,B) entry.

% modified from solve_Laplace_between_oppos_boundaries.m by removing the
% surface integral calculations at N_A belonging to free nodes (which holds 
% for both K_FF and K_FB,  as N_A=0 there) I also removed a bug (I had
% forgotten the minus sign when estimating PHi_U=-inv(K_FF)*K_FB*Phi_B.

%% preallocations:
Phi_Nodes_inv=zeros(size(Nodes,1),1);
Phi_Nodes_CG=zeros(size(Nodes,1),1); % CG: short for conjugate gradient
Grad_Phi_Nodes_inv=zeros(size(Nodes));
Grad_Phi_Nodes_CG=zeros(size(Nodes)); 
Grad_Phi_dist_inv=zeros(size(Nodes,1),1);  % max distance between calculated gradients from neighbouring elements per node, Phi calculated from stiffness matrix inversion-> not to be used
Grad_Phi_dist_CG=zeros(size(Nodes,1),1);  % max distance between calculated gradients from neighbouring elements per node, Phi calculated from conjugate gradient

if nargout>=6 % this will only occur if you request for the Discontinuous mesh version of Grad_Phi at the nodes (of a new discontinuous mesh that you need to define with the same element IDs as the continuous one.
    flag_discon_elements=1;
    %also preallocate
    Grad_Phi_Nodes_CG_discont=zeros(size(Elements_cmgui,1)*Nodes_per_elem_dir^3,3);
    Discon_mesh_Nodes=zeros(size(Elements_cmgui,1)*Nodes_per_elem_dir^3,3);
    Discon_mesh_Elements=zeros(size(Elements_cmgui,1),Nodes_per_elem_dir^3);
    if nargout==6  % if no fig_ind_in was specialised
        fig_ind_in=[];
    end
elseif nargout==3 % do nothing
    flag_discon_elements=0;   
else % if nargout is not 0 or 3 then something is wrong: throw error
    disp('error in solve_Laplace_between_oppos_boundaries.. in the number of arguments you have requested as additional output variables should be either 0 or 3');
    disp('outputs should be either: [Phi_Nodes_CG,Grad_Phi_Nodes_CG,Grad_Phi_dist_CG,Grad_Phi_Nodes_CG_discont,Discon_mesh_Nodes,Discon_mesh_Elements] or : [Phi_Nodes_CG,Grad_Phi_Nodes_CG,Grad_Phi_dist_CG]');
    return;
end
    
%% prep: find nodes on boundaries:
% boundary_1:
Bound_1_Elements=Boundaries(ismember(Boundaries(:,Nodes_per_elem_dir^2+2),boundary_1_ID_Phi1),1); % changed "find" into "ismember" to allow for lists of patch ID being boundary 1 or 2 (e.g. when you use two patch IDs to define endo (due to different BCs I want to prescribe to different groups of endo elements)
Bound_2_Elements=Boundaries(ismember(Boundaries(:,Nodes_per_elem_dir^2+2),boundary_2_ID_Phi0),1);
Bound_1_Nodes=sort(unique(reshape(Boundaries(ismember(Boundaries(:,Nodes_per_elem_dir^2+2),boundary_1_ID_Phi1),2:Nodes_per_elem_dir^2+1),[],1)));
Bound_2_Nodes=sort(unique(reshape(Boundaries(ismember(Boundaries(:,Nodes_per_elem_dir^2+2),boundary_2_ID_Phi0),2:Nodes_per_elem_dir^2+1),[],1)));

%% step 1: reorder nodes: nodes that don't lie on boundaries 1,2 go first and nodes on boundaries follow (gia na exw etoimi ti statiki sympiknwsi)
Nodes_ind=(1:size(Nodes,1)).'; % col vector
Free_nodes=Nodes_ind(~ismember(Nodes_ind,[Bound_1_Nodes;Bound_2_Nodes]));
Bound_nodes=[Bound_1_Nodes;Bound_2_Nodes];
Reordered_nodes=[Free_nodes;Bound_nodes]; % first the nodes that don't lie on boundaries, then the rest

%% step 2: build PHI_B (phi vector at the boundary nodes
Phi_Bound1=ones(size(Bound_1_Nodes));
Phi_Bound2=zeros(size(Bound_2_Nodes)); %checking if using equal and opposite values at boundaries solves the weird gradient problem (gradient changed direction very close to the boundary where \Phi is 1 (using \Phi @bound2=-1 instead of =0 - no improvement)
Phi_B=[Phi_Bound1;Phi_Bound2];

%% step 3: construct the "stiffness matrix" in the nodes (specifically the K_FF and K_FB, which is needed for the final residual vector R_fin=R_F-K_FB*PHI_B)
% -- Nov 5 2019: this step currently takes an hour to build 1/3 for a cubic Lagrange interpolation-- need to speed up.
% 15 Nov 2019: timing this for i_case=JUL_new and quadratic Lagrange interpolation it took me 20 mins just to build the two matrices.......!!!!!!!!!!!!!!!!!!!!!!!

%F:free nodes (nodes where no BCs are prescribed), B: bounded nodes (nodes which belong on surfaces for which boundary conditions are prescribed)


% preallocations
K_mat_FF=zeros(length(Free_nodes),length(Free_nodes)); % preallocate %F:free nodes (nodes where no BCs are prescribed), B: bounded nodes (nodes which belong on surfaces for which boundary conditions are prescribed)
K_mat_FB=zeros(length(Free_nodes),length(Bound_nodes)); % preallocate %F:free nodes (nodes where no BCs are prescribed), B: bounded nodes (nodes which belong on surfaces for which boundary conditions are prescribed)
tic

% notation : Node A corresponds to row, Node B to column
%% Try to vectorise these - create vectors that contain all the required info
% 3D elem:
[GPcoords, GPweights]=Gauss_Points_In_3D_Element(GP_per_elem_dir,GP_per_elem_dir,GP_per_elem_dir);
GPweights_tot=GPweights(:,1).*GPweights(:,2).*GPweights(:,3);
% 2D elem:
[surfGPcoords, surfGPweights]=Gauss_Points_In_2D_Element(GP_per_elem_dir,GP_per_elem_dir);
surfGPweights_tot=surfGPweights(:,1).*surfGPweights(:,2);
% build a shapefunction and shapefunction derivative vector order per GP
% for quick vectorisation

for n_GP=1:GP_per_elem_dir^3
    ksi_location_3D=GPcoords(n_GP,:);
    thetaN_thetaksi_per_GP(:,:,n_GP)=Shapefunction_derivative(ksi_location_3D,Nodes_per_elem_dir); % this will need to be a 3D st
end

%% new BUILDING STIFFNESS MATRIX:
for m_el=1:size(Elements_cmgui,1) % ayto ithela na to kanw parfor - exw thema me to parfor - dexetai mono vectors me indices to m_el
    Elem_Nodes=Nodes(Elements_cmgui(m_el,:),:);
    % thelw na ftiakxw gia kathe internal node a,b:
    % \sum_GP dot(\gradN(a,:),\gradN(b,:)) dV_GP w_GP
    element_K=zeros(Nodes_per_elem_dir^3,Nodes_per_elem_dir^3); % initialisation of summation over GPs
    for n_GP=1:GP_per_elem_dir^3
%         ksi_location_3D=GPcoords(n_GP,:);
%         thetaN_thetaksi=Shapefunction_derivative(ksi_location_3D,Nodes_per_elem_dir);
        thetaN_thetaksi=squeeze(thetaN_thetaksi_per_GP(:,:,n_GP));
        grad_N_mat=thetaN_thetaksi*inv(Elem_Nodes.'*thetaN_thetaksi); % size grad_N:[Nodes_per_elem_dir^3,3]
        grad_N_grad_N_trans_mat=grad_N_mat*grad_N_mat.';
        dV=det(Elem_Nodes.'*thetaN_thetaksi);
        element_K=element_K+dV*GPweights_tot(n_GP,:)*grad_N_grad_N_trans_mat;
    end
    % find free nodes indices of element nodes
    local_FREE_nodes_ind=[]; global_FREE_nodes_ind=[];
    local_FREE_nodes_ind=find(ismember(Elements_cmgui(m_el,:),Free_nodes)==1); % the local nodes indices in the element that belong to free nodes
    for n_nod=1:length(local_FREE_nodes_ind)
        global_FREE_nodes_ind(n_nod)=find(Free_nodes==Elements_cmgui(m_el,local_FREE_nodes_ind(n_nod)));
    end
    % and the bounded nodes indices
    local_BOUND_nodes_ind=[]; global_BOUND_nodes_ind=[];
    local_BOUND_nodes_ind=find(ismember(Elements_cmgui(m_el,:),Bound_nodes)==1); % the local nodes indices in the element that belong to free nodes
    for n_nod=1:length(local_BOUND_nodes_ind)
        global_BOUND_nodes_ind(n_nod)=find(Bound_nodes==Elements_cmgui(m_el,local_BOUND_nodes_ind(n_nod)));
    end
    % then split element_K to element_K_FF and element_K_FB
    K_mat_FF(global_FREE_nodes_ind,global_FREE_nodes_ind)=K_mat_FF(global_FREE_nodes_ind,global_FREE_nodes_ind)+element_K(local_FREE_nodes_ind,local_FREE_nodes_ind);
    K_mat_FB(global_FREE_nodes_ind,global_BOUND_nodes_ind)=K_mat_FB(global_FREE_nodes_ind,global_BOUND_nodes_ind)+element_K(local_FREE_nodes_ind,local_BOUND_nodes_ind);
    
    % store directly the element_K_FF to K_FF and element_K_FB to K_FB 
    % now there is another way of going about this and that is to first
    % build the whole K matrix and then re-arrange it to K_FF and K_FB but 
    
end

%% old BUILDING STIFFNESS MATRIX:
% % parfor n_A=1:length(Free_nodes) % I enclose the generation of K_FF and K_FB under a big parfor because I think that the initialisation of parfor maybe costly and it's best to avoid opening and closing
% %     % call again function Gauss_Points_In_3D_Element for each parfor
% %     % iteration to avoid excessive communication overhead
% %     
% %     %% K_FF
% %     v_mat_FF=zeros(1,length(Free_nodes));% preallocate pre-forloop
% %     for n_B=n_A:length(Free_nodes) % for n_B=1:length(Free_nodes):: if you include surface integral estimation then you must loop over all nodes because K is no longer symmetric (at the surface integral term)  % palia to for loop mou pigaine apo for n_A:length(Free_nodes) no need to calcuate the symmetric entries twice -- omws eixa ksexasei to component logw tou \int_A dot(grad\Phi,surf_norm)*N_A dA pou dinei synistwsa kai se n_B nodes pou den anikoun sto surface logw tou grad\Phi.-> ayto mou kostizei 12 lepta.
% %         elements_contain_nodes=[];
% %         if Free_nodes(n_B)~=Free_nodes(n_A)
% %             elements_contain_nodes=find(sum(ismember(Elements_cmgui,[Free_nodes(n_A),Free_nodes(n_B)]),2)==2); %vres ta elements pou periexoun kai tous dyo komvous
% %         else
% %             elements_contain_nodes=find(sum(ismember(Elements_cmgui,Free_nodes(n_A)),2)==1);
% %         end
% %         for m_el=1:length(elements_contain_nodes)
% %             local_ind_nA_in3Delem=find(Elements_cmgui(elements_contain_nodes(m_el),:)==Free_nodes(n_A)); % the local node index of n_A inside the element 
% %             local_ind_nB_in3Delem=find(Elements_cmgui(elements_contain_nodes(m_el),:)==Free_nodes(n_B)); % the local node index of n_B inside the element 
% %             Elem_Nodes=Nodes(Elements_cmgui(elements_contain_nodes(m_el),:),:);
% %             for n_GP=1:GP_per_elem_dir^3
% %                 ksi_location_3D=GPcoords(n_GP,:);
% % %                 N_at_ksi=Shapefunction_lagrange(ksi_location,Nodes_per_elem_dir);
% %                 thetaN_thetaksi=Shapefunction_derivative(ksi_location_3D,Nodes_per_elem_dir);
% %                 grad_NA=thetaN_thetaksi(local_ind_nA_in3Delem,:)*inv(Elem_Nodes.'*thetaN_thetaksi);
% %                 grad_NB=thetaN_thetaksi(local_ind_nB_in3Delem,:)*inv(Elem_Nodes.'*thetaN_thetaksi);
% %                 dV=det(Elem_Nodes.'*thetaN_thetaksi);
% % %                 K_mat_FF(n_A,n_B)=K_mat_FF(n_A,n_B)+dot(grad_NA,grad_NB)*dV*GPweights_tot(n_GP,:);% athroisma kata GP isodynamei me oloklirwsi sto element, athroisma kata m_el simainei pws lamvanw oles tis synistwses sto K_AB (stiffness matrix anamesa se komvo A kai B) ypopsin
% %                 v_mat_FF(n_B)=v_mat_FF(n_B)+dot(grad_NA,grad_NB)*dV*GPweights_tot(n_GP,:);
% %             end
% % % %           
% %         end
% %         
% %     end
% %          
% %     K_mat_FF(n_A,:)=v_mat_FF;
% %     
% %     %% K_FB
% %     v_mat_FB=zeros(1,length(Bound_nodes)); % preallocate for forloop
% %     for n_B=1:length(Bound_nodes) % this part is not symmetric - must loop over all Bound_Nodes
% %         elements_contain_nodes=[];
% %         if Bound_nodes(n_B)~=Free_nodes(n_A)
% %             elements_contain_nodes=find(sum(ismember(Elements_cmgui,[Free_nodes(n_A),Bound_nodes(n_B)]),2)==2); %vres ta elements pou periexoun kai tous dyo komvous
% %         else
% %             elements_contain_nodes=find(sum(ismember(Elements_cmgui,Free_nodes(n_A)),2)==1);
% %         end
% %         for m_el=1:length(elements_contain_nodes)
% %             local_ind_nA_in3Delem=find(Elements_cmgui(elements_contain_nodes(m_el),:)==Free_nodes(n_A)); % the local node index of n_A inside the element 
% %             local_ind_nB_in3Delem=find(Elements_cmgui(elements_contain_nodes(m_el),:)==Bound_nodes(n_B)); % the local node index of n_B inside the element 
% %             Elem_Nodes=Nodes(Elements_cmgui(elements_contain_nodes(m_el),:),:);
% %             for n_GP=1:GP_per_elem_dir^3
% %                 ksi_location_3D=GPcoords(n_GP,:);
% % %                 N_at_ksi=Shapefunction_lagrange(ksi_location,Nodes_per_elem_dir);
% %                 thetaN_thetaksi=Shapefunction_derivative(ksi_location_3D,Nodes_per_elem_dir);
% %                 grad_NA=thetaN_thetaksi(local_ind_nA_in3Delem,:)*inv(Elem_Nodes.'*thetaN_thetaksi);
% %                 grad_NB=thetaN_thetaksi(local_ind_nB_in3Delem,:)*inv(Elem_Nodes.'*thetaN_thetaksi);
% %                 dV=det(Elem_Nodes.'*thetaN_thetaksi);
% %                 if dV<=0 
% %                     disp('error in solve_Laplace_between_oppos_boundaries_new - negative dV!');
% %                 end
% % %                 K_mat_FF(n_A,n_B)=K_mat_FF(n_A,n_B)+dot(grad_NA,grad_NB)*dV*GPweights_tot(n_GP,:);% athroisma kata GP isodynamei me oloklirwsi sto element, athroisma kata m_el simainei pws lamvanw oles tis synistwses sto K_AB (stiffness matrix anamesa se komvo A kai B) ypopsin
% %                 v_mat_FB(n_B)=v_mat_FB(n_B)+dot(grad_NA,grad_NB)*dV*GPweights_tot(n_GP,:);
% %             end
% % 
% %             
% %         end
% %         
% %     end
% %     
% %      
% % %     malakia_ole(nA,1:length(Bound_nodes))=v_mat_FF;
% %     K_mat_FB(n_A,:)=v_mat_FB;
% % end
% % now fill the rest of the places of the symmetric K_mat: -- no need for
% this now, as I go element by element.
% K_mat_FF=K_mat_FF+K_mat_FF.'-diag(diag(K_mat_FF)); % complete the symmetric matrix,  afaireis to diag dioti to prostheses 2 fores --no need for this now I am building each node separately

toc
disp(['time required to build K_mat : ',num2str(toc)])

inv(K_mat_FF);
toc
disp(['time required to invert K_mat : ',num2str(toc)])



%%  step 4:construct the final residual forces vector: R_fin=R_F- K_FB* PHI_B
% no need


%% solve and find PHI_F=inv(K_FF) -- here was the bug, had forgotten the (-) before..

Phi_U=-inv(K_mat_FF)*(K_mat_FB*Phi_B); % tried K_mat_FF\(K_mat_FB*Phi_B) and had same result
toc
disp(['time required to for inv(K_mat_FF)*(K_mat_FB*Phi_B) is : ',num2str(toc)])
Phi_U_CG=-pcg(K_mat_FF,(K_mat_FB*Phi_B),10^(-7),1000); % use this 
toc
disp(['time required for pcg on K_FF^-1*(K_FB*Phi_B) is : ',num2str(toc)])

%% reorder \Phi into node indices:

Phi_Nodes_inv(Free_nodes,:)=Phi_U;
Phi_Nodes_inv(Bound_nodes,:)=Phi_B;
Phi_Nodes_CG(Free_nodes,:)=Phi_U_CG;
Phi_Nodes_CG(Bound_nodes,:)=Phi_B;

%% find gradient \Phi at nodes: -- errors can be big at nodes( gradients between elements: use discontinuous mesh instead?
for n_nod=1:size(Nodes,1) % for quadratic elements a for loop takes 30 mins, whereas the parfor takes 5 mins
    dist_inv=0; dist_CG=0; % this for comparison between gradients of node B calculated at different elements containing node B
    elements_contain_nodes=[]; grad_Phi_mel_inv=[]; grad_Phi_mel_CG=[];
    elements_contain_nodes=find(sum(ismember(Elements_cmgui,n_nod),2)>0); 
    grad_Phi_mel_inv=zeros(length(elements_contain_nodes),3); grad_Phi_mel_CG=zeros(length(elements_contain_nodes),3); % preallocate
        
    for m_el=1:length(elements_contain_nodes)
        local_nod_ind_in3Delem=find(Elements_cmgui(elements_contain_nodes(m_el),:)==n_nod);
        Elem_Nodes=Nodes(Elements_cmgui(elements_contain_nodes(m_el),:),:);
        Phi_Elem_Nodes=Phi_Nodes_inv(Elements_cmgui(elements_contain_nodes(m_el),:),:);
        Phi_Elem_Nodes_CG=Phi_Nodes_CG(Elements_cmgui(elements_contain_nodes(m_el),:),:);
%         [ksi_location,Location_0,error,flag_corrected_nodal_coords]=find_ksi_location_from_XYZ_coords_with_Newton_solve_in_3D_elem(Nodes(n_nod,1),Nodes(n_nod,2),Nodes(n_nod,3),Elem_Nodes,10^(-9),Nodes_per_elem_dir,10,10,1); --% this can be done faster by a bespoke function I created
        ksi_location=find_ksi_location_of_local_node_ID_in_3D_element_cmgui_ordering(local_nod_ind_in3Delem,Nodes_per_elem_dir);

%         if error==1
%             disp('error in finding ksi_location');
%             disp(['ksi_location : ',num2str(ksi_location),' Location_0 : ',num2str(Location_0),'original location : ',num2str(Nodes(n_nod,:)),' flag_corrected_nodal_coords ; ',num2str(flag_corrected_nodal_coords)]);
%         end
        thetaN_thetaksi=Shapefunction_derivative(ksi_location,Nodes_per_elem_dir);
        grad_Phi_mel_inv(m_el,:)=(((Phi_Elem_Nodes.'*thetaN_thetaksi)*inv(Elem_Nodes.'*thetaN_thetaksi))/norm((Phi_Elem_Nodes.'*thetaN_thetaksi)*inv(Elem_Nodes.'*thetaN_thetaksi))).'; % row vector
        grad_Phi_mel_CG(m_el,:)=(((Phi_Elem_Nodes_CG.'*thetaN_thetaksi)*inv(Elem_Nodes.'*thetaN_thetaksi))/norm((Phi_Elem_Nodes_CG.'*thetaN_thetaksi)*inv(Elem_Nodes.'*thetaN_thetaksi))).';% row vector
        
        if flag_discon_elements==1 % create discontinous mesh
            Discon_mesh_global_nod_ind=(elements_contain_nodes(m_el)-1)*Nodes_per_elem_dir^3+local_nod_ind_in3Delem;
            Grad_Phi_Nodes_CG_discont(Discon_mesh_global_nod_ind,:)=grad_Phi_mel_CG(m_el,:);
            Discon_mesh_Nodes(Discon_mesh_global_nod_ind,:)=Nodes(Elements_cmgui(elements_contain_nodes(m_el),local_nod_ind_in3Delem),:);
            % sanity check:
            if norm(Nodes(n_nod,:)-Discon_mesh_Nodes(Discon_mesh_global_nod_ind,:))>10^(-10)
                disp('error in solve_Laplace_between_oppos_boundaries.. function these should be exactly equal');
            end
            Discon_mesh_Elements(elements_contain_nodes(m_el),local_nod_ind_in3Delem)=Discon_mesh_global_nod_ind;
        end
        %% ayto den einai swsto giati sygkrinei mono ena element me ta proigoumena oxi ola ta gradients sto node gia kathe element metaksy tous ((m_el-1)! sto synolo sygkriseis)
% %         if m_el>1 % sygkrine ta gradient se kathe element pou periexei to node kai report it an yparxoun apokliseis..
% %             dist_prox=abs(norm(grad_Phi_mel(:,m_el)-grad_Phi_mel(:,m_el-1)));
% %             dist_prox_CG=abs(norm(grad_Phi_mel_CG(:,m_el)-grad_Phi_mel_CG(:,m_el-1)));
% %             if dist_prox>dist
% %                 dist=dist_prox;
% %             end
% %             if dist_prox_CG>dist_CG
% %                 dist_CG=dist_prox_CG;
% %             end
% %         end
    end

    for m_el=1:(length(elements_contain_nodes)-1) % compare gradPhi at element m_el with all the following ones
        dist_CG_mat=[];
        dist_CG_mat=grad_Phi_mel_CG((m_el+1):length(elements_contain_nodes),:)-[grad_Phi_mel_CG(m_el,1)*ones(length(elements_contain_nodes)-m_el,1),grad_Phi_mel_CG(m_el,2)*ones(length(elements_contain_nodes)-m_el,1),grad_Phi_mel_CG(m_el,3)*ones(length(elements_contain_nodes)-m_el,1)]; % rows: length(elements_contain_nodes)-1 , cols:3
        if max(sqrt(sum(dist_CG_mat.^2,2)))>dist_CG
            dist_CG=max(sqrt(sum(dist_CG_mat.^2,2)));
        end
        
        dist_inv_mat=[];
        dist_inv_mat=grad_Phi_mel_inv((m_el+1):length(elements_contain_nodes),:)-[grad_Phi_mel_inv(m_el,1)*ones(length(elements_contain_nodes)-m_el,1),grad_Phi_mel_inv(m_el,2)*ones(length(elements_contain_nodes)-m_el,1),grad_Phi_mel_inv(m_el,3)*ones(length(elements_contain_nodes)-m_el,1)]; % rows: length(elements_contain_nodes)-1 , cols:3
        if max(sqrt(sum(dist_inv_mat.^2,2)))>dist_inv
            dist_inv=max(sqrt(sum(dist_inv_mat.^2,2)));
        end
    end
    
        
    Grad_Phi_Nodes_inv(n_nod,:)=mean(grad_Phi_mel_inv,1).'/norm(mean(grad_Phi_mel_inv,1));
    Grad_Phi_dist_inv(n_nod,:)=dist_inv; % to check quality afterwards. If there are huge gaps in gradient estimation betweeen elements something needs to be done!
    Grad_Phi_Nodes_CG(n_nod,:)=mean(grad_Phi_mel_CG,1).'/norm(mean(grad_Phi_mel_CG,1)); % the grad phi estimations from Phi estimation with conjugate gradient method --uSE THIS!
    Grad_Phi_dist_CG(n_nod,:)=dist_CG; % to check quality afterwards. If there are huge gaps in gradient estimation betweeen elements something needs to be done!
end


%% plot results
fig_ind=[];
if flag_plot>0 % initialise fig_ind depending on input or lack thereof
    fig_ind=fig_ind_in;
    if isempty(fig_ind_in)==1
        fig_ind=0;
    end
end
if flag_plot==2 % plot resutls from both pcg and inv(K) --not generally needed only for my report on how I solved this
    fig_ind=fig_ind+1;
    figure(fig_ind);
    plot_cmgui_ordered_mesh(Elements_cmgui,Nodes,Nodes_per_elem_dir,3);
    hold on;
    quiver3(Nodes(:,1),Nodes(:,2),Nodes(:,3),Grad_Phi_Nodes_inv(:,1),Grad_Phi_Nodes_inv(:,2),Grad_Phi_Nodes_inv(:,3),1);
    title('grad(Phi) -inv(K)')
    hold off;
    fig_ind=fig_ind+1;
    figure(fig_ind);
    scatter3(Nodes(:,1),Nodes(:,2),Nodes(:,3),20*ones(size(Phi_Nodes_inv)),Phi_Nodes_inv,'f');
    colorbar;
    hold on;
    title('Phi -inv(K)')
    hold off;
    fig_ind=fig_ind+1;
    figure(fig_ind);
    scatter3(Nodes(:,1),Nodes(:,2),Nodes(:,3),20*ones(size(Phi_Nodes_CG)),Phi_Nodes_CG,'f');
    colorbar;
    hold on;
    title('Phi -pcg(..,..,10^-7, 1000)')
    hold off;
    fig_ind=fig_ind+1;
    figure(fig_ind);
    plot_cmgui_ordered_mesh(Elements_cmgui,Nodes,Nodes_per_elem_dir,3);
    hold on;
    quiver3(Nodes(:,1),Nodes(:,2),Nodes(:,3),Grad_Phi_Nodes_CG(:,1),Grad_Phi_Nodes_CG(:,2),Grad_Phi_Nodes_CG(:,3),1);
    title('grad(Phi) -pcg(..,..,10^-7, 1000)');
    hold off;
elseif flag_plot==1 % just plot results from pcg
    fig_ind=fig_ind+1;
    figure(fig_ind);
    scatter3(Nodes(:,1),Nodes(:,2),Nodes(:,3),20*ones(size(Phi_Nodes_CG)),Phi_Nodes_CG,'f');
    colorbar;
    hold on;
    title('Phi -pcg(..,..,10^-7, 1000)')
    hold off;
    fig_ind=fig_ind+1;
    figure(fig_ind);
    plot_cmgui_ordered_mesh(Elements_cmgui,Nodes,Nodes_per_elem_dir,3);
    hold on;
    quiver3(Nodes(:,1),Nodes(:,2),Nodes(:,3),Grad_Phi_Nodes_CG(:,1),Grad_Phi_Nodes_CG(:,2),Grad_Phi_Nodes_CG(:,3),1);
    title('grad(Phi) -pcg(..,..,10^-7, 1000)');
    hold off;    
end
if nargout==7
    fig_ind_out=fig_ind;
end
if nargout>3
    varargout{1}=Grad_Phi_Nodes_CG_discont;
    varargout{2}=Discon_mesh_Nodes;
    varargout{3}=Discon_mesh_Elements;
    varargout{4}=fig_ind_out;
end
    