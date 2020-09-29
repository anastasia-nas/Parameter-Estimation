function Wext_pV_per_DF=calc_Wext_with_pV(diastolic_frames_list,positive_pressures_at_diastolic_frames_without_ref,Ref_cavity_volume,Def_cavity_volume_struct)
% 
% created by Anastasia on 18-04-2020
%-------------------------------------------------------------------------------------------------------------------------
% Purpose: 
%-------------------------------------------------------------------------------------------------------------------------

%-------------------------------------------------------------------------------------------------------------------------
% Notes on methodology:
%-------------------------------------------------------------------------------------------------------------------------
% Wext is calculated as \int_V0^Vend p*dV
%-------------------------------------------------------------------------------------------------------------------------
% Input:
%-------------------------------------------------------------------------------------------------------------------------
% diastolic_frames_list: this needs to be all frames from reference frame on e.g. [27,28,1,2,3] - EXCLUDING the reference frame
% positive_pressures_at_diastolic_frames_without_ref: POSITIVE for INFLATION! corresponding pressure to each frame: e.g. [1500,1600,1700,1900,2000]-excluding pressure 
%                                 at reference frame which is assumed to be zero!
% Ref_Nodes: mesh nodes at reference frame
% Elements_cmgui:elements with with cmgui ordering (no priority at corner nodes)
% Def_cavity_volume_struct: a structure containing the cavity volume of the deformed mesh at each frame (e.g. Def_cavity_volume_struct.Frame_25,Def_cavity_volume_struct.Frame_27, Def_cavity_volume_struct.Frame_1 etc..
%                                   can be calculated with "calculate_cavity_volume_Lagrange_meshes_with_hole.m"
% Ref_cavity_volume: cavity volume at reference frame (calculated with calculate_cavity_volume_Lagrange_meshes_with_hole.m)
%-------------------------------------------------------------------------------------------------------------------------
% Output:
%-------------------------------------------------------------------------------------------------------------------------
% Wext_pV_per_DF: pV based external work calculation per diastolic frame

%-------------------------------------------------------------------------------------------------------------------------
% SOS points:
%-------------------------------------------------------------------------------------------------------------------------

%% ==============================================START==============================================================

for n_DF=1:length(diastolic_frames_list)
    def_volume(n_DF,:)=Def_cavity_volume_struct.(['Frame_',num2str(diastolic_frames_list(n_DF))]);
    press(n_DF,:)=-positive_pressures_at_diastolic_frames_without_ref(n_DF); % I reverse sign here because I want positive pressure in input to correspond to inflation (opposite sign to da vector!)        
end
pressures_data=[0;press]; % including reference frame
volumes_data=[Ref_cavity_volume;def_volume];

for n_DF=1:length(diastolic_frames_list) % note n_DF runs 1 step ahead of pressures data and volumes data because they also include the reference frame. (that's why the indices n_DF go +1 inside the pressure and volumem data structures
    Wext_pV_per_DF(n_DF)=sum(0.5*(pressures_data(1:(n_DF))+pressures_data(2:(n_DF+1))).*(volumes_data(2:(n_DF+1))-volumes_data(1:n_DF))); % n_DF in this contexts corresponds to previous step and n_DF+1 to current step cause pressure & volume data include the reference
end



