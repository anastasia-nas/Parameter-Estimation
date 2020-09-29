function [diast_pressures_Pa,Diastolic_Frames_Indices,Diastolic_Frames_Volumes_ml,flag_reference_frame_negative_press]=read_corrected_pressures_from_catheter_data_new_cases(diast_frames_info_folder,chosen_offset,reference_frame)
% created 1- June_2019
% function that reads corrected pressures (by creating an offset necessary
% to remove negative pressures in the "passive inflation" window as I define it
% PURPOSE: to use a function instead of a section that is copied from script to script, for more consistency and to avoid bugs.
% INPUT: diast_frames_info_folder: folder where the synced_frame_pressure_volume_data_offset_*_msec.txt files are stored from
% my matlab script
%NOTE: if chosen reference frame is not part of the diastolic indices as
%selected by the pV syncing script (this can only occur when I'm choosing as a reference frame a frame prior to the one with the minimun pressure, which only happens when I'm testing
% for the sensitivity of the pipeline to the reference frame (+/- 1 reference frame)), I just add it to the beginning of thediastolic indices to proceed with the analysis

flag_reference_frame_negative_press=0; % this flag will return an error if the reference frame has negative pressure and should be removed from analysis
list2=dir([diast_frames_info_folder,'/synced_frame_pressure_volume_data_offset_*_msec.txt']);
if max(size(list2))>1
    disp('WARNING 2 ! synced data exist for more than 2 offset values - make sure you use the right one!');
    synced_data_file_dirpath=[diast_frames_info_folder,'/synced_frame_pressure_volume_data_offset_',num2str(chosen_offset),'_msec.txt']; 
else
    synced_data_file_dirpath=[diast_frames_info_folder,'/',list2.name];
end
fid=fopen(synced_data_file_dirpath);
data_per_frame_cell=textscan(fid,'%f  %f  %f  %f','HeaderLines',1);
fclose(fid);

data_per_frame=cell2mat(data_per_frame_cell); % 1st col: frame id, 2nd col: time, 3rd col:press (mmHg), 4th volume (ml)
% se kPa einai ta pressures vasika:
% pressures_mmHg=data_per_frame(:,3).';
frame_indices=data_per_frame(:,1).';
pressures_kPa=data_per_frame(:,3).';
volumes_ml=data_per_frame(:,4).';

% read the relevant diastolic frames:
diastolic_frames_file_dirpath=[diast_frames_info_folder,'/diastolic_frames_from_offset_',num2str(chosen_offset),'_msec.txt']; 
fid=fopen(diastolic_frames_file_dirpath);
diast_indices=cell2mat(textscan(fid,'%d','HeaderLines',1)); 
fclose(fid);
reference_frame_ind=find(diast_indices==reference_frame);
if isempty(reference_frame_ind)==1% if reference frame doesn't belong to diastolic indices (this should happen when I choose as reference frame a frame before the one corresponding to min pressure, which occurs when I test for the effect of the choice of reference frame)
    diast_indices_pre=[reference_frame;diast_indices]; % I do this so that the reference frame is added to the diastolic frames -- added 11 October 2019
    
else
    diast_indices_pre=diast_indices(reference_frame_ind:end);
end
% just pressure data for each diastolic frame:
% diast_pressures_mmHg=pressures_mmHg(diast_indices);
diast_pressures_kPa_pre=pressures_kPa(diast_indices_pre); % don't use this one before you get rid of the negative pressures!!!
Diastolic_Frames_Volumes_pre=volumes_ml(diast_indices_pre);
% twra petakse ta frames sta opoia exeis arnitiki piesi:
coun_d=0;
for n_d=1:max(size(diast_indices_pre))
    if diast_pressures_kPa_pre(n_d)>=0
        coun_d=coun_d+1;
        Diastolic_Frames_Indices(coun_d)=diast_indices_pre(n_d);
        diast_pressures_PRE(coun_d)=diast_pressures_kPa_pre(n_d);
        Diastolic_Frames_Volumes_ml(coun_d)=Diastolic_Frames_Volumes_pre(n_d);
    else
       disp('SOS:-----------------------------------------')
       disp('in script read_corrected_pressures_from_catheter_data_new_cases.m')
       disp(['identified negative pressure in diastolic frame: ',num2str(diast_pressures_kPa_pre(n_d)),' !'])
       disp('don''t use this frame as reference frame');
       
    end
end

if reference_frame<Diastolic_Frames_Indices(1)
    flag_reference_frame_negative_press=1;
end
% diast_pressures_Pa=diast_pressures_mmHg*133.322368;
% diast_pressures_kPa=diast_pressures_Pa/1000; % you won't need this..

diast_pressures=diast_pressures_PRE-diast_pressures_PRE(1); % thewrw pws sto reference frame exw 0 pressure.
diast_pressures_Pa=diast_pressures*1000;
    