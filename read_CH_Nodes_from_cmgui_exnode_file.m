function CH_Nodes_of_Fields_in_file=read_CH_Nodes_from_cmgui_exnode_file(file_location)


% created by Anastasia on 23-05-2019
%-------------------------------------------------------------------------------------------------------------------------
% Purpose: 
%-------------------------------------------------------------------------------------------------------------------------
% read cmgui nodes from .exnode file
%-------------------------------------------------------------------------------------------------------------------------
% Notes on methodology:
%-------------------------------------------------------------------------------------------------------------------------
% each node contains values according to the cmgui ordering:
% (d/ds1,d/ds2,d2/ds1ds2,d/ds3,d2/ds1ds3,d2/ds2ds3,d3/ds1ds2ds3) per
% coordinate direction
%-------------------------------------------------------------------------------------------------------------------------
% Input:
%-------------------------------------------------------------------------------------------------------------------------
% file_location: location of .exnode file

%-------------------------------------------------------------------------------------------------------------------------
% Output:
%-------------------------------------------------------------------------------------------------------------------------
% CH_Nodes_of_Fields_in_file: structure composed of CH Nodes per field, for each field it's a matrix size( #Nodes, 8*3) , 8: nodal values per coordinate, 3: number of coordinates in 3D
% ordered first along the nodal value and derivatives: (d/ds1,d/ds2,d2/ds1ds2,d/ds3,d2/ds1ds3,d2/ds2ds3,d3/ds1ds2ds3) and then per coordinate direction: 8,8,8 per row 
%-------------------------------------------------------------------------------------------------------------------------
% SOS points:
%-------------------------------------------------------------------------------------------------------------------------
% output  contains field name e.g. CH_Nodes_of_Fields_in_file.coordinate
% nodal values are ordered according to cmgui ordering for cubic hermite: (d/ds1,d/ds2,d2/ds1ds2,d/ds3,d2/ds1ds3,d2/ds2ds3,d3/ds1ds2ds3)
%=========================================================================================================================








% for testing:: --------------------------------------------------------------------------------------------------------------------
% clear all; close all; clc;
% file_location='/staff/an11/new_cases_meshed/new_cases_hole/Cubic_Hermite_Meshes_From_Pablos_Toolkit/case_BAL/MNaNMesh9124_07.exnode'
% --------------------------------------------------------------------------------------------------------------------------------------

fid=fopen(file_location);
no_of_fields=[];
coun_read_field_in_node=[];
while feof(fid)==0 % checks if an operator has set the end of file indicator
    line_ex=fgetl(fid); % reads next line
    % find number of fields in file:
   
    % read number of fields and coordinates per field:
    if isempty(regexp(line_ex,'#Fields=','match'))==0
        no_of_fields=str2num(line_ex(regexp(line_ex,'#Fields=','end')+1));
    end
    
    % identify the possibility that a new field is introduced and read how
    % many coordinates correspond to the field
   
    for n_field=1:no_of_fields
        field_flag_str=[num2str(n_field),') '];
        if isempty(regexp(line_ex,field_flag_str))==0 && isempty(regexp(line_ex,'coordinate'))==0 && isempty(regexp(line_ex,'#Components'))==0 % these strings identify  the line that explains each field characteristics in exnode file
            field_components(n_field,:)=str2num(line_ex(regexp(line_ex,'#Components=','end')+1));% the number of field coordinates
            positions_with_comma_in_line=regexp(line_ex,',','start');
            field_title(n_field,:)=line_ex(regexp(line_ex,'coordinates,','end')+2:positions_with_comma_in_line(2)-1);
        end
    end
    
    % if we are now reading nodes: (this occurs after all fields have been
    % introduced
    if isempty(regexp(line_ex,'Node:'))==0
        coun_read_field_in_node=0; % set this to zero after a new node has been identified.
        node_ID_pre=textscan(line_ex,'Node: %f');
        node_ID=node_ID_pre{1};
        field=[]; % diagrafw tyxon proigoumenes apothikeyseis gia na min ginoun lathi
        
    end
    if isempty(coun_read_field_in_node)==0 % an exw arxisei na diavazw nodes (opote exw teleiwse me ta introductory stuff kai eite diavazw lines me "Node:" eite lines pou yparxoun fields kai nodal values.
        
        if isempty(regexp(line_ex,'Node:'))==1 % na min vriskomai se line pou leei to node ID ("Node:  ")
            coun_read_field_in_node=coun_read_field_in_node+1; % ayto metraei line pou yparxoun nodal values (kai oxi line pou na exei to intro "Node:")
            % find which field and which component of the field corresponds to coun_read_field_in_node (the current field we are reading)
            for n_field=1:no_of_fields
               if coun_read_field_in_node<= sum(field_components(1:n_field,:)) && coun_read_field_in_node>=sum(field_components(1:(n_field-1),:))+1
                   field_number=n_field; 
                   coord_no=coun_read_field_in_node-sum(field_components(1:(n_field-1),:));
               end
            end
            mpou=sscanf(line_ex,'%f').';
            field(coord_no,:)=reshape(mpou,[1,length(mpou)]);% turn into row in case it isn't
            % store all field coordinates together
%             if coord_no==field_components(field_number,:) % an i synistwsa einai i max synistwsa tou field apothikeyse to
%                 CH_Nodes_of_Fields_in_file.(field_title)(node_ID,:)=reshape(field.',[1,size(field,1)*size(field,2)]);
%             end

            % decided to store each coordinate separately, as will be easier to read and interpolate.
            CH_Nodes_of_Fields_in_file.([field_title,'_',num2str(coord_no)])(node_ID,:)=field(coord_no,:);
        end
        
    end
    % read values per node: identify "Node:" in line
    % general format is: 
    % field (x,y,z,f,s,n,p whatever) and associated values per row
    % so for CH mesh there are 8 values per row.
     
    
end
fclose(fid);