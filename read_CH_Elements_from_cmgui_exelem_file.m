function [CH_Elements,CH_Elements_scale_factors]=read_CH_Elements_from_cmgui_exelem_file(file_location)
% file_location=CH_mesh_exelem_file;

% created by Anastasia on 23-05-2019
%-------------------------------------------------------------------------------------------------------------------------
% Purpose: 
%-------------------------------------------------------------------------------------------------------------------------
% read cmgui nodes from .exelem file
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
% I assume that scale factor indices follow this order (prioritisation of nodal derivatives per node:
% 1.  #Values= 8
%        Value indices: 1 2 3 4 5 6 7 8 
%        Scale factor indices: 1 2 3 4 5 6 7 8 
%       2.  #Values= 8
%        Value indices: 1 2 3 4 5 6 7 8 
%        Scale factor indices: 9 10 11 12 13 14 15 16 
%       3.  #Values= 8
%        Value indices: 1 2 3 4 5 6 7 8 
%        Scale factor indices: 17 18 19 20 21 22 23 24 
%       4.  #Values= 8
%        Value indices: 1 2 3 4 5 6 7 8 
%        Scale factor indices: 25 26 27 28 29 30 31 32 
%       5.  #Values= 8
%        Value indices: 1 2 3 4 5 6 7 8 
%        Scale factor indices: 33 34 35 36 37 38 39 40 
%       6.  #Values= 8
%        Value indices: 1 2 3 4 5 6 7 8 
%        Scale factor indices: 41 42 43 44 45 46 47 48 
%       7.  #Values= 8
%        Value indices: 1 2 3 4 5 6 7 8 
%        Scale factor indices: 49 50 51 52 53 54 55 56 
%       8.  #Values= 8
%        Value indices: 1 2 3 4 5 6 7 8 
%        Scale factor indices: 57 58 59 60 61 62 63 64 
%=========================================================================================================================








% % for testing:: --------------------------------------------------------------------------------------------------------------------
% clear all; close all; clc;
% file_location='/staff/an11/new_cases_meshed/new_cases_hole/Cubic_Hermite_Meshes_From_Pablos_Toolkit/case_BAL/MNaNMesh9124_07.exelem';
% % --------------------------------------------------------------------------------------------------------------------------------------

fid=fopen(file_location);
no_of_fields=[];
intro_stuff_is_over=0;
coun_lines_after_nodes_read=0.1; % for initialisation to allow counters to work but make sure value can never be 2 unless the "coun_lines_after_nodes_read=1" line is activated
coun_lines_after_scale_factors_read=0.1; % for initialisation to allow counters to work but make sure value can never be 2 unless the "coun_lines_after_scale_factors_read=1" line is activated
% first_Nodes_line=0; % if #Nodes=  8 line has been encountered in text -> removed it as I don't use it
while feof(fid)==0 % checks if an operator has set the end of file indicator
    line_ex=fgetl(fid); % reads next line
    % find number of fields in file:
   
    % read number of fields and coordinates per field:
    if isempty(regexp(line_ex,'#Fields=','match'))==0
        no_of_fields_pre=textscan(line_ex,'#Fields= %d');
        no_of_fields=no_of_fields_pre{1};
    end
    
    if isempty(regexp(line_ex,'#Nodes=','match'))==0
        no_of_nodes=textscan(line_ex,'#Nodes= %d');
%         first_Nodes_line=1; % we have read the nodes number per element -> removed it as I don't use it
    end
    
    % identify the possibility that a new field is introduced and read how
    % many coordinates correspond to the field
   
    for n_field=1:no_of_fields
        field_flag_str=[num2str(n_field),') '];
        if isempty(regexp(line_ex,field_flag_str))==0 && isempty(regexp(line_ex,'coordinate'))==0 && isempty(regexp(line_ex,'#Components'))==0 % these strings identify  the line that explains each field characteristics in exnode file
            field_components(n_field,:)=str2num(line_ex(regexp(line_ex,'#Components=','end')+1));
            positions_with_comma_in_line=regexp(line_ex,',','start');
            field_title(n_field,:)=line_ex(regexp(line_ex,'coordinates,','end')+2:positions_with_comma_in_line(2)-1);
        end
    end
    
    % if we are now reading elements: (this occurs after all the
    % introductory stuff is over
    if isempty(regexp(line_ex,'Element:'))==0
        nodes_line_read=0; % not yet read "Nodes:" line for this element
        scale_factors_line_read=0; % not yet read "Scale factors:" line for this element
        intro_stuff_is_over=1;
        element_ID_pre=textscan(line_ex,'Element: %d ');
        element_ID=element_ID_pre{1};
%         field=[]; % diagrafw tyxon proigoumenes apothikeyseis gia na min ginoun lathi
        
    end
    if intro_stuff_is_over==1 % an exw arxisei na diavazw elements (opote exw teleiwse me ta introductory stuff kai eite diavazw lines me 'Element:' eite lines pou yparxoun nodes per element
        coun_lines_after_nodes_read=coun_lines_after_nodes_read+1; % to demonstrate next line is now read
        coun_lines_after_scale_factors_read=coun_lines_after_scale_factors_read+1;
        if isempty(regexp(line_ex,'Nodes:'))==0 
            nodes_line_read=1;   
            coun_lines_after_nodes_read=1; % set to 1 after "Nodes:" is read
        end
        if isempty(regexp(line_ex,'Scale factors'))==0 
            scale_factors_line_read=1;
            coun_lines_after_scale_factors_read=1;
        end
        if isempty(regexp(line_ex,'Element:'))==1 % na min vriskomai se line pou na min leei "Node:" or "Scale factors"
            
            if nodes_line_read==1 && coun_lines_after_nodes_read==2 %to check this line follows immediately after Nodes: and that the current line is not nodes
                
                CH_Elements(element_ID,:)=sscanf(line_ex,'%d').';
                nodes_line_read=0; % reset flag to zero
            end
            if scale_factors_line_read==1 && coun_lines_after_scale_factors_read==2 %to check this line follows immediately after Nodes: and that the current line is not nodes
                
                CH_Elements_scale_factors(element_ID,:)=sscanf(line_ex,'%d').';
                scale_factors_line_read=0; % reset flag to zero
            end
            
        end
        
    end

     
    
end

fclose(fid);