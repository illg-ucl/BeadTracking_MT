function [pathname_input  filename_input fileformat_input current_image head_fluorecent] = GuiFileInputSequenceEntry(hObject,handles)
%
% I/O input file from hard disk
%
% Input: 
%   hObject -- Object  of matlab
%   handles -- handles of matlab
% Output:
%   pathname_input -- file of path name
%   filename_input -- file of file name
%   fileformat_input -- file of format
%   current_image  -- current image to be shown, which corresponds to
%   handles.current_frame_number
%   head_fluorecent-- head information
%
% Note:
%   which will be supported by the following format:
%   - sif
%   - avi
%   - tif
%    

% get file name with special directory
[filename_input,pathname_input] = uigetfile({'*.*'}, 'my file','MultiSelect','on',handles.default_pathname);
if isequal(filename_input,0)
    set(handles.text_status,'string','no files have been selected here');
    return;
end

% protection: check file format and ignor uncorrect format
fileformat_input = filename_input(end-2:end); % 'sif'/'tif'/'avi'
if ~( (isequal(fileformat_input,'sif')) || (isequal(fileformat_input,'tif')) || (isequal(fileformat_input,'avi')) )
    set(handles.text_status,'string','file format can NOT be accessed!'), pause(0.05);    
    return;
end

% create three folders automatically to save the results
subCreateDirectory(pathname_input);

% I/O sif file (hard-disk)   
if isequal(filename_input(end-2:end),'sif') 

    % sif format
    %
    bright_field_image = handles.filetype;
    [data_fluorecent,head_fluorecent] = sifread2009_singleframe(strcat(pathname_input,filename_input),handles.current_frame_number,bright_field_image);
   
    
    % file info from its header
    [ ...
    ReturnCode, ...         % whether successful precedure
    SeriesLength, ...       % total frame number
    ImageSize, ...          % size of one frame
    TotalAcquisitionSize ...%(total frame number)*(size of one frame) 
    ] ...
    = GetAndorSifSize(strcat(pathname_input,filename_input),0); 

    % reflesh GUI interface 
    set(handles.edit_head_of_sequence,'string',num2str(handles.current_frame_number));
    set(handles.edit_tail_of_sequence,'string',num2str(SeriesLength));

    % current frame number
    current_frame_number = handles.current_frame_number;
        
    % this is used for the whole video
    start_slice_no = current_frame_number;
    end_slice_no   = current_frame_number; 

    % current image
    [imageData]     = read_sif_data_direct(strcat(pathname_input,filename_input), SeriesLength, ImageSize, start_slice_no, end_slice_no);
    % MAYBE this is a bug from Andor SOLIS and Matlab
    current_image   = imrotate(imageData{current_frame_number}.sliceData,90); 

   
elseif isequal(filename_input(end-2:end),'tif') 

    % tif format
    [stack, nbImages] = tiffread2(strcat(pathname_input,filename_input));
    imshow( double(stack(handles.current_frame_number).data), [] );
    
     % attribute from current sif file
    handles.im_para.width    = stack(1).width;
    handles.im_para.height   = stack(1).height;
    handles.im_para.length   = nbImages;
    handles.im_para.interval = 0.04; % internal const parameter because tiff header can not provide this information
    guidata(hObject,handles);
    set(handles.text_status,'string',['current frame number is: ' num2str(handles.current_frame_number) ...
        ' width is ' num2str(handles.im_para.width) ...
        ' height is ' num2str(handles.im_para.height) ...
        ' length is ' num2str(handles.im_para.length) ...
        ' interval is ' num2str(handles.im_para.interval)]);    
        
elseif isequal(filename_input(end-2:end),'avi') 

    % avi format
    % file Input
    [mov readerobj] = AviInput(strcat(pathname_input,filename_input));
    if isempty(mov)
        return;
    end
    
    % reflesh output parameters
    pathname_input   = pathname_input;   % 1
    filename_input   = filename_input;   % 2
    fileformat_input = fileformat_input; % 3
    current_image    = mov;              % 4
    head_fluorecent  = readerobj;        % 5

    
else

    set(handles.text_status,'string','file format can NOT be accessed!'); pause(0.1);            

end


    function subCreateDirectory(pathname_input)
    %
    % create three directories for saving the output files
    %

        current_path = pwd;
        cd(pathname_input);    
        if ~exist('export_all_intensity','dir') 
            mkdir(strcat(pathname_input,'export_all_intensity'));        
        end
        if ~exist('export_current_image','dir')
            mkdir(strcat(pathname_input,'export_current_image'));        
        end
        if ~exist('export_video','dir')
            mkdir(strcat(pathname_input,'export_video'));
        end
        cd(current_path);

    end % this function


end % this function