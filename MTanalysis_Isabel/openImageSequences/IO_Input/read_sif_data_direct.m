function [imageData] = read_sif_data_direct(filename_fluorescent, SeriesLength, ImageSize, start_slice_no, end_slice_no)
%     
% execution for extacting one or more frame from current sif file
%
% filename_fluorescent
%        -- file name
% SeriesLength
%        -- sequence length
% ImageSize
%        -- slice size
% start_slice_no
%        -- start index of sequence
% end_slice_no
%        -- end   index of sequence
% 
% imageData
%        -- structure for saving data
%

% protection
imageData = [];
if end_slice_no<start_slice_no
    return;
end
if end_slice_no>SeriesLength
    return;
end
if start_slice_no<0
    return;
end

% initial parameter
% this is fixed for sif  format
image_byte_const = 4; 

% whole sequence
for image_data_shift_frame=(start_slice_no-1):(end_slice_no-1)

    % total sif file size 
    s = dir(filename_fluorescent); 
    observation_file_size = s.bytes;

    % the size of header from current file
    estimate_image_data_size = ImageSize*SeriesLength*image_byte_const;
    computed_head_size = observation_file_size - estimate_image_data_size;
    computed_head_size = computed_head_size;
%     disp(' ')
%     disp([' computed_head_size ' num2str(computed_head_size) ' from ' filename_fluorescent]);
%     disp(' ')

    % open the file
    f=fopen(filename_fluorescent,'r');
    if f < 0
        error('Could not open the file.');
    end
    if ~isequal(fgetl(f),'Andor Technology Multi-Channel File')
    fclose(f);
    error('Not an Andor SIF image file.');
    
end

% access file 
N = computed_head_size+image_data_shift_frame*ImageSize*image_byte_const;
fseek(f,(N),'bof');
t = fread(f,ImageSize,'single=>single');

% show image
imageData{image_data_shift_frame+1}.sliceData = reshape(t,[sqrt(ImageSize) sqrt(ImageSize)]);

% %     figure(1)
% %     imshow(imrotate(imageData,90),[]);
% disp(' ');
% disp(['extracted data for image: ' num2str(image_data_shift_frame)]);
% disp(' ');

% close the file
fclose(f);

end

% disp(' ');
% disp(' finished ... ');
% disp(' ');

end 
