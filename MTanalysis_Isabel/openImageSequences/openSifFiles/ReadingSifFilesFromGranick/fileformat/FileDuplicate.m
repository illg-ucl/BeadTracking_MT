function fileinfo=FileDuplicate(sourcefile,targetdirectory)

%fileinfo=FileDuplicate(sourcefile)
%
%FileDuplicate is a front end which unifies the commands necessary to
%duplicate a file. 
%
%INCLUDE:   FileDetails.m
%
%INPUTS:    SOURCEFILE: The name of the file to be duplicated,
%                       ('C:\path\thisfile.ext')
%           TARGETDIRECTORY:    The path to copy the file to.  (Optional)
%                               If not specified, the copied file will be
%                               placed in 'path\modified\'
%OUTPUTS:   FILEINFO:   A structure containing pertinent information about
%                       the file. Some are universal, regardless of file
%                       format, while others are only relevant to this file
%                       format. Refer to FileDetails for more details.
%
%Written by Stephen Anthony 6/27/2005 U. Illinois Urbana-Champaign
%Last modified by Scott Parker on 11/9/09

%Determine the type of file based upon extension
[directory,name,ext]=fileparts(sourcefile);

% Determine if a targetdirectory has been specified
if ~exist('targetdirectory','var')    
    %Generate the target directory if it does not exist
    targetdirectory=fullfile(directory,'modified');
    if exist(targetdirectory,'dir')~=7
        mkdir(targetdirectory);
    end
end

%Create the duplicate file
targetfile=fullfile(targetdirectory,[name ext]);
copyfile(sourcefile,targetfile)

%Determine the pertinent information, including that this is a copy
fileinfo=FileDetails(targetfile);
fileinfo.copy='Duplicate File';

