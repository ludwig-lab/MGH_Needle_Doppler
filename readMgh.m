function output = readMgh(filename, varargin)
% reads an mgh file in a FRAME BY FRAME manner.
% without option structure opt, the function reads the first Bscan.
% able to read multiple frames, defined by nFrames in the option structure.
% 
% varargin opt can provide: 
%     iBscan (bScanIndex to read)
%     metadata = [0 mghType nAlines(PerFrame) nZpixels nFrames]
%     singleFlag = 1 : read into single precision for gpu processing
%     nFrames (number of frames to read after iFrame including the iFrame
%
% %%% usage example :
% readOpt.iFrame = iFrame;
% readOpt.nFrames = nFrames;
% output = readMgh(fullfile2open, readOpt);
% %%% output is of a size (nZpixels)x(nAlines)x(nFrames), where 
% %%% output(:,:,1) is the "iFrame"th Bscan of the volume.
%
% << snam@mit.edu 2015930 >>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

complexFlag = 0;
singleFlag = 0;

if nargin == 1
    iFrame = 1;
    nFrames = 1;
    metadata = readMetaMgh(filename);
elseif nargin == 2; % optional structure given
    opt = varargin{1};
    if ( isfield(opt,'iFrame') ) 
        iFrame = opt.iFrame; 
    else iFrame = 1; 
    end;
    if ( isfield(opt,'nFrames') ) 
        nFrames = opt.nFrames; 
    else nFrames = 1; 
    end;
    if ( isfield(opt,'iBscan') ) 
        iFrame = opt.iBscan; 
    end;
    if ( isfield(opt,'metadata') )
        metadata = opt.metadata;
    else
        metadata = readMetaMgh(filename);
    end
    if ( isfield(opt,'singleFlag') )
        singleFlag = opt.singleFlag;
    else
        singleFlag = 0;
    end        
else
    error('wrong input variables')
end


mghType = metadata(2);
nAlinesPerFrame = metadata(3); 
nZpixels = metadata(4);
% nFrames = metadata(5);


switch mghType
    case 0
        datatype = 'uint8';
        nBytes = 1;
    case 1
        datatype = 'uint16';
        nBytes = 2;
    case 2
        datatype = 'float32';
        nBytes = 4;
    case 3
        datatype = 'float32';
        nBytes = 4;
        complexFlag = 1;
    otherwise
            error('error! readMGH: wrong file type');         
end
if singleFlag
    datatype = [datatype, '=>single'];
end

fileID = fopen(filename,'r');
% set offset in fseek
metalen_bytes = 1048576; %1024*1024
offset = metalen_bytes + (iFrame-1)*(1+complexFlag)*nAlinesPerFrame*nZpixels*nBytes;
% offset should be given in bytes. recall 1byte = 8bits
fseek(fileID,offset,'bof');
b = fread(fileID,(1+complexFlag)*nZpixels*nAlinesPerFrame*nFrames,datatype);

if complexFlag
    output = reshape(b(1:2:end),nZpixels,nAlinesPerFrame,nFrames) + ...
        sqrt(-1)*reshape(b(2:2:end),nZpixels,nAlinesPerFrame,nFrames);
else
    % ImageJ and MATLAB reads binary files in a different order, this is to
    % correnct for the flip and rotational effec.
    c = reshape(b,nAlinesPerFrame,nZpixels,nFrames); 
    c = flip(c,1);
    c = rot90(c,3);
    output = c;
end

fclose(fileID);


end