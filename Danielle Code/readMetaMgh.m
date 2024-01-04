function metadata = readMetaMgh(mghFilename)
% reads the metadata of .mgh file format.
%
% metadata = [0 mghType nAlinesPerFrame nZpixels nFrames]
% mghType: 0=uint8, 1=uint16, 2=float32, 3=float32(complex interleaved)
% refer to readMgh.m for more detail
%
% << snam@mit.edu 20150305 >>
%%%%%%%%%%%%%%
 

fileID = fopen(mghFilename,'r');
metadata = fread(fileID, 5, 'int32');
fclose(fileID);

end