%This code was developed to use pre-stored phase differences measured from 
%a needle probe and create surface-referenced maps of intensity, phase 
%retardation and optic axis as the needle moves.  

%V1.0, Danielle J. Harper and Benjamin J. Vakoc, 2023
%Further information can be found in this paper:
   % https://doi.org/10.48550/arXiv.2305.14390

clear all;

%open intensity file
ii_file = '[p.needle_2][s.salmon][06-16-2021_14-15-13]ii.mgh';

%open averaged oa, rr and ii images 
oa_avg = imread('OAavg_salmon.tif');
rr_avg = imread('rravg_salmon.tif');
ii_avg = imread('iiavg_salmon.tif');

%open distance from Doppler
load('totaldisr_salmon.mat')
z = totaldisr(1:511:end)*600;
pixel = round(z/max(z)*1399)+1;
pixel(pixel<1)=1;
%needle visualization
k = [-10 0; 10 0; 10 1425; -10 1385] ; 

fid = fopen(ii_file,'rb');
%fidrr = fopen(rr_file,'rb');
fseek(fid,1024*1024,'bof');

ang = 88:0.05:92;
v = [];
clear rr_Int;
oa_Int = uint8(255*(rand(3600,200,3)*0.5+0.5));
rr_Int = uint8(255*(rand(3600,200)*0.5+0.5));
ii_Int = uint8(255*(zeros(3600,200)));

height = zeros(3600,1,1);
writemode = 'overwrite';

 writerObj = VideoWriter('salmon_rr.avi','Uncompressed AVI');
 writerObj.FrameRate=50;
 open(writerObj);
 
% Need to change from the default renderer to zbuffer to get it to work right.
% openGL doesn't work and Painters is way too slow.
set(gcf, 'renderer', 'zbuffer');


for i=1:700
    
    fseek(fid,1024*1024 + 1024*1024*(i-1+200),'bof');
    img = fread(fid,[1024 1024],'uint8=>single').';
    %imgrr = fread(fidrr,[1024 1024],'uint8=>single').';
    
    for j=1:8
        
        simg = img(500:600,(1:128) + (j-1)*128);
        R = radon(edge(simg),ang);
        [~,mi] = max(std(R));
        v = [v ang(mi)];
        
    end
   % temp = cumsum(tand((v-90))*128);
    temp = pixel(i+200);
    pos(i) = round(temp(end));
 %    pos(i) = z(i)+1;
    %rr(i) = mean(mean(imgrr(500:600,:)));
%     figure(1);
%     plot(temp);
%     drawnow;
    pos = (pos + 1401);
    
    figure(2);
    height(pos(i):(pos(i)+554)) = height(pos(i):(pos(i)+554)) + 1;
    oa_Int(pos(i):(pos(i)+554),2:end,:) = oa_Int(pos(i):(pos(i)+554),1:(end-1),:);
    oa_Int(pos(i):(pos(i)+554),1,:) = oa_avg(470:1024,i+200,:);
    
    rr_Int(pos(i):(pos(i)+554),2:end,:) = rr_Int(pos(i):(pos(i)+554),1:(end-1));
    rr_Int(pos(i):(pos(i)+554),1,:) = rr_avg(470:1024,i+200);
    
    ii_Int(pos(i):(pos(i)+554),2:end,:) = ii_Int(pos(i):(pos(i)+554),1:(end-1));
    ii_Int(pos(i):(pos(i)+554),1,:) = ii_avg(470:1024,i+200);
    
    
    oa_sInt = oa_Int((pos(i)-1400):(pos(i)+400),:,:);
    oa_sInt(1400:1410,:,:) = 0;
    
    rr_sInt = rr_Int((pos(i)-1400):(pos(i)+400),:);
    rr_sInt(1400:1410,:) = 0;
    
    ii_sInt = ii_Int((pos(i)-1400):(pos(i)+400),:);
    ii_sInt(1400:1410,:) = 0;
    
    height_s = height((pos(i)-1400):(pos(i)+400));
    height_s(height_s>199)=199;
    %if (mod(i,10)==0)
        oa_sIntB = zeros(length(height_s),200,3,'uint8');
        for j=1:length(height_s)
            if (height_s(j)>1)
                oa_sIntB(j,1:200,1) = uint8(interp1(1:height_s(j),single(oa_sInt(j,1:height_s(j),1)),linspace(1,sqrt(height_s(j)),200).^2));
                oa_sIntB(j,1:200,2) = uint8(interp1(1:height_s(j),single(oa_sInt(j,1:height_s(j),2)),linspace(1,sqrt(height_s(j)),200).^2));
                oa_sIntB(j,1:200,3) = uint8(interp1(1:height_s(j),single(oa_sInt(j,1:height_s(j),3)),linspace(1,sqrt(height_s(j)),200).^2));
            else
                oa_sIntB(j,1:200,1:3) = uint8(0);
            end
        end
        
        rr_sIntB = zeros(length(height_s),200,'uint8');
        for j=1:length(height_s)
            if (height_s(j)>1)
                rr_sIntB(j,1:200,1) = uint8(interp1(1:height_s(j),single(rr_sInt(j,1:height_s(j),1)),linspace(1,sqrt(height_s(j)),200).^2));
               
            else
                rr_sIntB(j,1:200) = uint8(0);
                
            end
        end
        
                ii_sIntB = zeros(length(height_s),200,'uint8');
        for j=1:length(height_s)
            if (height_s(j)>1)
                ii_sIntB(j,1:200,1) = uint8(interp1(1:height_s(j),single(ii_sInt(j,1:height_s(j),1)),linspace(1,sqrt(height_s(j)),200).^2));
               
            else
                ii_sIntB(j,1:200) = uint8(0);
                
            end
        end
        
     cmap=colormap(bone(256));
     cmap = cmap(:,[2 1 3]);
    %choose which plot to make here - ii, rr or oa
       % imshow(imrotate(imfuse(flip(ii_sIntB,2),ii_sIntB,'Montage'),0));
        imshow(imrotate(imfuse(flip(rr_sIntB,2),rr_sIntB,'Montage'),0),cmap);
       % imshow(imrotate(imfuse(flip(oa_sIntB,2),oa_sIntB,'Montage'),0));
        %imwrite(imrotate(imfuse(flip(oa_sIntB,2),oa_sIntB,'Montage'),90),'out.tif','WriteMode',writemode);writemode = 'append';
        
        patch([0 0 400 400],[1400 1800 1800 1400],'y','FaceAlpha',0.1)
        patch(200+k(:,1),k(:,2),[0.5,0.5,0.5])
       %   patch(200+k(:,1),pixel(i+200)-117+k(:,2),'c')
        
        drawnow;
        i
        
        thisFrame = getframe(gca);
	% Write this frame out to a new video file.
  	writeVideo(writerObj, thisFrame);
	%myMovie(frameIndex) = thisFrame;
    %end
    
    
    
end
close(writerObj)
fclose(fid);
        
