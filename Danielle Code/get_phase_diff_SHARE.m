%This code was developed to open pre-stored OCT phase data from a needle 
%probe and calculate the phase difference between subsequent A-lines. 
%Smoothing is applied. 

%V1.0, Alejandra Gómez-Ramírez and Danielle J. Harper, 2021
%Further information can be found in this paper:
   % https://doi.org/10.48550/arXiv.2305.14390
clc, close all, clear all

addpath ('E:\My Drive\4_UW-Madison\2_4_MGH+FDA+UW OCT\Needle_PS-OCT_on_WS\MGH_Needle_Doppler\Danielle Code') %code folder
Example_Data_files_folder = 'E:\My Drive\4_UW-Madison\2_4_MGH+FDA+UW OCT\needle_code_share\Example_data';
addpath (Example_Data_files_folder) % Files folder

%initial frame
readOpt.iFrame = 1;
%number of frames
readOpt.nFrames = 1024;

%open phase file. Only need one pol. state, but can try the other one to
%see if it looks better
phase1= readMgh('[p.needle_2][s.salmon][06-16-2021_14-15-13]phase1.mgh',readOpt);
%phase2= readMgh('[p.needle_2][s.salmon][06-16-2021_14-15-13]phase2.mgh',readOpt);

%back to the true phase values from -pi to pi (saved as 8 bit, 0-255 previously):
ph1=(single(phase1).*(2*pi)./255)-pi; %1024x1024x1024
%ph2=(single(phase2).*(2*pi)./255)-pi;


%number of frames
f=readOpt.nFrames;

%figure('Name','Calibration line'),plot(unwrap(ph1(473,:,a)));

%% phase difference of only odd columns
diff_ph1=zeros(1024,511,f); %1024x511x1024
for j=1:f
    for i=1:511
        diff=ph1(:,(2*(i-1))+1,j)-ph1(:,(2*(i-1))+3,j);
        diff_ph1(:,i,j)=diff;
    end
end

%phase difference of only the region of interest (sliced)
num_rows=272-220+1; %Number of rows

num_cols=215-205; %Calibration line (fiber tip)
%diffph1sl=diffph1(478:586,:,:); %slicing
diff_ph1_sliced=diff_ph1(220:272,:,:); %slicing

%concatenating 2 frames together
diff_ph_con_=zeros(num_rows,1022,f);
for l=1:2:f
    if l+1>f
        diff_ph_con_(:,1:511,l)=diff_ph1_sliced(:,:,l);
        break
    end
    diff_ph_con_(:,:,l)=horzcat(diff_ph1_sliced(:,:,l),diff_ph1_sliced(:,:,l+1));
end
diff_ph_con=diff_ph_con_(:,:,1:2:end); %because the for has a step of 2, even number frames are 0

%number of frames rounded
num_frm_con=ceil(f/2);

%Filter
i=0;j=0;k=0;z=0;
for k=1:num_frm_con
    for j=1:1022
        for i=1:num_rows
            v=diff_ph_con(i,j,k);
            if v>2.5 
                diff_ph_con(i,j,k)=0;
                z=z+1;
            end  
        end
    end
end
%%
%Moving median values every 73 columns
movm=zeros(num_rows,1022);
movma=zeros(num_rows,1058);
movmb=zeros(num_rows,1094);
diffphconM=zeros(num_rows,1022,num_frm_con);
diffphconMa=zeros(num_rows,1058,num_frm_con);
diffphconMb=zeros(num_rows,1094,num_frm_con);
for k=1:num_frm_con
    if num_frm_con==1
        movm=movmedian(diff_ph_con(:,:,k),73,2);
        diffphconM(:,:,k)=movm;
        break
    end
    if k+1>num_frm_con
        diffphconMa(:,:,k)=horzcat(diff_ph_con(:,987:1022,k-1),diff_ph_con(:,:,k));
        movma=movmedian(diffphconMa(:,:,k),73,2);
        diffphconM(:,:,k)=movma(:,37:end,:);
        break
    end
    if k~=1
        diffphconMb(:,:,k)=horzcat(diff_ph_con(:,987:1022,k-1),diff_ph_con(:,:,k),diff_ph_con(:,1:36,k+1));
        movmb=movmedian(diffphconMb(:,:,k),73,2);
        diffphconM(:,:,k)=movmb(:,37:1058,:);
        continue
    end
    diffphconMa(:,:,k)=horzcat(diff_ph_con(:,:,k),diff_ph_con(:,1:36,k+1));
    movma=movmedian(diffphconMa(:,:,k),73,2);
    diffphconM(:,:,k)=movma(:,1:1022,:);
    
 
end
 diffphconM1=diffphconM +(pi);
 figure('Name','movmedian+2pi'),imshow3D(diffphconM1,[0 2*pi]);
 colormap(redblue);
%%
%Average of each A-scan in region of interest
block=zeros(1,1022,num_frm_con);
block2=zeros(1,1022*num_frm_con);
for i=1:num_frm_con
    block(:,:,i)=mean(diffphconM(:,:,i),1);
    a=1022*(i-1)+1;
    b=1022*(i);
    block2(:,a:b)=block(:,:,i);
end
block2=block2-mean(block2(:,1:end));
figure('Name','Block'),plot(unwrap(block2(:,:)));

%Distance
w_length=0.0013; %mm %1.3um
n=1.39; %refractive index of tissue
dis=(block2.*w_length)./(4*pi*n);
%correct for drift with line below if needed. Region where needle is not touching tissue should be
%zero
dis=dis-1.6e-7;
time=(1:1022*num_frm_con)*40*(10^-6);
figure('Name','Relative distance'),plot(time,dis(:,:));
title('Relative Distance')
xlabel('time [s]') 
ylabel('Distance []')
sum=0;
totaldisr=zeros(1,1022*num_frm_con);
% for u=1:1022*fr
%     sum=sum+dis(1,u);
%     totaldisr(:,u)=sum;
% end
dis2 = movmedian(dis,5000);
for u=1:1022*num_frm_con
    sum=sum+dis2(1,u);
    totaldisr(:,u)=sum;
end
figure('Name','Absolute Distance'),plot(time,totaldisr(:,:));
hold on

title('Absolute Distance')
xlabel('time [s]') 
ylabel('Distance [mm]')

%save to .mat file
% Define the filename you want to save the variable under
filename = 'totaldisr_salmon.mat';

% Combine the folder and filename to create a full file path
fullFilePath = fullfile(Example_Data_files_folder, filename);

% Save the 'totaldisr' variable to the specified file path
save(fullFilePath, 'totaldisr');
%%

figure(11)
% Create axes
axes1 = axes('Parent',figure(11),...
    'Position',[0.623353293413173 0.167330677290833 0.449101796407186 0.697211155378486]);
hold(axes1,'on');
alpha 0.2

plot(time,totaldisr(:,:),'--o','Color',[0.77,0.72,0.97],'LineWidth',4,'MarkerSize',10,'MarkerFaceColor',[0.47,0.36,0.94],'MarkerEdgeColor',[0.47,0.36,0.94]);
%figure(11), legend({'Observed','Predicted'},'LineWidth',2);
figure(11), title('Doppler-tracked Absolute Distance','FontSize',40);
figure(11), xlabel('time [s]','FontSize',36);
figure(11), ylabel('Distance [mm]','FontSize',36);
xlim([0,25])
ylim([-5,55])
set(gca,'FontSize',36);
set(gcf, 'Units', 'normalized', 'outerposition', [0, 0, 1, 0.7],'Color',[1 1 1],'OuterPosition',[0 0 1 1])
set(get(gca,'YLabel'),'Rotation',0)
ylh = get(gca,'ylabel');
gyl = get(ylh);                
ylp = get(ylh, 'Position');
ylp(1) = -5;
set(ylh, 'Rotation',90, 'Position',ylp, 'VerticalAlignment','middle', 'HorizontalAlignment','center')
ax = gca;
num_cols = ax.Position;
ax.Position(1) = 0.1802;
pbaspect([1 1 1]);
set(axes1,'LineWidth',6,'TickLength',[0.025 0.04],'XTick',...
    [0 5 10 15 20 25],'FontSize',36);
%%
%Distance calibration line
w_length=0.0013; %mm %1.3um
n=1.3;
disc=(diffphconM(num_cols,:,:).*w_length)./(4*pi*n);
disc=disc-mean(disc(:,1:10000));
time=(1:1022*num_frm_con)*40*(10^-6);
figure('Name','Relative distance c'),plot(time,disc(:,:));
title('Relative Distance')
xlabel('time [s]') 
ylabel('Distance [mm]')
sum=0;
totaldisc=zeros(1,1022*num_frm_con);
for u=1:1022*num_frm_con
    sum=sum+disc(1,u);
    totaldisc(:,u)=sum;
end
figure('Name','Absolute Distance c'),plot(time,totaldisc(:,:));
title('Absolute Distance')
xlabel('time [s]') 
ylabel('Distance [mm]')


% %Tiff file
% filename = 'phasedifference.tiff';
% for k=1:fr
%     B=diffphconM1(:, :, k);
%     %figure,mesh(B)
%     B=uint8(B.*(255/(2*pi)));
%     B=ind2rgb(B,redblue);
%     if k==1
%         imwrite(B,filename);
%     else
%         imwrite(B,filename,'WriteMode','append');
%     end
%     
% end

