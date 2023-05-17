% app for cropping multi slices 
close all
clear all
visy = input('visium(y)? ','s');% not empty y or n (and others)
if visy=='y'
cd('/data/runs/samples')  
suffix=''; starter=1;
else
suffix = input('Type suffix of file: ','s');
starter = input('From which file to start? ');
end
fprintf('Please Select '); fprintf(2, 'Folder/s \n');
listy=uipickfiles('num',[],'FilterSpec',['*',suffix,'.tif'],'out','struct'); % loop images 
%% run folders
idx=1;

for lis=1:length(listy)
  if visy=='y'
% image_file_names = dir([listy(lis).name filesep ['*',suffix,'_gray.tif']]); %  for dic and dapi 
image_file_names = dir([listy(lis).name filesep ['*gray.tif']]); % for dapi only 
idx=1;
  else % ~ visy 
image_file_names = dir([listy(lis).name filesep ['*',suffix,'.tif']]); % get the contents of the image_folder
  end
  image_file_names = natsortfiles({image_file_names.name});
  im_tif = string(image_file_names);
  [~, file_numCol] = size(image_file_names); %count number of image tif files in the directory
  im_path = fullfile(listy(lis).name, 'cropped'); % create multicrop file
if ~exist(im_path)
    mkdir(im_path)
end

for fol=starter:file_numCol
sprintf([num2str(fol),' image out of ',num2str(file_numCol),'\nfolder ',num2str(lis),' out of ',num2str(length(listy)), ' folder/s'])
trial=listy(lis).name;
cd(trial);
name_img=char(im_tif(fol));
disp('Extracting each channel values')
warning('off','all')
bfr= BioformatsImage(name_img);
warning('on','all')
%  I1= getPlane(bfr, 1, 1, 1);
%  disp("read ch1")
%  I2= getPlane(bfr, 1, 2, 1);
%  disp("read ch2")
%  I3= getPlane(bfr, 1, 3, 1);
%  disp("read ch3")
%  if bfr.sizeC==4
%  I4= getPlane(bfr, 1, 4, 1);
%  end

% num_sl = input('How many slices in image?');
I1=imread(name_img,'index',idx);
disp("reading multichannels... please wait...")
% I1=multi_img;
if visy~='y'
I2=imread(name_img,'index',2);
I3=imread(name_img,'index',3);
if bfr.sizeC==4 || bfr.sizeZ==4  || bfr.sizeT==4 
I4=imread(name_img,'index',4);
end
end

fig=figure('Name',name_img,'NumberTitle','off');
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
imshow(imadjust(I1),[]);
[h, w, ~] = size(I1);
sl=1;
imgzoompan(gcf, 'ImgWidth', w, 'ImgHeight', h);
while 1 % Trueuntil cancel
if visy=='y'
load([trial,'/mapped_xy.mat'], 'mapped_xy')
load([trial,'/10xMtxFile.mat'])
end
title(['Crop # ',num2str(sl)])
% prompt = {'Enter Slice Name (without format):'};
% dlgtitle = 'Name';
% dims = [1 35];
% definput = {name_img(1:(end-4))};
commandwindow;
answer=input('Enter Slice Name (without format):','s');
% answer = char(inputdlg(prompt,dlgtitle,dims,definput));
shg;
if isempty(answer)==1
       close all
       break;
else
answer=char(answer);
[~,rect] = imcrop;
cd(im_path);
x1=floor(rect(2));
xf=floor(rect(2)+rect(4));
y1=floor(rect(1));
yf=floor(rect(1)+rect(3));
if x1<=0 
    x1=1;
end
if y1<=0 
    y1=1;
end
title('Writing cropped image in directory...')
R1=I1(x1:xf,y1:yf);
if visy~='y'
imwrite(R1,[answer,'_cr','.tif']);
R2=I2(x1:xf,y1:yf);
imwrite(R2,[answer,'_cr','.tif'],'WriteMode','append');
R3=I3(x1:xf,y1:yf);
imwrite(R3,[answer,'_cr','.tif'],'WriteMode','append');
if bfr.sizeT==4 || bfr.sizeZ==4  || bfr.sizeC==4
R4=I4(x1:xf,y1:yf);
imwrite(R4,[answer,'_cr','.tif'],'WriteMode','append');
end


else % visy 
    
ix_path = fullfile([listy(lis).name,'/cropped'], num2str(answer)); % create multicrop file
if ~exist(ix_path)
    mkdir(ix_path)
end
cd(ix_path)
imwrite(R1,[answer,'.tif']);    
new_id=mapped_xy(:,1)>=x1 & mapped_xy(:,1)<=xf & mapped_xy(:,2)>=y1 & mapped_xy(:,2)<=yf;
% new_mapped_xy=[mapped_xy(new_id,2)-y1,mapped_xy(new_id,1)-x1];
% mapped_xy=mapped_xy(new_id,:);
mapped_xy=[mapped_xy(new_id,1)-x1,mapped_xy(new_id,2)-y1];
barcodesx=barcodesx(new_id,1);
datax=datax(:,new_id);
save('10xMtxFile.mat','barcodesx','datax','geneidx')
save('mapped_xy.mat', 'mapped_xy'); 
% save('new_id.mat','new_id'); 
% save('x1xfy1yf.mat', 'x1','xf','y1','yf');
end
end 


sl=sl+1;


end
end
% addpath(genpath(im_path)) 

      
end
disp('Finished')