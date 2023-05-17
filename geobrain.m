%% GeoBrain for Alignment to Brain coronal atlas
% geobrain3 : fdisco: /bigdata/microscope_images/presentation_demo/FDISCO/Tif
%% 1) set parameters
clearvars -except tv av st
warning('off','all')
set(0,'DefaultFigureWindowStyle','normal')
uuu=cd;
% cd('/data/Aligment/AP_histology-master/AllenCCF'); % for table view
% pause(0.5)
% cd('/bigdata')
try
    transcendentalDoc = '';
    clc
    S = {...
        struct('name','Geo Step','type','enum','values',{{'All','Resize','Adjust','Rotate','Flip & Order','Register','Manual','Extract'}},'doc',transcendentalDoc);...
        struct('name','File Suffix','type','str','default','ed.');...
        struct('name','Rotation','type','int','default',0);...
        struct('name','Experiment','type','enum','values',{{'4X','20X','AISH','VIS','FDISCO'}});...
        struct('name','Species','type','enum','values',{{'Mouse','Bat','ZebraFish'}});...
        struct('name','Hemisphere','type','enum','values',{{'Both','Right','Left'}});...
        struct('name','Transform(F|R)','type','enum','values',{{'pwl','affine'}});...
        };
    % 4x==>
    
    cont=0;
    vis=0;
    ish=0;
    
    
    % %
    Param = Settings_GUI(S);
    stepy=string(Param(1));
    step_org=stepy;
    suffix=char(Param(2));
    rotit=cell2mat(Param(3));
    expy=cell2mat(Param(4));
    if expy=="20X"
        cont=1;
    elseif expy=="VIS"
        vis=1;
    elseif expy=="AISH"
        ish=1;
    end
    species=string(Param(5));
    if species=="Bat"
        baty=1;gf=0;
    elseif species=="ZebraFish"
        baty=1;gf=1;
    else % mouse
        baty=0;gf=0;
    end
    hemi=string(Param(6));
    if hemi=="Right"
        hemi='l';
    elseif hemi=="Left"
        hemi='r';
    else % both hemi
        hemi=[];
    end
    trnsfrm=char(Param(7));
    if vis==1 % visium file
        cd('/data/runs/samples/')
        fprintf('Please Select '); fprintf(2, 'Visium Mother Folder/s \n');
        listy=uipickfiles('num',[],'FilterSpec',[],'out','struct'); % loop images
        
    else % not visium
        fprintf('Please Select '); fprintf(2, 'Folder/s \n');
        listy=uipickfiles('num',[],'FilterSpec',['*',suffix,'tif'],'out','struct'); % loop images
        cd(listy(1).name)
    end
catch
    close all
    cd(uuu);
    return ;
end

close all
clc
warning('on','all')
%starter1 = input('From which FOLDER to start? ');
%starter2 = input('From which Image to start? ');%%%%%%%%%%%%5
%rot = input('What is your prefered rotation? ');
if cont==1
    suffix2 = input('Type suffix of subimage file: ','s');
    hemi=[];
end
if stepy=="Extract"
    start_slice=input('from which slice to extract: ');
else
    start_slice=1;
end % auto=1; % to count auto not manual draw
%% 1) Load CCF
if baty==0
    if exist('tv','var') || exist('av','var') || exist('st','var')
        disp('Loaded Libraries')
    else
        
        disp('Loading Libraries...')
        % Load CCF atlas
        allen_atlas_path = '/data/Aligment/AP_histology-master/AllenCCF';
        %addpath(genpath('/data/Technion_analysis/Amygdala/FISH/Aligment'))
        tv = readNPY([allen_atlas_path filesep 'template_volume_10um.npy']);
        av = readNPY([allen_atlas_path filesep 'annotation_volume_10um_by_index.npy']);
        st = loadStructureTree([allen_atlas_path filesep 'structure_tree_safe_2017.csv']);
    end
elseif gf==0 % bat
    disp('Loading Libraries...')
    %     cd('/data/Aligment/Bat mask/Batamir/')
    tv=  cell2mat(struct2cell(load('/data/Aligment/Bat mask/Batamir/coronal_labels.mat', 'coronal_labels'))); % annotated
    av= cell2mat(struct2cell(load('/data/Aligment/Bat mask/Batamir/coronal_labels.mat', 'coronal_labels'))); % annotated
    st = readtable('/data/Aligment/Bat mask/Batamir/Atlas_labels.csv');% regions names
else % baty==1 and gf
    disp('Loading Libraries...')
    %     cd('/data/Aligment/Bat mask/Batamir/')
    %     tv=  cell2mat(struct2cell(load('/data/Aligment/GoldFish Mask/zebrafish_atlas.mat', 'gfa'))); % annotated
    tv=  cell2mat(struct2cell(load('/data/Aligment/goldfish atlas/stacked.mat', 'im_pad'))); % annotated
    %     av= cell2mat(struct2cell(load('/data/Aligment/GoldFish Mask/zebrafish_atlas.mat', 'gfa'))); % annotated
    av=  cell2mat(struct2cell(load('/data/Aligment/goldfish atlas/stacked.mat', 'im_pad'))); % annotated
    st = readtable('/data/Aligment/GoldFish Mask/zebrafish_atlas.csv');% regions names
end

%% set paths for slide and slice images tif & nd2
for lis=1:length(listy)  % go over multidirectiores
    
    
    trial=listy(lis).name;
    rot=rotit;
    cd(trial)
    
    if expy=="FDISCO"
        image_file_names = dir([trial filesep [suffix,'*.tif']]); % get the contents of the image_folder
        image_file_names = natsortfiles({image_file_names.name});
        im_tif = string(image_file_names);
        
    elseif ish==0
        image_file_names = dir([trial filesep ['*',suffix,'tif']]); % get the contents of the image_folder
        image_file_names = natsortfiles({image_file_names.name});
        im_tif = string(image_file_names);
    else
        image_file_names = dir([trial filesep ['*',suffix,'jpg']]); % get the contents of the image_folder
        image_file_names = natsortfiles({image_file_names.name});
        im_tif = string(image_file_names);
        
    end
    
    % keep record file
    number_of_files=length(im_tif)
    info_geobrain=[hemi,' , ',num2str(number_of_files),' , ', date];
    save('info_geobrain.mat','info_geobrain') % save info about geobrain
    % % % % im_path = fullfile(trial, 'processed');
    % % % % if ~exist(im_path)
    % % % %     mkdir(im_path)
    % % % % end
    slice_path = fullfile(trial, 'slices');
    if ~exist(slice_path)
        mkdir(slice_path)
    end
    addpath(genpath(slice_path))
    [~, file_numCol] = size(image_file_names); %count number of image tif files in the directory
    %im_path = 'path to folder with images goes here';
    %slice_path = [im_path filesep 'slices'];
    %     savefig(gcf,'volir.fig') % take property name
    
    %% padarray to same size
    cd(trial)
    if vis~=1
        if stepy=="Resize"||stepy=="All"
            padquit=0;
            if cont==0
                hi=0;we=0;
                for pady=1:file_numCol
                    ss = imfinfo(im_tif(pady));
                    we=max(we,ss(1).Width);
                    hi=max(hi,ss(1).Height);
                    if stepy=="All" && pady==1 % chcek only the first time and only when in all mode
                        sslast=imfinfo(im_tif(file_numCol));lastw=sslast.Width;
                        ssmid=imfinfo(im_tif(round(file_numCol/2)));midw=ssmid.Width;
                        if lastw==midw && lastw==we
                            padquit=1;
                            break;
                        end
                    end
                    clear ss
                    %         end
                end
                %   if vis~=1
                if padquit==0 % do pads from start
                    pb = CmdLineProgressBar('Padding images with zeros to stack... ');
                    
                    for pady=1:file_numCol
                        if ish ==1
                            pb.print(pady,file_numCol)
                            namex=char(im_tif(pady));
                            pad_im=im2gray(imread(namex));[mp,np]=size(pad_im);
                            K_pad = padarray(pad_im, [floor((hi-mp)/2) floor((we-np)/2)], 0,'post');
                            K_pad = padarray(K_pad, [ceil((hi-mp)/2) ceil((we-np)/2)], 0,'pre');
                            imwrite(K_pad,namex);
                            pad_im2=im2gray(imread([namex(1:(end-4)),'1.jpg']));[mp,np]=size(pad_im2);
                            K_pad2 = padarray(pad_im2, [floor((hi-mp)/2) floor((we-np)/2)], 0,'post');
                            K_pad2 = padarray(K_pad2, [ceil((hi-mp)/2) ceil((we-np)/2)], 0,'pre');
                            imwrite(K_pad2,[namex(1:(end-4)),'1.jpg']);
                        else
                            pb.print(pady,file_numCol)
                            pad_im=imread(im_tif(pady),'index',1);
                            [mp,np]=size(pad_im);
                            pad_im2=imread(im_tif(pady),'index',2);
                            try
                            pad_im3=imread(im_tif(pady),'index',3);
                            end
                            K_pad = padarray(pad_im, [floor((hi-mp)/2) floor((we-np)/2)], 0,'post');
                            K_pad = padarray(K_pad, [ceil((hi-mp)/2) ceil((we-np)/2)], 0,'pre');
                            imwrite(K_pad,char(im_tif(pady)));
                            K2_pad = padarray(pad_im2, [floor((hi-mp)/2) floor((we-np)/2)], 0,'post');
                            K2_pad = padarray(K2_pad, [ceil((hi-mp)/2) ceil((we-np)/2)], 0,'pre');
                            imwrite(K2_pad,char(im_tif(pady)),'WriteMode','append');
                            try % only if there is 3rd channel 
                                 K3_pad = padarray(pad_im3, [floor((hi-mp)/2) floor((we-np)/2)], 0,'post');
                                K3_pad = padarray(K3_pad, [ceil((hi-mp)/2) ceil((we-np)/2)], 0,'pre');
                                imwrite(K3_pad,char(im_tif(pady)),'WriteMode','append');   
                            end
                        end
                    end
                else  % dont do pads
                    %  do nothing
                end
                
                
                
                
                
                
                %   end
            end
            
            if stepy~="All"
                inputy = input('continue to Adjust(y)?:','s');
                if inputy=="y"
                    stepy="Adjust";
                else
                    clc
                    return
                end
            end
        end
    end
    %% check if num.tif in slices needs to be resized and minorized
    cd(trial)
    if vis~=1
        if stepy~="All"&& stepy~="Adjust" && stepy~="Resize"
            if cont==0
                disp('downscaling images by 0.1')
                slice_info=imfinfo(im_tif(1));% real
                cd(slice_path)
                rgb_info=imfinfo('1.tif');
                w1=slice_info.Width;w2=rgb_info.Width;
                h1=slice_info.Height;
                if w1<=w2
                    pb = CmdLineProgressBar('Resizing... ');
                    for rsz=1: file_numCol
                        pb.print(rsz,file_numCol)
                        curr_s=imresize(imread([num2str(rsz),'.tif']),0.1);%%%%%%%%%%%%%% resize HERE !!!
                        imwrite(curr_s, [num2str(rsz),'.tif']);
                    end
                end
            end
        end
    end
    % %%%%%%%%%%%%%%%%% if cleary==1
    % % disp('clearing previous variables...'); pause(0.5);
    % % clearing
    % %
    % % end
    %% edit: rotate and crop and mask
    if stepy=="Adjust"|| stepy=="All"
        %% 2) Preprocess slide images to produce slice images
        cd(slice_path)
        slice_images=1; % true if in each image there is one slice
        if cont==1 && trnsfrm=="affine"|| vis==1 && trnsfrm=="affine" || expy=="FDISCO" % 20 x one section at  a time
            resize_factor=1;
        else % 4x stacked alignment
            disp('downscaling resolution and resizing image')
            resize_factor=0.1;
        end
        AP_process_histology(trial,resize_factor,slice_images,suffix,vis,expy);
        if stepy~="All"&& vis~=1
            inputy = input('continue to Rotate(y)?:','s');
            if inputy=="y"
                stepy="Rotate";
            else
                clc
                return
            end
        end
    end
    %%% save flip rot order thre in slices...
    cd(slice_path)
    %%%
    %% rotate 2D plane
    if stepy=="All" || stepy=="Rotate"
        if vis~=1
            ang_rotate=zeros(1,file_numCol);
            %             try
            %                 load([slice_path,'orders.mat']','all_order')
            %             catch
            all_order=1:file_numCol;
            %             end
            %             try
            %                 load([slice_path,'flipers.mat'],'all_flip')
            %             catch
            all_flip=uint8(zeros(file_numCol,1));
            %             end
            if stepy~="All"
                while 1
                    starter = input('From which image to rotate(num/all=1000/enter to quit)?:');
                    if isempty(starter)
                        break
                    else
                        AP_rotate_histology(slice_path,starter,0,all_flip,all_order,expy);
                    end
                end
            else
                AP_rotate_histology(slice_path,1000,0,all_flip,all_order,expy);
            end
            % (optional) Rotate, center, pad, flip slice images
            if stepy~="All"
                inputy = input('continue to Register(y)?:','s');
                if inputy=="y"
                    stepy="Register";
                else
                    clc
                    return
                end
            end
        else
            J=imread('1.tif');
            ftilt=figure;
            hIm = imshow(J);
            camroll(rot)
            try
                tilter;%% go and tilt
                uiwait(ftilt)
                CurrAng=round(CurrAng);
                rot=rot+CurrAng;
            end
            
            
            J=imrotate(J,rot,'crop');
            imwrite(J,'1.tif');
            
            save('ang_rotate.mat','rot')
        end
    end
    %% flipper and orders
    if stepy=="Flip & Order"
        %             try
        %                 load([slice_path,'orders.mat']','all_order')
        %             catch
        all_order=1:file_numCol;
        %             end
        %             try
        %                 load([slice_path,'flipers.mat'],'all_flip')
        %             catch
        all_flip=uint8(zeros(file_numCol,1));
        %             end
        AP_rotate_histology(slice_path,1000,1,all_flip,all_order,expy);
        
        % (optional) Rotate, center, pad, flip slice images
        if stepy~="All"
            inputy = input('continue to Register(y)?:','s');
            if inputy=="y"
                stepy="Register";
            else
                clc
                return
            end
        end
        
    end
    %% 3) Align angle CCF to slices
    
    % Find CCF slices corresponding to each histology slice
    if stepy=="All" || stepy=="Register"
        disp("Loading images for registration, please wait...")
        AP_grab_histology_ccf(tv,av,st,slice_path,baty,gf,expy);
        if stepy=="Register"
            warny = input('continue to Auto[PREVIOUS MASKS WILL BE DELETED](y)?:','s');
        end
        if stepy=="All" || warny=="y"
            % manual (auto is now a decoy)
            disp("Making auto alingnment decoy mask")
            if baty~=1
                auto_array= [3.25161696225555,0.00213521415871705,0;0.000963895600246671,3.27120525079391,0;1779.84929313089,1448.69801886149,1];
            else % baty==1
                auto_array= [3.25161696225555,0.00213521415871705,0;0.000963895600246671,3.27120525079391,0;1779.84929313089,1448.69801886149,1];
            end
            cd(slice_path)
            for decoy=1:file_numCol
                
                atlas2histology_tform{decoy, 1}=auto_array;
            end
            save('atlas2histology_tform.mat','atlas2histology_tform')
        end
        % Align CCF slices and histology slices
        % (first: automatically, by outline)
        % % % disp('Autoaligning.... this may take a while....please wait...')
        % % % AP_auto_align_histology_ccf(slice_path);
        if stepy~="All"
            inputy = input('continue to Manual(y)?:','s');
            if inputy=="y"
                stepy="Manual";
            else
                clc
                return
            end
        end
    end
    
    %% upscale images again
    % %     cd(trial)
    % %     if trnsfrm=="xxxxx"
    % %     if vis~=1
    % %         if cont==0
    % %             disp('upscaling images')
    % %             slice_info=imfinfo(im_tif(1));% real
    % %             cd(slice_path)
    % %             rgb_info=imfinfo('1.tif');
    % %             w1=slice_info.Width;w2=rgb_info.Width;
    % %             h1=slice_info.Height;
    % %             if w1>w2
    % %                 pb = CmdLineProgressBar('Resizing... ');
    % %                 for rsz=1: file_numCol
    % %                     pb.print(rsz,file_numCol)
    % %                     curr_s=imresize(imread([num2str(rsz),'.tif']),[h1,w1]);
    % %                     imwrite(curr_s, [num2str(rsz),'.tif']);
    % %                 end
    % %             end
    % %         else
    % %                     w1=1; h1=1;
    % %         end
    % %     else
    % %         w1=1; h1=1;
    % %     end
    % %     else
    % %      w1=1; h1=1;
    % %     end
    %% manual align outline
    warning('off','all')
    if stepy=="All" || stepy=="Manual"
        w1=1; h1=1;
        
        AP_manual_align_histology_ccf(tv,av,st,slice_path,cont,w1,h1,vis,baty,gf,trnsfrm,im_tif,expy);
        
        
        stepy="Extract";
    end
    %%
    if expy=="FDISCO"
        tvx=permute(tv,[2 3 1]);
        disp('Aligning whole image with downsampled image')
        twin_file_names = dir([slice_path filesep '*tif']); % get the contents of the image_folder
        twin_file_names = natsortfiles({twin_file_names.name});
        l_dir=length(twin_file_names);
        slice_im_dirx=1:50:l_dir;
        if slice_im_dirx(end)~=l_dir
            slice_im_dirx=  [slice_im_dirx,l_dir];
        end
        twin_file_names=twin_file_names(slice_im_dirx);
        atlas_idx=zeros(length(slice_im_dirx),2);
        atlas_idx(:,1)=slice_im_dirx';
        im_twin = string(twin_file_names);
        all_sainity=zeros(length(im_twin),1);
        try
            load([slice_path,'/hist_pts.mat'], 'hist_pts')
            load([slice_path,'/atlas_pts.mat'], 'atlas_pts')
        end
        [~, f_length] = size(twin_file_names); %count number of image tif files in the directory
        if vis~=1 % NOT visium
            load([slice_path,'/ang_rotate.mat'], 'ang_rotate');% also put it in visium
            load([slice_path,'/orders.mat'], 'all_order')
            try
                load([slice_path,'/flipers.mat'], 'all_flip')
                %  all_flip=all_flip(all_order);
            end
            
            try % updtae the real images names
                im_tif = string(image_file_names);
                im_tifn=im_tif;
                im_tif=im_tif(all_order);
                ang_rotate=ang_rotate(all_order);
            end
        end
        for curr_slice =start_slice:f_length
            try
                f_name=char(im_tifn(slice_im_dirx(curr_slice)));f_name=f_name(1:(end-4));f_name
            catch
                f_name=char(im_tif(curr_slice));f_name=f_name(1:(end-4));f_name
            end
            % % % % %      cd(slice_path)
            %     f_name=char(im_twin(curr_slice));f_name=f_name(1:(end-4));f_name
            if vis==1
                ang=rot;
            else
                ang=rot+ang_rotate(curr_slice);% rotate to match 4x
            end
            cd(trial)
            if ish==0
                orgin=imread(im_tif(curr_slice),'index',1);
            else
                orgin=im2gray(imread(im_tif(curr_slice)));
            end
            tilted=imrotate(orgin,ang,'crop');
            try
                if all_flip(curr_slice,1)==1 % check if to flip
                    tilted=flip(tilted ,2);
                end
            end
            %write first page of updated rotated oredered fliped  edited
            cd(slice_path)
            
            
            %%crop the 20x position from 4x
            % Loop through images, get fluorescence in ROI
            [me,ne]=size(tilted);
            ccf_slice_fn = [slice_path filesep 'histology_ccf.mat'];
            load(ccf_slice_fn);
            % Warp AV slice to align
            curr_slice_av_unaligned = histology_ccf(curr_slice).av_slices;
            curr_slice_av_unaligned(isnan(curr_slice_av_unaligned)) = 1;
            atlas_idx(curr_slice,2)=round(mean(mean(histology_ccf(curr_slice).plane_ap)));
        end
        %%
        ntv_idx=zeros(l_dir,2);
        ntv_idx(:,1)=1:l_dir;
        for curr_slice=1:l_dir
            app_idx=find(curr_slice./slice_im_dirx>=1);
            rem_d=rem(curr_slice,50)/50; % get the remainder of divison then get its ratio fron 50 ;
            if app_idx(end)==length(slice_im_dirx)
                app_sl=round((rem_d*atlas_idx(app_idx(end),2))+((1-rem_d)*atlas_idx(app_idx(end),2)));% 0.02 =1 / 50
            else
                app_sl=round((rem_d*atlas_idx(app_idx(end),2))+((1-rem_d)*atlas_idx(app_idx(end)+1,2)));% 0.02 =1 / 50
            end
            allxmask(:,:,curr_slice)=tvx(:,:,app_sl);
            ntv_idx(curr_slice,2)= app_sl;
        end
    end % fdisco
    %     sliceViewer(allxmask)
    %% align with 20x
    if stepy=="All" || stepy=="Extract"


        
        % % % % % % % % % % % % % % %
        if cont==1
            %
            disp('Aligning whole image with 20X magnified subimage')
            
            twin_file_names = dir([trial filesep ['*',suffix2,'.tif']]); % get the contents of the image_folder
            % twin_file_names = natsortfiles({twin_file_names.name});
            im_twin = twin_file_names.name
            
            % [~, f_length] = size(twin_file_names); %count number of image tif files in the directory
            load('ang_rotate.mat', 'ang_rotate');
            try
                load([slice_path,'/hist_pts.mat'], 'hist_pts')
                load([slice_path,'/atlas_pts.mat'], 'atlas_pts')
            end
            %         for x=1:f_length
            %             x
            cd(slice_path)
            % f_name=char(im_tif(x));
            % % big=im2gray(imread('1.tif'));
            % % big = imrotate(imread(f_name),ang,'crop');
            big=im2gray(imread('1.tif'));
            %      ch1=char(im_twin(x));
            %      [~,name1,~] = fileparts(ch1);
            %cimg= [num2str(name),'.tif'];
            %cimg=TIFF_coloring(cimg,Amp);
            ang=rot+ang_rotate(1);% rotate to match 4x
            cd(trial)
            %imwrite(cimg, [num2str(name),'.tif']);
            orgin=imread(twin_file_names.name,'index',1); orgin2=imread(twin_file_names.name,'index',2);
            orgin3=imread(twin_file_names.name,'index',3);
            orgin4=imread(twin_file_names.name,'index',4);
            %         ch2=char(im_twin(x));
            %      [~,name2,~] = fileparts(ch2);
            tilted=imrotate(orgin,ang,'crop');
            tilted2=imrotate(orgin2,ang,'crop');
            tilted3=imrotate(orgin3,ang,'crop');
            tilted4=imrotate(orgin4,ang,'crop');
            % % %      img=imresize(tilted,(1/(1.63/0.32)));%0.02 or 0.005?
            % img=uint8(round(double(img)/(2^8)));% convert to uint8
            %Read two images into the workspace, and convert them to grayscale for use with normxcorr2. Display the images side-by-side.
            %small =rgb2gray(img);
            
            %montage({big,small})
            % % % % % % % % figure;
            % % % % % % % % montage({imadjust(big), imadjust(img)},[]);
            % % % % % % % % title('Montage of subimage and whole image')
            % % % % % % % %
            % % % % % % % % %Perform cross-correlation, and display the result as a surface.
            % % % % % % % % disp('Conducting Normalized 2-D cross-correlation...');
            % % % % % % % % c = normxcorr2(img,big);
            % % % % % % % %
            % % % % % % % % %surf(c)
            % % % % % % % % %shading flat
            % % % % % % % %
            % % % % % % % % %Find the peak in cross-correlation.
            % % % % % % % % [ypeak,xpeak] = find(c==max(c(:)));
            % % % % % % % %
            % % % % % % % % %Account for the padding that normxcorr2 adds.
            % % % % % % % % yoffSet = ypeak-size(img,1);
            % % % % % % % % xoffSet = xpeak-size(img,2);
            % % % % % % % %
            % % % % % % % % %Display the matched area by using the drawrectangle function. The 'Position' name-value pair argument specifies the upper left coordinate, width, and height of the ROI as the 4-element vector [xmin,ymin,width,height]. Specify the face of the ROI as fully transparent.
            % % % % % % % %
            % % % % % % % % imshow(big,[]);
            % % % % % % % % hold on;
            % % % % % % % % % % rectangle('Position',[xoffSet,yoffSet,size(img,2),size(img,1)],...
            % % % % % % % % % %   'EdgeColor', 'r',...
            % % % % % % % % % %   'LineWidth', 3,...
            % % % % % % % % % %   'LineStyle','-')
            % % % % % % % % roiy=drawrectangle(gca,'Position',[xoffSet,yoffSet,size(img,2),size(img,1)], ...
            % % % % % % % %     'FaceAlpha',0);
            % % % % % % % % title('In blue - subimage position in main image')
            
            % % % % pause(2);
            % % % % close all
            %         end
            % crop the 20x position from 4x
            % Loop through images, get fluorescence in ROI
            ccf_slice_fn = [slice_path filesep 'histology_ccf.mat'];
            load(ccf_slice_fn);
            % Load histology/CCF alignment
            ccf_alignment_fn = [slice_path filesep 'atlas2histology_tform.mat'];
            load(ccf_alignment_fn);
            
            slice_im_dir = dir([slice_path filesep '*.tif']);
            slice_im_fn = natsortfiles(cellfun(@(path,fn) [path filesep fn], ...
                {slice_im_dir.folder},{slice_im_dir.name},'uni',false));
            
            slice_fluor = nan(length(slice_im_fn),1);
            
            
            for curr_slice = start_slice:length(slice_im_fn)
                
                try
                    load([slice_path,'/','cpmpfp.mat']);
                    disp('Loaded prevoius points')
                    [mp,fp] = cpselect(big,imadjust(tilted),mp,fp,'Wait',true);
                catch
                    [mp,fp] = cpselect(big,imadjust(tilted),'Wait',true);
                end
                
                % Load slice image
                curr_slice_im = imread(slice_im_fn{curr_slice});
                
                % Warp AV slice to align
                curr_slice_av_unaligned = histology_ccf(curr_slice).av_slices;
                curr_slice_av_unaligned(isnan(curr_slice_av_unaligned)) = 1;
                if trnsfrm=="affine"
                    tform = affine2d;
                    tform.T = atlas2histology_tform{curr_slice};
                    tform_size = imref2d([size(curr_slice_im,1),size(curr_slice_im,2)]);
                    curr_slice_av = imwarp(curr_slice_av_unaligned, ...
                        tform,'nearest','OutputView',tform_size);
                else % pwl trnsfrm
                    resize_factor = 0.1;
                    atlas_xy=atlas_pts{curr_slice,1};
                    hist_xy=hist_pts{curr_slice,1}*(1/resize_factor);
                    tform_size = imref2d([size(tilted,1),size(tilted,2)]);
                    disp("Transforming orignal...")
                    tform = fitgeotrans(atlas_xy,hist_xy,trnsfrm);
                    curr_slice_av= imwarp(histology_ccf(curr_slice).av_slices,tform,'nearest','OutputView',tform_size);
                    curr_slice_av(isnan(curr_slice_av)) = 1;
                end
                % % % yendSet=(yoffSet+size(img,1));
                % % % xendSet=(xoffSet+size(img,2));
                % % % [maxy,maxx]=size(curr_slice_av);
                % % %
                % % %     if xoffSet<=0
                % % %     xoffSet=1;
                % % %     end
                % % % if yoffSet<=0
                % % %    yoffSet=1;
                % % % end
                % % %  if yendSet>maxy
                % % %     yendSet=maxy;
                % % %  end
                % % %   if xendSet>maxx
                % % %     xendSet=maxx;
                % % %   end
                % % %
                % % % cropmat=curr_slice_av(yoffSet:yendSet,xoffSet:xendSet);
                % % zmask=imresize(cropmat,[size(tilted,1),size(tilted,2)],'nearest');
                % % zmask(zmask<2)=1;
                %%%
                % xmask=zmask>2;
                % set(gcf,'name','Please select in PARALLEL a minimum of 4 CONTOUR points','numbertitle','off')
                disp('Transforming 20x...')
                t = fitgeotrans(mp,fp,'nonreflectivesimilarity');
                Rfixed = imref2d(size(tilted));
                zmask = imwarp(curr_slice_av,t,'nearest','OutputView',Rfixed);
                
                disp('Fusing mask and subimage..')
                % % imshowpair(tilted,zmask,'blend');
                
                %%%
                tilted(zmask==1)=0;
                tilted2(zmask==1)=0;
                tilted3(zmask==1)=0;
                tilted4(zmask==1)=0;
                disp("writing image channels")
                % %  CC=imfuse(zmask,tilted,'blend');
                % % figure; imshow(CC,[])
                cd(trial);
                % %  save([im_twin(1:(end-4)),'_maskcrop.mat'],'zmask');
                save([slice_path,'/','cpmpfp.mat'],'mp','fp');
                disp('Saved points');
                imwrite(tilted,[im_twin(1:(end-4)),'_geo.tif']);
                disp('Done Channel 1');
                imwrite(tilted2,[im_twin(1:(end-4)),'_geo.tif'],'WriteMode','append');
                disp('Done Channel 2');
                imwrite(tilted3,[im_twin(1:(end-4)),'_geo.tif'],'WriteMode','append');
                disp('Done Channel 3');
                imwrite(tilted4,[im_twin(1:(end-4)),'_geo.tif'],'WriteMode','append');
                disp('Done Channel 4');
                %  make plot mask lines and save them;
                maskline = boundarymask(zmask);
                if baty~=1
                    colrmask;
                    %  make plot mask lines and save them;
                    disp('Saving Mask');
                    save([im_twin(1:(end-4)),'_maskcrop.mat'],'zmask','maskline','RGBmask');
                else % baty
                    save([im_twin(1:(end-4)),'_maskcrop.mat'],'zmask','maskline');
                end
            end % cuur_slice
            
        else  % align with 4x/AISH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% cont=0
            %%
            disp('Aligning whole image with downsampled image')
            twin_file_names = dir([slice_path filesep '*tif']); % get the contents of the image_folder
            twin_file_names = natsortfiles({twin_file_names.name});
            im_twin = string(twin_file_names);
            all_sainity=zeros(length(im_twin),1);
            try
                load([slice_path,'/hist_pts.mat'], 'hist_pts')
                load([slice_path,'/atlas_pts.mat'], 'atlas_pts')
            end
            [~, f_length] = size(twin_file_names); %count number of image tif files in the directory
            if vis~=1 % NOT visium
                load([slice_path,'/ang_rotate.mat'], 'ang_rotate');% also put it in visium
                load([slice_path,'/orders.mat'], 'all_order')
                try
                    load([slice_path,'/flipers.mat'], 'all_flip')
                    %  all_flip=all_flip(all_order);
                end
                
                try % updtae the real images names
                    im_tif = string(image_file_names);
                    im_tifn=im_tif;
                    im_tif=im_tif(all_order);
                    ang_rotate=ang_rotate(all_order);
                end
            end
            for curr_slice =start_slice:f_length
                try
                    f_name=char(im_tifn(curr_slice));f_name=f_name(1:(end-4));f_name
                catch
                    f_name=char(im_tif(curr_slice));f_name=f_name(1:(end-4));f_name
                end
                % % % % %      cd(slice_path)
                %     f_name=char(im_twin(curr_slice));f_name=f_name(1:(end-4));f_name
                if vis==1
                    ang=rot;
                else
                    ang=rot+ang_rotate(curr_slice);% rotate to match 4x
                end
                cd(trial)
                if ish==0
                    orgin=imread(im_tif(curr_slice),'index',1);
                else
                    orgin=im2gray(imread(im_tif(curr_slice)));
                end
                tilted=imrotate(orgin,ang,'crop');
                try
                    if all_flip(curr_slice,1)==1 % check if to flip
                        tilted=flip(tilted ,2);
                    end
                end
                %write first page of updated rotated oredered fliped  edited
                cd(slice_path)
                
                
                %%crop the 20x position from 4x
                % Loop through images, get fluorescence in ROI
                [me,ne]=size(tilted);
                ccf_slice_fn = [slice_path filesep 'histology_ccf.mat'];
                load(ccf_slice_fn);
                % Warp AV slice to align
                curr_slice_av_unaligned = histology_ccf(curr_slice).av_slices;
                curr_slice_av_unaligned(isnan(curr_slice_av_unaligned)) = 1;
                if trnsfrm=="xxxx"%affine - out of loop
                    % Load histology/CCF alignment
                    ccf_alignment_fn = [slice_path filesep 'atlas2histology_tform.mat'];
                    load(ccf_alignment_fn);
                    
                    slice_im_dir = dir([slice_path filesep '*.tif']);
                    slice_im_fn = natsortfiles(cellfun(@(path,fn) [path filesep fn], ...
                        {slice_im_dir.folder},{slice_im_dir.name},'uni',false));
                    
                    slice_fluor = nan(length(slice_im_fn),1);
                    disp("extracting other channels indexes values")
                    
                    
                    % Load slice image
                    curr_slice_im = imread(slice_im_fn{curr_slice});
                    
                    
                    
                    tform = affine2d;
                    tform.T = atlas2histology_tform{curr_slice};
                    tform_size = imref2d([size(curr_slice_im,1),size(curr_slice_im,2)]);
                    zmask = imwarp(curr_slice_av_unaligned, ...
                        tform,'nearest','OutputView',tform_size);
                else % "pwl"
                    resize_factor = 0.1;
                    atlas_xy=atlas_pts{curr_slice,1};
                    hist_xy=hist_pts{curr_slice,1}*(1/resize_factor);
                    tform_size = imref2d([me,ne]);
                    disp("Transforming orignal...")
                    try % check if transform is possible
                        tform = fitgeotrans(atlas_xy,hist_xy,trnsfrm);
                        zmask= imwarp(histology_ccf(curr_slice).av_slices,tform,'nearest','OutputView',tform_size);
                        zmask(isnan(zmask)) = 1;
                        sainity=1
                    catch
                       
                        sainity=0
                        zmask=histology_ccf(curr_slice).av_slices;
                    end
                end
                %                 zmask=curr_slice_av;%(yoffSet:size(curr_slice_av,1),xoffSet:size(curr_slice_av,2));
                if sainity~=-1
                    disp("writing image channels")
                    [mz,nz]=size(zmask);
                    
%                     ref_rect=zeros(mz,nz);
ref_rect=ones(mz,nz);
                    % zmask=round(zmask);
                    
                    if vis~=1
                        zmask(zmask<2)=1;
                        %                     CC=imfuse(zmask,tilted,'blend');figure;imshow(CC,[])
                        
                        tilted(zmask==1)=0;
                        % %                     try
                        % %                         if all_flip(curr_slice,1)==1 % check if to flip
                        % %                             tilted=flip(tilted ,2);
                        % %                         end
                        % %                     end
                        ref_rect(zmask>=2)=1;stats=regionprops(ref_rect,'BoundingBox');
                        % figure;imshow(imadjust(ref_rect),[]);impixelinfo;
                        thisBB = stats.BoundingBox;
                        BB1=round(thisBB(1));BB2=round(thisBB(2));
                        BB13=round(thisBB(1)+thisBB(3));BB24=round(thisBB(2)+thisBB(4));
                        if thisBB(1)<0
                            thisBB(1)=1;
                        end
                        if thisBB(3)<0
                            thisBB(3)=1;
                        end
                        if BB13>ne
                            BB13=ne;
                        end
                        if BB24>me
                            BB24=me;
                        end
                        cd(trial)
                        
                        %wether its one channel or 2 or 3
                        
                        if ish==0
                            try
                            orgin2=imread(im_tif(curr_slice),'index',2);
                            catch 
                            orgin2=zeros(size(orgin));
                            end
                            try
                            orgin3=imread(im_tif(curr_slice),'index',3);
                            tilted3=imrotate(orgin3,ang,'crop');
                            end
                        else % AISH
                            orgin2=im2gray(imread([f_name,'1.jpg']));% expression
                            %                     orgin3=im2gray(imread([f_name,'2.jpg']));% ISH (1 is Nissl)
                        end
                        
                        tilted2=imrotate(orgin2,ang,'crop');
                      
                        try
                            if all_flip(curr_slice,1)==1 % check if to flip
                                tilted2=flip(tilted2 ,2);
                                try
                                tilted3=flip(tilted3 ,2);
                                end
                            end
                        end
                        
                        tilted2(zmask==1)=0;
                        tilted2=tilted2(BB2:BB24,BB1:BB13);
                        if ish==0
                            try
                            tilted3(zmask==1)=0;
                            tilted3=tilted3(BB2:BB24,BB1:BB13);
                            end
                        end
                        zmask=zmask(BB2:BB24,BB1:BB13);
                        tilted=tilted(BB2:BB24,BB1:BB13);
                        
                        % do it again for index2 and index 3
                        
                        
                        
                        
                        %  for etay cut   mask intpo HEMISPHERE
                        avg_LinPos=round(size(zmask,2)/2);
                        
                        if hemi=='l' %avg_LinPos<(size(zmask,1)/2) % delete what is below
                            zmask(:,(1:avg_LinPos))=0;
                            tilted(:,(1:avg_LinPos))=0;
                            if ish==0
                                tilted2(:,(1:avg_LinPos))=0;
                                try
                                tilted3(:,(1:avg_LinPos))=0;
                                end
                            end
                            %      maskline(:,(1:avg_LinPos))=0;
                        elseif hemi=='r'  %hemi='r'% delete what is above
                            zmask(:,(avg_LinPos:size(zmask,2)))=0;
                            tilted(:,(avg_LinPos:size(zmask,2)))=0;
                            if ish==0
                                tilted2(:,(avg_LinPos:size(zmask,2)))=0;
                                try
                                tilted3(:,(avg_LinPos:size(zmask,2)))=0;
                                end
                            end
                            %      maskline(:,(max_LinPos:size(zmask,1)))=0;
                        end
                        if ~isempty(hemi)
                            [me,ne]=size(tilted);[mz,nz]=size(zmask);
                            ref_rect=zeros(mz,nz);
                            ref_rect(zmask>=2)=1;stats=regionprops(ref_rect,'BoundingBox');
                            % figure;imshow(imadjust(ref_rect),[]);impixelinfo;
                            thisBB = stats.BoundingBox;
                            BB1=round(thisBB(1));BB2=round(thisBB(2));
                            BB13=round(thisBB(1)+thisBB(3));BB24=round(thisBB(2)+thisBB(4));
                            if thisBB(1)<0
                                thisBB(1)=1;
                            end
                            if thisBB(3)<0
                                thisBB(3)=1;
                            end
                            if BB13>ne
                                BB13=ne;
                            end
                            if BB24>me
                                BB24=me;
                            end
                            cd(trial)
                            
                            
                            zmask=zmask(BB2:BB24,BB1:BB13);
                            tilted=tilted(BB2:BB24,BB1:BB13);
                            % do it again for index2 and index 3
                            tilted2=tilted2(BB2:BB24,BB1:BB13);tilted2(zmask==1)=0;
                            if ish==0
                                try
                                tilted3=tilted3(BB2:BB24,BB1:BB13);tilted3(zmask==1)=0;
                                end
                            end
                        end
                        % % xmask=zmask>2;
                        % % % set(gcf,'name','Please select in PARALLEL a minimum of 4 CONTOUR points','numbertitle','off')
                        % %
                        % % [mp,fp] = cpselect(xmask,imadjust(tilted),'Wait',true);
                        % % t = fitgeotrans(mp,fp,'polynomial');
                        % % Rfixed = imref2d(size(tilted));
                        % % zmask = imwarp(zmask,t,'OutputView',Rfixed);
                        % % disp('Fusing mask and subimage..')
                        % % imshowpair(tilted,zmask,'blend');
                    else % visium
                        cd(trial)
                        zmask(zmask<2)=0;
                        fmidliner=figure;
                        imshow(imadjust(tilted),[]);
                        title('select one midpoint')
                        [linx,liny]=ginputcr(1,'Color', 'r', 'LineWidth', 3);
                        %                     hold on
                        %                     scatter(linx,liny)
                        %                     midliner;% draw midline
                        close all
                        %                     LinPos=round(LinPos);
                        %                     max_LinPos=round(max(LinPos(1,1),LinPos(2,1)));
                        %                     min_LinPos=round(min(LinPos(1,1),LinPos(2,1)));
                        avg_LinPos=[round(linx),round(liny)];%round((LinPos(1,1)+LinPos(2,1))/2);
                        if avg_LinPos<(size(zmask,1)/2) % delete what is below
                            zmask(:,(1:avg_LinPos))=0;
                            tilted(:,(1:avg_LinPos))=0;
                            %      maskline(:,(1:avg_LinPos))=0;
                        else % delete what is above
                            zmask(:,(avg_LinPos:size(zmask,1)))=0;
                            tilted(:,(avg_LinPos:size(zmask,1)))=0;
                            %      maskline(:,(max_LinPos:size(zmask,1)))=0;
                        end
                        % visium!!
                        tilted=imrotate(tilted,-rot,'crop'); % -rot to bring it back to orignial
                        zmask=imrotate(zmask,-rot,'crop'); % -rot to bring it back to orignial
                        %      maskline = boundarymask(zmask);% create mask of roi lines
                        % % %
                        % % % zmask(zmask<2)=1;
                        % % % zmask(zmask==1)=0;
                        % % % ref_rect=zmask;
                        % % % zmask( all(ref_rect==0,2), : ) = [];  %d rows
                        % % % zmask( :, all(ref_rect==0,1) ) = [];  % d columns
                        % % % % zmask( all(ref_rect2==0,2), : ) = [];  %d rows
                        % % % tilted(all(ref_rect==0,2), : )=[];
                        % % % tilted( :, all(ref_rect==0,1) ) = [];  % d columns
                        % % % %  tilted(all(ref_rect2==0,2), : )=[];
                        % % % zmask(zmask==0)=1;
                        % % % tilted(zmask==1)=0;
                        % % % ref_rect2=tilted;
                        % % % zmask( :, all(ref_rect2==0,1) ) = [];  % d columns
                        % % % tilted( :, all(ref_rect2==0,1) ) = [];  % d columns
                        
                    end
                    % C=imfuse(zmask,tilted,'blend');imshow(C);imshow(tilted,[]);imshow(zmask,[]);
                    
                    % % % % % % f_name=char(im_tif(curr_slice));f_name=f_name(1:(end-4));
                    % % rectangle('Position', [thisBB(1),thisBB(2),thisBB(3),thisBB(4)],...
                    % % 'EdgeColor','r','LineWidth',2 );
                    % % % % % % % % % % % % % % % % % % % % % % % % % % % figure;imshowpair(maskline,imadjust(tilted),'blend');
                    % save([num2str(f_name),'_maskcrop.mat'],'zmask');
                    % cd(trial);
                    
                    if vis~=1
                        imwrite(tilted,[num2str(f_name),'_geo.tif']);
                        disp('Done Channel 1');
                        if ish==0
                            imwrite(tilted2,[num2str(f_name),'_geo.tif'],'WriteMode','append');
                            disp('Done Channel 2');
                            try
                            imwrite(tilted3,[num2str(f_name),'_geo.tif'],'WriteMode','append');
                            disp('Done Channel 3');
                            end
                        else % AISH
                            imwrite(tilted2,[num2str(f_name),'1_geo.tif']);
                            disp('Done Channel 2');
                            %                     imwrite(tilted3,[num2str(f_name),'2_geo.tif']);
                            %                     disp('Done Channel 3');
                        end
                    end
                    maskline = boundarymask(zmask);
                    se = strel('disk',1); %dilate lines of maskline
                    if baty~=1
                        colrmask;% rgbmask output ||||||||||| THIS IS  A FUNCTION!!!!!!!!!!!!
                        %  make plot mask lines and save them;
                        if vis==1 % delete extra areas
                            %                     maskline=imrotate(maskline,rot,'crop');
                            col_nums = find(~any(maskline,1))';   %default for 2D array
                            row_nums = find(~any(maskline,2));
                            gaps=diff(row_nums)~=1;
                            fgap=logical(gaps==1);
                            gapr1=row_nums(fgap);% find where are the zeros
                            %                     gapr2=row_nums(fgap+1);
                            gaps=diff(col_nums)~=1;
                            fgap=logical(gaps==1);
                            gapc1=col_nums(fgap);
                            %                     gapc2=col_nums(fgap+1);
                            maskline(row_nums,:)=[];
                            maskline(:,col_nums)=[];
                            RGBmask(row_nums,:,:)=[];
                            RGBmask(:,col_nums,:)=[];
                            masklined = ~imdilate(maskline,se);
                            RGBmask=imfuse(masklined,RGBmask,'blend');
                            binaryImage=~masklined;
                            [B,~]=bwboundaries(binaryImage);
                            save([num2str(f_name),'_maskcrop.mat'],'zmask','masklined','RGBmask','gapc1','gapr1','B');
                        else % regular non vis
                            masklined = ~imdilate(maskline,se);
                            RGBmask=imfuse(masklined,RGBmask,'blend');
                            save([num2str(f_name),'_maskcrop.mat'],'zmask','masklined','RGBmask');
                            %                     imshow(RGBmask)
                        end
                    else % baty
                        masklined = ~imdilate(maskline,se);
                        save([num2str(f_name),'_maskcrop.mat'],'zmask','masklined');
                    end % baty
                    disp('Saved Mask');
                    close all
                end
                all_sainity(curr_slice)=sainity;
            end  %curr slice
            save([slice_path,'/all_sainity.mat'],'all_sainity')
        end
    end
    close all
    stepy=step_org;
    
    %% reserve for quick modifications
    % % % %     if stepy=="XXX"
    % % % %         cd(trial)
    % % % %         chtif=char(im_tif)
    % % % %     load([chtif(1:(end-4)),'_maskcrop.mat'],'zmask','masklined','RGBmask','gapc1','gapr1');
    % % % %     end
end
disp('Finished ;)');
warning('on','all')

