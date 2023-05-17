%% sliceview app:
% get settings parameters
clc
clear all
warning('off','all')
uuu=cd;
try
    disp('opening Slice Viewer..')
    transcendentalDoc = '';
    clc
    S = {...
        struct('name','Number of Slices','type','enum','values',{{'Sneak Peek','Color Image','Batch Images','Multi Folder','imRGB','scatter & alignment','2D Alignment','3D Alignment','2D Cell Scatter','3D Cell Scatter','AISH','FDISCO'}});...
        struct('name','Number of channels','type','int','default',1);...
        struct('name','Channel Name','type','enum','values',{{'DAPI','GFP','RFP','CY5'}});...
        struct('name','Suffix','type','str','default','.');...
        struct('name','Enhancment% (L 90)','type','int','default',99.99);...
        };
    
    Param = Settings_GUI(S);
    finder=string(cell2mat(Param(1)));
    n_idx=cell2mat(Param(2));
    channel=string(cell2mat(Param(3)));
    if channel=="DAPI"
        idx=1;
    elseif channel=="GFP"
        idx=2;
    elseif channel=="RFP"
        idx=3;
    else  %CY5
        idx=4;
    end
    suffix=char(Param(4));
    en=cell2mat(Param(5));
    %convert=cell2mat(Param(6));
    %edit=cell2mat(Param(7));
catch ME
    cd(uuu)
    %     disp(ME)
end
%

%%    get files
clc
warning('on','all')
if finder=="Multi Folder"
    fprintf('Please Select '); fprintf(2, 'Folder/s \n');
    listy=uipickfiles('num',[],'FilterSpec',['*',suffix,'tif'],'out','struct'); % loop images
    sorty={listy.name};
    sorty=string(natsortfiles(sorty'));
elseif finder=="2D Alignment"
    fprintf('Please Select '); fprintf(2, 'Folder/s \n');
    listy=uipickfiles('num',[],'FilterSpec',['*',suffix,'tif'],'out','struct'); % loop images
    trial=listy.name;
    slice_path=[trial,'/slices'];
elseif finder=="3D Alignment"
    fprintf('Please Select '); fprintf(2, 'Folder/s \n');
    listy=uipickfiles('num',[],'FilterSpec',['*',suffix,'tif'],'out','struct'); % loop images
    trial=listy.name;
    slice_path=[trial,'/slices'];
elseif finder=="3D Cell Scatter"
    fprintf('Please Select '); fprintf(2, 'Folder/s \n');
    listy=uipickfiles('num',[],'FilterSpec',['*',suffix,'tif'],'out','struct'); % loop images
    trial=listy.name;
    slice_path=[trial,'/slices'];
elseif finder=="FDISCO"
    fprintf('Please Select '); fprintf(2, 'Image/s \n');
    suffix='C0';
    listy=uipickfiles('num',[],'FilterSpec',[suffix,'*.tif'],'out','struct'); % loop images
    sorty={listy.name};
    sorty=string(natsortfiles(sorty'));
else
    fprintf('Please Select '); fprintf(2, 'Image/s \n');
    listy=uipickfiles('num',[],'FilterSpec',['*',suffix,'tif'],'out','struct'); % loop images
    sorty={listy.name};
    sorty=string(natsortfiles(sorty'));
end

clc

%% mono , folder, and batch to view imshow3D
if finder=="Batch Images"
    cd(listy(1).folder);
    C={};x=0;y=0;
    disp('reading and stiching images please wait...')
    for i=1:length(listy)   % loop images
        Ii={imadjust(imread((sorty(i)),'Index',idx))};
        x=max(x,size(cell2mat(Ii),1));y=max(y,size(cell2mat(Ii),2));
        C = cat(3,C,Ii);
        sprintf(['stiching... ', num2str(round(i*100/length(listy))),' %% ' ])
    end
    CC=[];
    
    %close all
    for i=1:length(listy)   % loop images
        Iz=zeros(x,y);
        Ii=imadjust(imread((sorty(i)),'Index',idx));
        D=size(Ii);
        Iz((1:D(1)),(1:D(2))) = Ii;
        % add zeros to sides
        CC = cat(3,CC,Iz);
        sprintf(['padarray to concatanate... ', num2str(round(i*100/length(listy))),' %% ' ])
    end
    figure;
    imshow3D(CC,[]);
    clc
    
elseif  finder=="Multi Folder" % run each folder
    C={};
    disp('reading and stiching images please wait...')
    x=0;
    y=0;
    for i=1:length(listy)
        
        I = dir([listy(i).name filesep ['*',suffix,'tif']]); % get the contents of the image_folder
        %all_files(i).names=I.name;
        cd(listy(i).name)
        Ii=imadjust(imread(I.name,'Index',idx));
        x=max(x,size(Ii,1));y=max(y,size(Ii,2));
        
        C = cat(3,C,Ii);
        sprintf(['stiching... ', num2str(round(i*100/length(listy))),' %% ' ])
    end
    
    CC=[];
    
    for i=1:length(listy)
        Iz=zeros(x,y);
        
        I = dir([listy(i).name filesep ['*',suffix,'tif']]); % get the contents of the image_folder
        %all_files(i).names=I.name;
        cd(listy(i).name)
        Ii=imadjust(imread(I.name,'Index',idx));
        D=size(Ii);
        Iz((1:D(1)),(1:D(2))) = Ii;
        CC = cat(3,CC,Iz);
        sprintf(['padarray to concatanate... ', num2str(round(i*100/length(listy))),' %% ' ])
    end
    
    %save('C.mat','C');
    figure;
    imshow3D(CC,[]);
    clc
elseif finder=="Sneak Peek"
    
    for i=1:length(listy)  % go over multidirectiores
        cd(listy(i).folder);
        disp('Montage image....')
        i
        im_tif = string(listy(i).name);
        diry=strsplit(im_tif,'/');
        diry=char(diry(end));
%         x1=imadjust(imread(diry,'index',1));
        x1=imread(diry,'index',1);
        if n_idx==1
            fig=figure('Name',num2str(diry),'NumberTitle','off');
            imshow(x1)
            imcontrast;
        elseif n_idx==2
            x2=imadjust(imread(diry,'index',2));
            fig=figure('Name',num2str(diry),'NumberTitle','off');  montage({x1,x2})
        elseif n_idx==3
            x2=imadjust(imread(diry,'index',2));
            x3=imadjust(imread(diry,'index',3));
            fig=figure('Name',num2str(diry),'NumberTitle','off'); montage({x1,x2,x3})
        else   % n_idx==4
            x4=imadjust(imread(diry,'index',4));
            x2=imadjust(imread(diry,'index',2));
            x3=imadjust(imread(diry,'index',3));
            fig=figure('Name',num2str(diry),'NumberTitle','off');
            montage([x1 x2;x3 x4])
            
        end
        
        %     title(num2str(diry))
    end
    disp('Finished')
%     autoArrangeFigures()
elseif finder=="imRGB" || finder=="2D Cell Scatter"
    for lis=1:length(listy)
        cd(listy(lis).folder);
        im_tif = sorty(lis);
        diry=strsplit(im_tif,'/');
        diry=char(diry(end));
        img=imread(diry,'index',1);
        img2=imread(diry,'index',2);
        img3=imread(diry,'index',3);
            % % % % disp('ch1:green |ch2:blue |ch3:red')
            % % % % pause(1)
            
%             img4=imread(diry,'index',4);
            
            % % % % rgbImage = cat(3, grayImageR, grayImageG, grayImageB);
            % % % % colorThresholder(rgbImage)
            disp('Loading multichannel image...');
            mainy=figure('Name',num2str(diry),'NumberTitle','off');mainy.Position=[500 0 1000 1500];
            tiledlayout(2,2)
            if  finder=="2D Cell Scatter"
                B=importdata([diry(1:(end-4)),'_position_counts_per_dapi.txt']);
                b_gfp=B.data(:,4)>0;
                b_rfp=B.data(:,5)>0;
                b_cy5=B.data(:,6)>0;
                
                hold on
                axis off
                axis equal
                
            end
            ax1 = nexttile;
            %     subplot(2,2,1)
            f1=imshow(img,[]);
            if  finder=="2D Cell Scatter"
                hold on
                scatter(B.data(:,1),B.data(:,2),100,'MarkerEdgeColor','g','LineWidth',2)
            end
            title("Ch1")
            ax2 = nexttile;
            %     subplot(2,2,2)
            f2=imshow(img2,[]);
            if  finder=="2D Cell Scatter"
                hold on
                scatter(B.data(b_gfp,1),B.data(b_gfp,2),50,'MarkerEdgeColor','cyan','LineWidth',2)
            end
            title("Ch2")
            ax3 = nexttile;
            %     subplot(2,2,3)
            f3=imshow(img3,[]);
            if  finder=="2D Cell Scatter"
                hold on
                scatter(B.data(b_rfp,1),B.data(b_rfp,2),50,'MarkerEdgeColor','magenta','LineWidth',2)
            end
            title("Ch3")
            ax4 = nexttile;
            %         subplot(2,2,4)
            f4=imshow(img4,[]);
            if  finder=="2D Cell Scatter"
                hold on
                scatter(B.data(b_cy5,1),B.data(b_cy5,2),50,'MarkerEdgeColor','y','LineWidth',2)
            end
            title("Ch4")
            impixelinfo;
            zoom on
            annotation(mainy,'textbox',...
                [0.212580086580087 0.914807303984314 0.571699117322441 0.0791075033179162],...
                'String','Right mouse click UP and DOWN to contrast | press SPACE to save contrast and Finish',...
                'LineStyle','none',...
                'FontSize',18);
            temp1=imcontrast(f1);temp1.Position = [100 500 300 300];
            temp2=imcontrast(f2);temp2.Position = [1600 500 300 300];
            temp3=imcontrast(f3);temp3.Position = [100 100 300 300];
            temp4=imcontrast(f4);temp4.Position = [1600 100 300 300];
            
            linkaxes([ax1 ax2 ax3 ax4],'xy')
            pause; % till space
            cmax1 = str2num(temp1.Children(1).Children(3).Children.Children(2).Children.Children(2).Children(2).String);
            cmin1 = str2num(temp1.Children(1).Children(3).Children.Children(2).Children.Children(2).Children(5).String);
            cont1=[cmin1 cmax1];
            cmax2 = str2num(temp2.Children(1).Children(3).Children.Children(2).Children.Children(2).Children(2).String);
            cmin2 = str2num(temp2.Children(1).Children(3).Children.Children(2).Children.Children(2).Children(5).String);
            cont2=[cmin2 cmax2];
            cmax3 = str2num(temp3.Children(1).Children(3).Children.Children(2).Children.Children(2).Children(2).String);
            cmin3 = str2num(temp3.Children(1).Children(3).Children.Children(2).Children.Children(2).Children(5).String);
            cont3=[cmin3 cmax3];
            cmax4 = str2num(temp4.Children(1).Children(3).Children.Children(2).Children.Children(2).Children(2).String);
            cmin4 = str2num(temp4.Children(1).Children(3).Children.Children(2).Children.Children(2).Children(5).String);
            cont4=[cmin4 cmax4];
            close all
            Jimg=imadjust(img,[cont1/65535]);Jimg2=imadjust(img2,[cont2/65535]);Jimg3=imadjust(img3,[cont3/65535]);
            Jimg4=imadjust(img4,[cont4/65535]);
            J(:,:,3) = Jimg2 + Jimg/3;
            J(:,:,1) = Jimg3 + Jimg/3;
            J(:,:,2) = Jimg4+ Jimg/3;
            f_f=figure;imshow(J);
            load('chaname.mat', 'chaname');
            title([['\color{gray}',num2str(chaname(1))],' ',['\color{blue}',num2str(chaname(2))],' ',['\color{red}',num2str(chaname(3))],' ',['\color{green}',num2str(chaname(4))]]);
            clc
            % % if  finder=="2D Cell Scatter"
            % %    hold on
            % %    axis off
            % %    axis equal
            % % scatter(B.data(:,1),B.data(:,2),50,'+','MarkerEdgeColor','g','LineWidth',2)
            % % scatter(B.data(b_gfp,1),B.data(b_gfp,2),50,'MarkerEdgeColor','cyan','LineWidth',2)
            % % scatter(B.data(b_rfp,1),B.data(b_rfp,2),50,'MarkerEdgeColor','magenta','LineWidth',2)
            % % scatter(B.data(b_cy5,1),B.data(b_cy5,2),50,'MarkerEdgeColor','y','LineWidth',2)
            % % end
            zoom(f_f)
            
            uiwait(f_f)
    end
            %%
elseif finder=="FDISCO"
    
    for lis=1:length(listy)
        cd(listy(lis).folder);
        im_tif = sorty(lis);
        diry=strsplit(im_tif,'/');
        diry=char(diry(end));
        diryx=regexprep(diry,'C0','C1');
        img=imread(diry,'index',1);
        img2=imread(diry,'index',2);
        img3=imread(diryx,'index',3);
            % % % disp('ch2:green |ch3:blue |ch4:red')
            % % % pause(1)
            % % % dapiy=imadjust(imread(filename,'index',1));
            % % % grayImageG=imread(filename,'index',4)+dapiy/3;
            % % % grayImageB=imread(filename,'index',2)+dapiy/3;
            % % % grayImageR=imread(filename,'index',3)+dapiy/3;
            % % % rgbImage = cat(3, grayImageR, grayImageG, grayImageB);
            % % % colorThresholder(rgbImage)
            disp('Loading multichannel image...');
            mainy=figure('Name',num2str(diry),'NumberTitle','off');mainy.Position=[0 500 2000 500];
                B=importdata([diry(1:(end-4)),'_position_counts_per_dapi.txt']);
                b_gfp=B.data(:,4)==1;
                b_rfp=B.data(:,5)==1;
                b_overlap=B.data(:,6)==1;
                hold on
                axis off
                axis equal
                zoom(mainy)
            
%             tiledlayout(1,3)
            
            ax1 =subplot(1,3,1);
            f1=imshow(img,[]);
                hold on
                scatter(B.data(b_overlap,1),B.data(b_overlap,2),100,'MarkerEdgeColor','y','LineWidth',2)
            title('Ch3')
            
            ax2 = subplot(1,3,2);
            f2=imshow(img2,[]);
            %%
                hold on
                scatter(B.data(b_gfp,1),B.data(b_gfp,2),100,'MarkerEdgeColor','cyan','LineWidth',2)
                %     scatter(B.data(b_overlap,1),B.data(b_overlap,2),100,'MarkerEdgeColor','y','LineWidth',2)
            
            title('Ch2')
            %    
            ax3 = subplot(1,3,3);
            f3=imshow(img3,[]);
            %%
                hold on
                scatter(B.data(b_rfp,1),B.data(b_rfp,2),100,'MarkerEdgeColor','magenta','LineWidth',2)
                %     scatter(B.data(b_overlap,1),B.data(b_overlap,2),100,'MarkerEdgeColor','y','LineWidth',2)
            
            title('Ch3')
            impixelinfo;
            annotation(mainy,'textbox',...
                [0.212580086580087 0.914807303984314 0.571699117322441 0.0791075033179162],...
                'String','Right mouse click UP and DOWN to contrast | press SPACE to save contrast and Finish',...
                'LineStyle','none',...
                'FontSize',18);
            temp1=imcontrast(f1);temp1.Position = [100 0 300 300];
            temp2=imcontrast(f2);temp2.Position = [800 0 300 300];
            temp3=imcontrast(f3);temp3.Position = [1500 0 300 300];
            linkaxes([ax1 ax2 ax3],'xy')
            pause;
            cmax1 = str2num(temp1.Children(1).Children(3).Children.Children(2).Children.Children(2).Children(2).String);
            cmin1 = str2num(temp1.Children(1).Children(3).Children.Children(2).Children.Children(2).Children(5).String);
            cont1=[cmin1 cmax1];
            cmax2 = str2num(temp2.Children(1).Children(3).Children.Children(2).Children.Children(2).Children(2).String);
            cmin2 = str2num(temp2.Children(1).Children(3).Children.Children(2).Children.Children(2).Children(5).String);
            cont2=[cmin2 cmax2];
            cmax3 = str2num(temp3.Children(1).Children(3).Children.Children(2).Children.Children(2).Children(2).String);
            cmin3 = str2num(temp3.Children(1).Children(3).Children.Children(2).Children.Children(2).Children(5).String);
            cont3=[cmin3 cmax3];
            close all
            Jimg=zeros(size(img)); Jimg2=imadjust(img2,[cont2/65535]);Jimg3=imadjust(img3,[cont3/65535]);
            J = cat(3, Jimg3,Jimg, Jimg2);
            
            f_f=figure;imshow(J);zoom(f_f)
            title('\color{blue}GFP, \color{red}RFP');
            %%
            clc
                hold on
                axis off
                axis equal
                scatter(B.data(b_gfp,1),B.data(b_gfp,2),100,'MarkerEdgeColor','b','LineWidth',2)
                scatter(B.data(b_rfp,1),B.data(b_rfp,2),100,'MarkerEdgeColor','r','LineWidth',2)
                scatter(B.data(b_overlap,1),B.data(b_overlap,2),100,'MarkerEdgeColor','y','LineWidth',2)
            
            uiwait(f_f)
    end  
    
    %%
elseif finder=="2D Alignment"
    allen_atlas_path = '/data/Aligment/AP_histology-master/AllenCCF';
    st = loadStructureTree([allen_atlas_path filesep 'structure_tree_safe_2017.csv']);
    AP_view_aligned_histology_all(st,trial);
    
    %%
elseif finder=="3D Alignment"
    disp('Loading libraries ...please wait...')
    allen_atlas_path = '/data/Aligment/AP_histology-master/AllenCCF';
    tv = readNPY([allen_atlas_path filesep 'template_volume_10um.npy']);
    av = readNPY([allen_atlas_path filesep 'annotation_volume_10um_by_index.npy']);
    st = loadStructureTree([allen_atlas_path filesep 'structure_tree_safe_2017.csv']);
    AP_view_aligned_histology_volume(tv,av,st,slice_path,2);
    %%
elseif finder=="3D Cell Scatter"
    braingrid;
    %%
elseif finder=="AISH"
    %          regions=input('ROI/s(,): ','s');
    %          regions='COAa,COAp,AAA,BA,CEAc,CEAl,CEAm,IA,MEAad,MEAav,MEApd,MEApv,LA,BLAa,BLAp,BLAv,BMAa,BMAp,PA,';
    %          regions=strsplit(regions,',');
    %          regions=regions(1:end-1);
    percx=0;
    cd(listy(1).folder);
    load('tonoraml.mat', 'tonormal')
    xnormal=normalize(tonormal);
    cmap = interp1([0;1],[1 1 1; 1 0 0],linspace(0,1,256));

    for i=1:length(listy)  % go over multidirectiores
        v = rescale(xnormal,1,256);%/0.000154462827933214);%slc17a6
        cd(listy(i).folder);
        im_tif = string(listy(i).name);
        diry=strsplit(im_tif,'/');
        diry=char(diry(end));
        roiextr=readtable([diry(1:(end-11)),'_ed1_sum_mono.csv']);
        regions=roiextr.ROI_names;
        perc=roiextr.perc;
        %             perc=perc*perc(end)/0.0000838577887568599;% ad
        load([diry(1:(end-8)),'_maskcrop.mat'],'zmask')
        masklined=boundarymask(zmask);
        allen_atlas_path = '/data/Aligment/AP_histology-master/AllenCCF';
        st = loadStructureTree([allen_atlas_path filesep 'structure_tree_safe_2017.csv']);
        
        if ~isempty(regions)
            allmask_R=zeros(size(zmask));
            allmask_G=zeros(size(zmask));
            allmask_B=zeros(size(zmask));
            
            %                         colors = distinguishable_colors(length(regions));
            %                         colors=redblue(length(regions));
            % create heatmap color scale
            %                       cmap = jet(256);
            
%             perc=xnormal((sum(percx)+1):length(perc)+(sum(percx)));
%             v = rescale(perc,1,256);%/0.000154462827933214);%slc17a6
            %                       v = rescale(perc,1,256*perc(end)/0.000154462827933214); %gad2
            %                       v = rescale(perc,1,256*perc(end)/0.000248372433059417); %slc17a7
            numValues = length(perc);
            markerColors = zeros(numValues, 3);
            % Now assign marker colors according to the value of the data.
            v=v((sum(percx)+1):length(perc)+(sum(percx)));% update for vector v
            for k = 1 : numValues
                row = round(v(k));
                markerColors(k, :) = cmap(row, :);
            end
            % end of scalecolor
            for k=1:length(regions) %  subregions
                roi_area=regions(k);
                roi_id_path = st.structure_id_path(find(strcmp(st.safe_name,roi_area)));%search full name
                if isempty(roi_id_path)
                    roi_id_path = st.structure_id_path(find(strcmp(st.acronym,roi_area)));% search acronym
                end
                roi_idx = find(contains(st.structure_id_path,roi_id_path));% all subregions
                xmask=ismember(zmask,roi_idx);
                mask_R= markerColors(k,1).*xmask;
                mask_G= markerColors(k,2).*xmask;
                mask_B= markerColors(k,3).*xmask;
                allmask_R=allmask_R+mask_R; % merge all binary masks
                allmask_G=allmask_G+mask_G; % merge all binary masks
                allmask_B=allmask_B+mask_B; % merge all binary masks
                
            end
            allmask=cat(3,allmask_R,allmask_G,allmask_B);
            %                     allmask(allmask==0)=255;
            %                     masklinedx=imfuse(masklined,allmask,'blend');
            [B] = bwboundaries(masklined);
            [nx,mx,~]=size(masklined);
            %                     masklinedy = imcrop(masklined,[round(mx/2) round(nx/2) mx nx]);
        end
        %                 masklined=imfuse(masklined,irgb,'blend');
        % save to pdf
        figure('Name',num2str(diry),'NumberTitle','off');
        set(gcf,'color','w')
        imshow(allmask)
        hold on
        for k = 1:length(B)
            boundaryx = B{k};
            plot(boundaryx(:,2), boundaryx(:,1), 'w', 'LineWidth', 1)
        end
        axis([round(mx/2) mx round(nx/2) nx])
        colormap(cmap)
        colorbar
        namepdf=['slc17a6_',num2str(i)];
        disp('Saving...')
        eval(['export_fig ',namepdf,'.pdf -r 600']);
        percx=numValues;
    end
    disp('Finished')
    %%
elseif finder=="scatter & alignment"
    regions=input('ROI/s(,): ','s');
    regions=strsplit(regions,',');
    regions=regions(1:end-1);
    for i=1:length(listy)  % go over multidirectiores
        cd(listy(i).folder);
        im_tif = string(listy(i).name);
        diry=strsplit(im_tif,'/');
        diry=char(diry(end));
        
        gimg=imadjust(imread(diry,'index',1));%,[100/65535 1000/65535]);
        try
            bimg=imadjust(imread(diry,'index',2));%,[2000/65535 7000/65535]);
            rimg=imadjust(imread(diry,'index',3));%,[2000/65535 7000/65535]);
            irgb = cat(3, rimg, gimg, bimg);
        catch
            irgb=gimg;
        end
        try
            load([diry(1:(end-8)),'_maskcrop.mat'], 'masklined')
            load([diry(1:(end-8)),'_maskcrop.mat'], 'zmask')
        catch
            load([diry(1:(end-8)),'_maskcrop.mat'], 'zmask')
            masklined = boundarymask(zmask);
            %                 se = strel('disk',3); %dilate lines of maskline
        end
        allen_atlas_path = '/data/Aligment/AP_histology-master/AllenCCF';
        st = loadStructureTree([allen_atlas_path filesep 'structure_tree_safe_2017.csv']);
        if ~isempty(regions)
            allmask=zeros(size(zmask));
            
            %          colors = distinguishable_colors(length(regions))
            for k=1:length(regions) %  subregions
                roi_area=regions(k);
                roi_id_path = st.structure_id_path(find(strcmp(st.safe_name,roi_area)));%search full name
                if isempty(roi_id_path)
                    roi_id_path = st.structure_id_path(find(strcmp(st.acronym,roi_area)));% search acronym
                end
                roi_idx = find(contains(st.structure_id_path,roi_id_path));% all subregions
                xmask=ismember(zmask,roi_idx);
                %                             allmask = cat(3,0.5*ones(size(zmask)).*xmask, 1*ones(size(zmask)).*xmask, 0.8*ones(size(zmask)).*xmask);
                allmask=[allmask | xmask]; % merge all binary masks
                
            end
            
            masklined=imfuse(masklined,allmask,'blend');
        end
        masklined=imfuse(masklined,irgb,'blend');
        figure('Name',num2str(diry),'NumberTitle','off');
        imshow(masklined)
        try
            B=importdata([diry(1:(end-4)),'_position_counts_per_dapi.txt']);
            hold on
            b_gfp=B.data(:,4)==1;
            b_rfp=B.data(:,5)==1;
            b_overlap=B.data(:,6)==1;
            scatter(B.data(b_overlap,1),B.data(b_overlap,2),50,'y');
            sg=scatter(B.data(b_gfp,1),B.data(b_gfp,2),50,'cyan');
            sr=scatter(B.data(b_rfp,1),B.data(b_rfp,2),50,'magenta');
            zvalues=zeros(length(b_gfp),1);
            roi_id_name=strings([length(zvalues), 1]);
            %                 splited=split(st.structure_id_path,"/");
            for zz=1:length(zvalues)
                zvalues(zz)=zmask(floor(B.data(zz,1)),floor(B.data(zz,2)));
                try
                    roi_id_name(zz) = st.safe_name(st.id==zvalues(zz));
                catch
                    roi_id_name(zz)=  zvalues(zz);
                end
            end
            annotationg = dataTipTextRow('Region',roi_id_name);
            sg.DataTipTemplate.DataTipRows(end+1) = annotationg;
            %                 zvalues=zeros(length(b_rfp),1);
            %                 roi_id_name=strings([length(zvalues), 1]);
            for zz=(1+sum(b_gfp)):length(zvalues)
                zvalues(zz)=zmask(floor(B.data(zz,1)),floor(B.data(zz,2)));
                try
                    roi_id_name(zz) = st.safe_name(st.id==zvalues(zz));
                catch
                    roi_id_name(zz)=  zvalues(zz);
                end
            end
            annotationr = dataTipTextRow('Region',roi_id_name);
            sr.DataTipTemplate.DataTipRows(end+1) = annotationr;
            drawnow % in order to make scatter on the last plot
        end
    end
    disp('Finished')
    %         autoArrangeFigures()
else     % color Image and stich
    disp('reading and coloring multipage TIFF image , this may take a while, please wait...')
    
    for i=1:length(listy)  % go over multidirectiores
        cd(listy(i).folder);
        im_tif = string(listy(i).name);
        diry=strsplit(im_tif,'/');
        diry=char(diry(end));
        
        if  n_idx ==1
            figure;   imshow(imadjust(imread(diry,'index',1)))
        elseif n_idx==2
            
        elseif  n_idx==3
            disp('DAPI |GFP | RFP');
            [mosaic] = TIFF_coloring2(listy(i).folder,diry,en);
            
        elseif  n_idx==4
            load('chaname.mat', 'chaname');
            disp([num2str(chaname(1)),' | ',num2str(chaname(2)),' | ',num2str(chaname(3)),' | ',num2str(chaname(4))]);
            [mosaic] = TIFF_coloring(listy(i).folder,diry,en);
        end
        C=[];
        clc
        disp('stiching image....')
        i
        for i=1:n_idx  % loop images
            Ii=mosaic(:,:,n_idx);
            C = cat(3,C,Ii);
        end
        figure;
        imshow3D(mosaic,[]);
        
    end
    disp('Finished')
    
end

