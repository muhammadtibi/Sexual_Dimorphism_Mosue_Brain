%% spot: for counting dots
% try: cd('/bigdata/microscope_images/Lin/Analysis/72_01')
clear all
warning('off','all')
uuu=cd;
% cd('/bigdata')
% etay rabies :cd example: "/bigdata/microscope_images/Etay/Amygdala/Brains/Females/EA3_2/EA3_2_MEA_Tiff/3_2_L"
% preview = input('Preview Counter(y)? ','s');
listx = {'Rabies','smFISH','Rabies-3ch','FDISCO','Rabies-No DAPI'};
[exp,~] = listdlg('ListString',listx,'SelectionMode','single','Name','App Selection','ListSize',[150,85]);
% prlow = input('add percent [1:2]: ');
% try
transcendentalDoc = '';
if exp==2
    exp="smFISH";
    S = {...
        struct('name','File Suffix','type','str','default','geo.');...
        struct('name','Enhancment%[234] (L 90)','type','str','default','99.99,99.99,99.99');...
        struct('name','P-Max-Transform[234] ','type','str','default','8000,8000,8000');...
        struct('name','Manual Threshold','type','int','default',3);...
        struct('name','Magnification','type','int','default',20);...
        struct('name','Diameter(μ): High ','type','int','default',30);...
        struct('name','                 Low','type','int','default',10);...
        struct('name','Figure Output?','type','checkbox');...
        % struct('name','Preview?','type','checkbox');...
        };
    Param = Settings_GUI(S);
    suffix=char(Param(1));
    en=Param(2);
    str = regexprep(en,',',' ');
    en = str2num(cell2mat(str));en1=en(1);en2=en(2);en3=en(3);
    pmt=Param(3);
    str0 = regexprep(pmt,',',' ');
    pmt = str2num(cell2mat(str0));pmt1=pmt(1);pmt2=pmt(2);pmt3=pmt(3);
    manual_th=cell2mat(Param(4));
    mag=cell2mat(Param(5));
    if mag==10
        mag=0.65;
    elseif mag==20
        mag=0.32;
    elseif mag==40
        mag=0.16;
    else  %4
        mag=1.63;
    end
    d_high=cell2mat(Param(6));
    d_low=cell2mat(Param(7));
    size_high=round(pi*((d_high/2)^2)./mag);
    size_low=round(pi*((d_low/2)^2)./mag);
    fig_flag=cell2mat(Param(8));
    if fig_flag==1
        fig_fig='on';
    else
        fig_fig='off';
    end
    %      preview=cell2mat(Param(9));
elseif exp==3 %%% AISH
   exp="Rabies-3ch";
    S = {...
        struct('name','File Suffix','type','str','default','.');...
        struct('name','*factor GFP|BFP|RFP','type','str','default','2.4,2.4,2.4');...
        struct('name','Magnification','type','int','default',4);...
        struct('name','Diameter(μ): High ','type','int','default',15);...
        struct('name','              Low','type','int','default',6);...
        struct('name','Overlap(pixels)','type','int','default',15);...
        struct('name','Frame','type','int','default',5);...
        struct('name','out bound dist','type','int','default',100);...
        struct('name','Figure Output?','type','checkbox','default',1);...
        };

    Param = Settings_GUI(S);
    suffix=char(Param(1));
    prlow=Param(2);
    str = regexprep(prlow,',',' ');
    prlow = str2num(cell2mat(str));en1=prlow(1);en2=prlow(2);en3=prlow(3);
    % % bin=Param(3);
    % % str1 = regexprep(bin,',',' ');
    % % bin = str2num(cell2mat(str1));bin1=bin(1);bin2=bin(2);
    mag=cell2mat(Param(3));
    if mag==10
        mag=0.65;
    elseif mag==20
        mag=0.32;
    elseif mag==40
        mag=0.16;
    else  %4
        mag=1.63;
    end
    d_high=cell2mat(Param(4));
    d_low=cell2mat(Param(5));
    size_high=round(pi*((d_high/2)^2)./mag);
    size_low=round(pi*((d_low/2)^2)./mag);
    d_over=cell2mat(Param(6));
    frame=cell2mat(Param(7));
    disty=cell2mat(Param(8));
    fig_flag=cell2mat(Param(9));
    if fig_flag==1
        fig_fig='on';
    else
        fig_fig='off';
    end
elseif exp==4
    exp="FDISCO";


    S = {...
        struct('name','File Suffix','type','str','default','C0');...
        struct('name','Adaptive th. (GFP|RFP)','type','str','default','0.62,0.62');...
        struct('name','Magnification','type','int','default',4);...
        struct('name','Diameter(μ): High ','type','int','default',15);...
        struct('name','              Low','type','int','default',6);...
        struct('name','Overlap(pixels)','type','int','default',15);...
        struct('name','Frame','type','int','default',5);...
        struct('name','out bound dist','type','int','default',200);...
        struct('name','Window th. (GFP|RFP)','type','str','default','27,27');...
        struct('name','Figure Output?','type','checkbox','default',1);...
        };

    Param = Settings_GUI(S);
    suffix=char(Param(1));
    prlow=Param(2);
    str = regexprep(prlow,',',' ');
    prlow = str2num(cell2mat(str));
    en1=prlow(1);en2=prlow(2);
    % % bin=Param(3);
    % % str1 = regexprep(bin,',',' ');
    % % bin = str2num(cell2mat(str1));bin1=bin(1);bin2=bin(2);
    mag=cell2mat(Param(3));
    if mag==10
        mag=0.65;
    elseif mag==20
        mag=0.32;
    elseif mag==40
        mag=0.16;
    else  %4
        mag=1.63;
    end
    d_high=cell2mat(Param(4));
    d_low=cell2mat(Param(5));
    size_high=round(pi*((d_high/2)^2)./mag);
    size_low=round(pi*((d_low/2)^2)./mag);
    d_over=cell2mat(Param(6));
    frame=cell2mat(Param(7));
    disty=cell2mat(Param(8));
    prlow=Param(9);
    str = regexprep(prlow,',',' ');
    prlow = str2num(cell2mat(str));
    wind1=prlow(1);wind2=prlow(2);
    fig_flag=cell2mat(Param(10));
    if fig_flag==1
        fig_fig='on';
    else
        fig_fig='off';
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% rabies lin
elseif exp==5
    exp="Rabies-No DAPI";
    S = {...
        struct('name','File Suffix','type','str','default','.');...
        struct('name','*factor GFP|RFP','type','str','default','2.4,2.4');...
        struct('name','Magnification','type','int','default',4);...
        struct('name','Diameter(μ): High ','type','int','default',15);...
        struct('name','              Low','type','int','default',6);...
        struct('name','Overlap(pixels)','type','int','default',15);...
        struct('name','Frame','type','int','default',5);...
        struct('name','out bound dist','type','int','default',100);...
        struct('name','Figure Output?','type','checkbox','default',1);...
        };

    Param = Settings_GUI(S);
    suffix=char(Param(1));
    prlow=Param(2);
    str = regexprep(prlow,',',' ');
    prlow = str2num(cell2mat(str));en1=prlow(1);en2=prlow(2);
    % % bin=Param(3);
    % % str1 = regexprep(bin,',',' ');
    % % bin = str2num(cell2mat(str1));bin1=bin(1);bin2=bin(2);
    mag=cell2mat(Param(3));
    if mag==10
        mag=0.65;
    elseif mag==20
        mag=0.32;
    elseif mag==40
        mag=0.16;
    else  %4
        mag=1.63;
    end
    d_high=cell2mat(Param(4));
    d_low=cell2mat(Param(5));
    size_high=round(pi*((d_high/2)^2)./mag);
    size_low=round(pi*((d_low/2)^2)./mag);
    d_over=cell2mat(Param(6));
    frame=cell2mat(Param(7));
    disty=cell2mat(Param(8));
    fig_flag=cell2mat(Param(9));
    if fig_flag==1
        fig_fig='on';
    else
        fig_fig='off';
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else %exp==1  %%% Rabies! etay
    exp="2-IHC";
    S = {...
        struct('name','File Suffix','type','str','default','geo.');...
        struct('name','*factor GFP|RFP','type','str','default','2.4,2.4');...
        % struct('name','*Binnirization[12]','type','str','default','99.90,99.90');...
        struct('name','Magnification','type','int','default',4);...
        struct('name','Diameter(μ): High ','type','int','default',15);...
        struct('name','              Low','type','int','default',6);...
        struct('name','Overlap(pixels)','type','int','default',15);...
        struct('name','Frame','type','int','default',5);...
        struct('name','out bound dist','type','int','default',100);...
        struct('name','Figure Output?','type','checkbox','default',1);...
        % struct('name','Preview?','type','checkbox','default',1);...

        };

    Param = Settings_GUI(S);
    suffix=char(Param(1));
    prlow=Param(2);
    str = regexprep(prlow,',',' ');
    prlow = str2num(cell2mat(str));en1=prlow(1);en2=prlow(2);
    % % bin=Param(3);
    % % str1 = regexprep(bin,',',' ');
    % % bin = str2num(cell2mat(str1));bin1=bin(1);bin2=bin(2);
    mag=cell2mat(Param(3));
    if mag==10
        mag=0.65;
    elseif mag==20
        mag=0.32;
    elseif mag==40
        mag=0.16;
    else  %4
        mag=1.63;
    end
    d_high=cell2mat(Param(4));
    d_low=cell2mat(Param(5));
    size_high=round(pi*((d_high/2)^2)./mag);
    size_low=round(pi*((d_low/2)^2)./mag);
    d_over=cell2mat(Param(6));
    frame=cell2mat(Param(7));
    disty=cell2mat(Param(8));
    fig_flag=cell2mat(Param(9));
    if fig_flag==1
        fig_fig='on';
    else
        fig_fig='off';
    end

end % end of swithes exp
%%%%

preview=0;
if preview==1
    prev= input('Use previous preview(y)?','s');
    if prev=="y"
        [file,pathper] = uigetfile({'*_preview.mat'}, 'Please select preview file');
        load([pathper,'/',file], 'medians')
    else
        prevnum= input('How many dots per condition?');
        fprintf('Please Select preview'); fprintf(2, ' Image/s \n');
        listy=uipickfiles('num',[],'FilterSpec',['*',suffix,'tif'],'out','struct');
        clc
        v_rfp=[];v_gfp=[];b_rfp=[];b_gfp=[];
        for f=1:length(listy)
            cd(listy(f).folder);
            im_tif = string(listy(f).name);
            filename=strsplit(im_tif,'/');
            filename=char(filename(end));
            bfr= BioformatsImage(filename);
            ch_color=["DAPI","GFP","RFP","CY5"];

            for j=2:bfr.sizeT
                img=imread(filename,'index',j);
                figure('Name',num2str(filename),'NumberTitle','off'); imshow(img,[]);
                axis equal;
                axis off;
                hold on;
                [h, w, ~] = size(img);
                imgzoompan(gcf, 'ImgWidth', w, 'ImgHeight', h);
                temp=imcontrast;temp.Position = [0 0 300 300];
                title([char(ch_color(j)),': please slelect',num2str(prevnum) ,'LOW SIGNAL points(yellow)'])

                hA = gca;
                resetplotview(hA,'InitializeCurrentView');
                [xv, yv, buttonv] = ginputc(prevnum,'Color', 'y', 'LineWidth', 2,'LineStyle', ':');
                set(hA);
                title([char(ch_color(j)),': please slelect 10 HIGH SIGNAL points(yellow)'])
                [xv2,yv2 , buttonv2] = ginputc(prevnum,'Color', 'y', 'LineWidth', 2,'LineStyle', ':');
                set(hA);
                title([char(ch_color(j)),': please slelect 10 BACKGROUND points(red)'])
                [xb,yb , buttonb] = ginputcr(prevnum,'Color', 'r', 'LineWidth', 2,'LineStyle', ':');
                pxl_v2(f).images(j).ch=[xv2 yv2 buttonv2];
                pxl_v(f).images(j).ch=[xv yv buttonv];
                pxl_b(f).images(j).ch=[xb yb buttonb];
                close all
                for point=1:prevnum
                    dot_b=[];
                    dot_v=[];
                    dot_v2=[];
                    dot_b=[dot_b,img((round(pxl_b(f).images(j).ch(point,2))-1),round((pxl_b(f).images(j).ch(point,1)))-1)];
                    dot_b=[dot_b,img((round(pxl_b(f).images(j).ch(point,2))),round((pxl_b(f).images(j).ch(point,1)))-1)];
                    dot_b=[dot_b,img((round(pxl_b(f).images(j).ch(point,2))-1),round((pxl_b(f).images(j).ch(point,1))))];
                    dot_b=[dot_b,img((round(pxl_b(f).images(j).ch(point,2))+1),round((pxl_b(f).images(j).ch(point,1)))+1)];
                    dot_b=[dot_b,img((round(pxl_b(f).images(j).ch(point,2))),round((pxl_b(f).images(j).ch(point,1)))+1)];
                    dot_b=[dot_b,img((round(pxl_b(f).images(j).ch(point,2))+1),round((pxl_b(f).images(j).ch(point,1))))];
                    dot_b=[dot_b,img((round(pxl_b(f).images(j).ch(point,2))-1),round((pxl_b(f).images(j).ch(point,1)))+1)];
                    dot_b=[dot_b,img((round(pxl_b(f).images(j).ch(point,2))+1),round((pxl_b(f).images(j).ch(point,1)))-1)];
                    dot_b=[dot_b,img((round(pxl_b(f).images(j).ch(point,2))),round((pxl_b(f).images(j).ch(point,1))))];%middle

                    dot_v=[dot_v,img((round(pxl_v(f).images(j).ch(point,2))-1),round((pxl_v(f).images(j).ch(point,1)))-1)];
                    dot_v=[dot_v,img((round(pxl_v(f).images(j).ch(point,2))),round((pxl_v(f).images(j).ch(point,1)))-1)];
                    dot_v=[dot_v,img((round(pxl_v(f).images(j).ch(point,2))-1),round((pxl_v(f).images(j).ch(point,1))))];
                    dot_v=[dot_v,img((round(pxl_v(f).images(j).ch(point,2))+1),round((pxl_v(f).images(j).ch(point,1)))+1)];
                    dot_v=[dot_v,img((round(pxl_v(f).images(j).ch(point,2))),round((pxl_v(f).images(j).ch(point,1)))+1)];
                    dot_v=[dot_v,img((round(pxl_v(f).images(j).ch(point,2))+1),round((pxl_v(f).images(j).ch(point,1))))];
                    dot_v=[dot_v,img((round(pxl_v(f).images(j).ch(point,2))-1),round((pxl_v(f).images(j).ch(point,1)))+1)];
                    dot_v=[dot_v,img((round(pxl_v(f).images(j).ch(point,2))+1),round((pxl_v(f).images(j).ch(point,1)))-1)];
                    dot_v=[dot_v,img((round(pxl_v(f).images(j).ch(point,2))),round((pxl_v(f).images(j).ch(point,1))))];%middle

                    dot_v2=[dot_v2,img((round(pxl_v2(f).images(j).ch(point,2))-1),round((pxl_v2(f).images(j).ch(point,1)))-1)];
                    dot_v2=[dot_v2,img((round(pxl_v2(f).images(j).ch(point,2))),round((pxl_v2(f).images(j).ch(point,1)))-1)];
                    dot_v2=[dot_v2,img((round(pxl_v2(f).images(j).ch(point,2))-1),round((pxl_v2(f).images(j).ch(point,1))))];
                    dot_v2=[dot_v2,img((round(pxl_v2(f).images(j).ch(point,2))+1),round((pxl_v2(f).images(j).ch(point,1)))+1)];
                    dot_v2=[dot_v2,img((round(pxl_v2(f).images(j).ch(point,2))),round((pxl_v2(f).images(j).ch(point,1)))+1)];
                    dot_v2=[dot_v2,img((round(pxl_v2(f).images(j).ch(point,2))+1),round((pxl_v2(f).images(j).ch(point,1))))];
                    dot_v2=[dot_v2,img((round(pxl_v2(f).images(j).ch(point,2))-1),round((pxl_v2(f).images(j).ch(point,1)))+1)];
                    dot_v2=[dot_v2,img((round(pxl_v2(f).images(j).ch(point,2))+1),round((pxl_v2(f).images(j).ch(point,1)))-1)];
                    dot_v2=[dot_v2,img((round(pxl_v2(f).images(j).ch(point,2))),round((pxl_v2(f).images(j).ch(point,1))))];%middle
                    ch(j).fv(point+(f*prevnum-prevnum),1:9)=dot_v;ch(j).fb(point+(f*prevnum-prevnum),1:9)=dot_b; ch(j).fv2(point+(f*prevnum-prevnum),1:9)=dot_v2;
                end

            end

        end
        gfp_medv=median(median(ch(2).fv)); rfp_medv=median(median(ch(3).fv));
        gfp_medv2=median(median(ch(2).fv2)); rfp_medv2=median(median(ch(3).fv2));
        gfp_medb=median(median(ch(2).fb)); rfp_medb=median(median(ch(3).fb));
        medians=[gfp_medv rfp_medv gfp_medv2 rfp_medv2 gfp_medb rfp_medb];
        save([num2str(f),'_preview.mat'],'medians');
    end

end


%%%
% code for inspecting the data---
%  dot_b=zeros(length(listy)*(bfr.sizeT-1)*10,1);
%  dot_v=dot_b;
%  for f=1:length(listy)
%    for j=2:bfr.sizeT
%        for point=1:10
% dot_v(f*j*point)=img(round(pxl_b(1).images(2).ch(1,1)),round(pxl_b(1).images(2).ch(2,2)));
% dot_b(f*j*point)=img(round(pxl_b(1).images(2).ch(1,1)),round(pxl_b(1).images(2).ch(2,2)));
% end
% end
%  end
%                   P = {...
%                   struct('name','pixel low','type','str','default','(2170, 3286)  209');...
%                   struct('name','pixel high ','type','str','default','(2170, 3286)  5038');...
%                   struct('name','diameter low','type','int','default',5);...
%                   struct('name','diameter high','type','int','default',10);...
%                   };
%                   imtool(img);
%                   prev= Settings_GUI(P);
%                   pl=char(prev(1));
%                   pl = str2double(pl((end-3):end));
%                   ph=char(prev(2));
%                   ph = str2double(ph((end-3):end));
%                   dl=prev(3);
%                   dh=prev(4);
%                   pxl_v(f).images(j).ch=[pl ph dl dh];
%                   close all force
%                   end
%              end
%  else-----
%%%

% catch ME
% close all
% cd(uuu);
% return;
% end

close all
clc
warning('on','all')

%% pick files
% and count
methodc=1;
if exp=="smFISH"
    fprintf('Please pick a '); fprintf(2, 'Folder\n');
    listy=uipickfiles('num',[],'FilterSpec',['*',suffix,'tif'],'out','struct');
    clc

    for i=1:length(listy)
        cd(listy(i).name)
        I = dir([listy(i).name filesep ['*',suffix,'tif']]); % get the contents of the image_folder
        image=I(1).name;image
        path=[listy(i).name,'/'];
        cd(path);
        sprintf([num2str(i),' out of ',num2str(length(listy))])
        load('chaname.mat', 'chaname')
        ch1=char(chaname(1));    ch2=char(chaname(2));    ch3=char(chaname(3));     ch4=char(chaname(4));
        count_fish(path,image,ch1,ch2,ch3,ch4,manual_th,fig_fig,en1,en2,en3,pmt1,pmt2,pmt3,size_low,size_high)
        close all
        disp('Finished')

    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif exp=="Rabies-3ch"
    fprintf('Please Select '); fprintf(2, 'Image/s \n');
    listy=uipickfiles('num',[],'FilterSpec',['*',suffix,'*tif'],'out','struct');
    clc

    for j=1:length(listy)
        path=[listy(1).folder,'/'];
        cd(path);

        im_tif = string(listy(j).name);
        filename=strsplit(im_tif,'/');
        filename=char(filename(end)); filename
        %         if disty>0
        %             load([filename(1:(end-8)),'_maskcrop.mat'], 'zmask')
        %             % xmask=zmask==1;
        %         else
        zmask=0;
        %         end
        sprintf([num2str(j),' out of ',num2str(length(listy))])
        [numobjC1,numobjC2,numobjC1_2] = count_rabies_3ch(filename,path,fig_flag,size_low,size_high,en1,en2,en3,d_over,methodc,frame,disty,zmask);%abs[] prev


    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif   exp=="FDISCO"
    fprintf('Please Select '); fprintf(2, 'Tiff Folder \n');
    cd('/bigdata/microscope_images/Etay/FDISCO');
    %     listy=uipickfiles('num',[],'FilterSpec',['*0.tif'],'out','struct');
    trial=uigetdir();
    clc
    zdir = dir([trial filesep ['*Ch0.tif']]);
    znms=zdir(1).name;
    cd(trial);
    maxZ=length(imfinfo(znms))% ALL
    cd('/bigdata/microscope_images/Muhammad/deocy_numerator');
    fprintf('Please Select '); fprintf(2, 'Z Slice/s \n');
    listyx=uipickfiles('num',[],'FilterSpec',[],'out','struct');
    cd(trial);

    zdir1 = dir([trial filesep ['*Ch1.tif']]);
    znms1=zdir1(1).name;
        [namlistyx,idlistyx]=natsort({listyx.name});
listyx=listyx(idlistyx,:);
    %%  (optinal) create pseudozmasks
    cd(trial);
    maskq=input('(optinal) Create pseudozmasks(y)? ','s');
    if maskq=="y"
        se = strel('disk',20);
        listyx=listyx(1:maxZ);
%         namlistyx=natsort({listyx.name})';
        for is= 1:size(listyx,1)
            is
            img=imread(znms,'index',is);
            xmask=imbinarize(img,'global');
            xmask=imopen(xmask,se);
            xmask = imfill(xmask,'holes');
            %             xmask = imclearborder(xmask,4); % problematic
            %             xmask = bwareafilt(xmask, 1); % Extract largest blob.- lower cell count

            zmask = imerode(xmask,se);
            filename=strsplit(listyx(is).name,'/');
            filename=char(filename(end));
            save([filename,'_maskcrop.mat'],'zmask')

        end
    end
    %%
    dis=1;
    cd(trial);
    for  jz=1:1:size(listyx,1)

        im_tif = string(listyx(jz).name);
        filename=strsplit(im_tif,'/');
        filename=char(filename(end));
        filename
        if str2double(filename)<=maxZ
        if disty>0
            load([filename,'_maskcrop.mat'], 'zmask')
            % xmask=zmask==1;
        else
            zmask=0;
        end
        sprintf([num2str(jz),' out of ',num2str(length(listyx)/dis)])
        % % countfdisco for compiled tiff
        [numobjC1,numobjC2,numobjC1_2] = count_fdisco(znms,trial,fig_flag,size_low,size_high,en1,en2,d_over,methodc,frame,disty,zmask,wind1,wind2,str2double(filename));%abs[] prev
        %         [numobjC1,numobjC2,numobjC1_2] = count_fdisco(filename,path,fig_flag,size_low,size_high,en1,en2,d_over,methodc,frame,disty,zmask,wind1,wind2);%abs[] prev
        % % countfdisco2 for seperate tiff
        end
    end % images

elseif   exp=="Rabies-No DAPI"
    fprintf('Please Select '); fprintf(2, 'Image/s \n');
    listy=uipickfiles('num',[],'FilterSpec',['*',suffix,'*tif'],'out','struct');
    clc

    for j=1:length(listy)
        path=[listy(1).folder,'/'];
        cd(path);

        im_tif = string(listy(j).name);
        filename=strsplit(im_tif,'/');
        filename=char(filename(end)); filename
        %         if disty>0
        %             load([filename(1:(end-8)),'_maskcrop.mat'], 'zmask')
        %             % xmask=zmask==1;
        %         else
        zmask=0;
        %         end
        sprintf([num2str(j),' out of ',num2str(length(listy))])
        [numobjC1,numobjC2,numobjC1_2] = count_rabies_nodapi(filename,path,fig_flag,size_low,size_high,en1,en2,d_over,methodc,frame,disty,zmask);%abs[] prev


    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else   % 2-IHC
    fprintf('Please Select '); fprintf(2, 'Image/s \n');
    listy=uipickfiles('num',[],'FilterSpec',['*',suffix,'tif'],'out','struct');
    clc

    for  j=1:length(listy)
        path=[listy(1).folder,'/'];
        cd(path);

        im_tif = string(listy(j).name);
        filename=strsplit(im_tif,'/');
        filename=char(filename(end)); filename
        if disty>0
            load([filename(1:(end-8)),'_maskcrop.mat'], 'zmask')
            % xmask=zmask==1;
        else
            zmask=0;
        end
        sprintf([num2str(j),' out of ',num2str(length(listy))])
        try
%             if preview==1
%                 [numobjC1,numobjC2,numobjC1_2] = count_rabies(filename,path,fig_flag,en1,en2,bin1,bin2,size_low,size_high);%pct%
%             else %
                [numobjC1,numobjC2,numobjC1_2] = count_rabies_abs(filename,path,fig_flag,size_low,size_high,en1,en2,d_over,methodc,frame,disty,zmask);%abs[] prev
%             end
        end
    end
end
disp('Finished')

