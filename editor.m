%% editor
clear all
starter = input('From which Image to start?');%%%%%%%%%%%%5
% rot = input('What is your prefered rotation? ');
suffix = input('What is your image suffix?','s');
col_ch = input('show all channels(y)?','s');
rot=0; 
fprintf('Please Select '); fprintf(2, 'Folder/s \n');
listy=uipickfiles('num',[],'FilterSpec',['*',suffix,'.tif'],'out','struct'); % loop images
clc
%%
for lis=1:length(listy)
    trial=listy(lis).name;
    cd(trial)
    image_file_names = dir([trial filesep ['*',suffix,'.tif']]); % get the contents of the image_folder
    image_file_names = natsortfiles({image_file_names.name});
    im_tif = string(image_file_names);
    image_ex_names = dir([trial filesep '*_ed.tif']); % get the contents of the edited folder
    image_geo_names = dir([trial filesep '*_geo.tif']); % get the contents of the edited folder
    try
        image_ex_names = natsortfiles({image_ex_names.name});
        ex_tif = string(image_ex_names);
        for ee=1:length(ex_tif)
            image_co_names=find(im_tif==ex_tif(ee));
            im_tif([image_co_names])=[];% exclude _ed folders
        end
    end
    try
        image_geo_names = natsortfiles({image_geo_names.name});
        geo_tif = string(image_geo_names);
        for ee=1:length(geo_tif)
            image_co_names=find(im_tif==geo_tif(ee));
            im_tif([image_co_names])=[];% exclude _ed folders
        end
    end
    file_numCol=length(im_tif);
    reditx=1;
    for u= starter:file_numCol
        disp("Reading image...")
        ch=char(im_tif(u));
        [~,name,~] = fileparts(ch);
        img=imread(ch,'index',1);imgx=img;
        xx=1;
        try
            img2=imread(ch,'index',2);img2x=img2;
            xx=2
        end
        try
            img3=imread(ch,'index',3);img3x=img3;
            xx=3
        end
        try
            img4=imread(ch,'index',4);img4x=img4;
            xx=4
        end
        % J=imadjust(J);
        sprintf([num2str(u-starter+1),' image out of ',num2str(file_numCol-starter+1),'\nfolder ',num2str(lis),' out of ',num2str(length(listy)), ' folder/s'])
        k=0;
        while k<1
            clear J Jimg Jimg2 Jimg3 Jimg4
            
            if col_ch=="y"
                
                if xx==3
                    disp('Loading multichannel image...');
                    mainy=figure('Name',num2str(name),'NumberTitle','off');mainy.Position=[0 500 2000 500];
                    tiledlayout(1,3)
                    %     subplot(1,3,1)
                    ax1 = nexttile;
                    f1=imshow(img,[]);
                    title('Ch1')
                    %     subplot(1,3,2)
                    ax2 = nexttile;
                    f2=imshow(img2,[]);
                    title('Ch2')
                    %     subplot(1,3,3)
                    ax3 = nexttile;
                    f3=imshow(img3,[]);
                    title('Ch3')
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
                    Jimg=imadjust(img,[cont1/65535]);Jimg2=imadjust(img2,[cont2/65535]);Jimg3=imadjust(img3,[cont3/65535]);
                    J = cat(3, Jimg3, Jimg, Jimg2);
                    figure;imshow(J);title('\color{green}DAPI, \color{blue}GFP, \color{red}RFP')
                    
                else % four channels
                    disp('Loading multichannel image...');
                    mainy=figure('Name',num2str(name),'NumberTitle','off');mainy.Position=[500 0 1000 1500];
                    tiledlayout(2,2)
                    ax1 = nexttile;
                    %     subplot(2,2,1)
                    f1=imshow(img,[]);
                    title("Ch1")
                    ax2 = nexttile;
                    %     subplot(2,2,2)
                    f2=imshow(img2,[]);
                    title("Ch2")
                    ax3 = nexttile;
                    %     subplot(2,2,3)
                    f3=imshow(img3,[]);
                    title("Ch3")
                    ax4 = nexttile;
                    %         subplot(2,2,4)
                    f4=imshow(img4,[]);
                    title("Ch4")
                    annotation(mainy,'textbox',...
                        [0.212580086580087 0.914807303984314 0.571699117322441 0.0791075033179162],...
                        'String','Right mouse click UP and DOWN to contrast | press SPACE to save contrast and Finish',...
                        'LineStyle','none',...
                        'FontSize',18);
                    temp1=imcontrast(f1);temp1.Position = [100 500 300 300];
                    temp2=imcontrast(f2);temp2.Position = [1600 500 300 300];
                    temp3=imcontrast(f3);temp3.Position = [100 100 300 300];
                    temp4=imcontrast(f4);temp4.Position = [1600 100 300 300];
                    %          [h, w, ~] = size(img);
                    % imgzoompan(gcf, 'ImgWidth', w, 'ImgHeight', h);
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
                    figure;imshow(J);
                    load('chaname.mat', 'chaname');
                    title([['\color{gray}',num2str(chaname(1))],' ',['\color{blue}',num2str(chaname(2))],' ',['\color{red}',num2str(chaname(3))],' ',['\color{green}',num2str(chaname(4))]]);
                    
                end
            else
                
                mainy=figure('Name',[num2str(name)],'NumberTitle','off');
                J=im2gray(img);
                imshow(J,[]);
                title('Ch: Right mouse click UP and DOWN to contrast| press SPACE to save contrast and Finish');
                temp=imcontrast;temp.Position = [100 100 300 300];
                pause;
                cmax = str2num(temp.Children(1).Children(3).Children.Children(2).Children.Children(2).Children(2).String);
                cmin = str2num(temp.Children(1).Children(3).Children.Children(2).Children.Children(2).Children(5).String);
                cont1=[cmin cmax];
                J=imadjust(img,[cont1/65535]);
                
            end
            
            prompt = {'1.Rotation(CCW):','4.Cropping Out','5.Cropping In','2.Flip(hv)','3.Stich'};
            dlgtitle = 'Editor';
            dims = [1 35];
            definput = {num2str(rot),'ok','ok','hv','ok'};
            answer = string(inputdlg(prompt,dlgtitle,dims,definput));
            if isempty(answer)==1
                close all
                k=1;
            else
                rot=str2double(answer(1));
                J=imrotate(J,rot);
                if answer(4)=="h"
                    J= flip(J ,2);
                elseif answer(4)=="v"
                    J= flip(J ,1);
                elseif answer(4)=="hv"||answer(4)=="vh"
                    J= flip(J ,2);
                    J= flip(J ,1);
                end
                if answer(5)=="ok" %stich
                    %                     Jst=J;
                    close all
                    fst=figure;
                    hIm=imshow(J,[]);
                    sticher;
                    uiwait(fst);
                    try
                        CurrAng=round(CurrAng);
                    catch
                        CurrAng=0;
                    end
                    infox = regionprops(bw1,'Boundingbox') ;BBy = infox(1).BoundingBox;
                    [crim,~] = imcrop(J,BBy);
                    [crBW,~] = imcrop(bw1,BBy);
                    crim(crBW==0)=0;%real subimage
                    J(bw1==1)=0;
                    xs=CurrPos(:,1);ys=CurrPos(:,2);
                    [ms,ns]=size(J);
                    BWp=poly2mask(xs,ys,ms,ns);
                    infox = regionprops(BWp,'Boundingbox') ;BBx = infox(1).BoundingBox;
                    trim=imtranslate(crim,[BBx(1) BBx(2)],'FillValues',0,'OutputView','full');
                    roim=imrotate(trim,-CurrAng,'crop');
                    BWp=imrotate(BWp,-CurrAng,'crop');
                    J(BWp==1)=0;
                    Jst=J(:,:,1);
                    J=imfuse(roim,Jst,'diff','Scaling','independent');
                end %stich
                if answer(2)=="ok"
                    mask_zero;
                end
                if answer(3)=="ok"
                    [croppedimage,rectout] = imcrop(J);
                    J=croppedimage;
                end
                close all
                figure('Name',[num2str(name)],'NumberTitle','off');
                imshow(J)
                mTitle = 'Another Round';
                mQuestion = 'Do you wish to re-edit? ';
                choice = string(questdlg(mQuestion,mTitle,'New','re-edit', 'Next','Next'));
                close all
                
                if choice=="Next"
                    k=1;
                    reidtx=0;
                elseif choice=="re-edit"
                    k=1;
                    reidtx=1;
                else % new
                    reidtx=-1;
                end
                
                close all
                if reidtx==0||reidtx==1
                    cd(trial)  % start entering to each channel
                    disp('Writing multichannel image')
                    
                    R1=imrotate(imgx,rot);
                    if answer(4)=="h"
                        R1= flip(R1 ,2);
                    elseif answer(4)=="v"
                        R1= flip(R1 ,1);
                    elseif answer(4)=="hv"||answer(4)=="vh"
                        R1= flip(R1 ,2);
                        R1= flip(R1 ,1);
                    end
                    if answer(5)=="ok"
                        [crim,~] = imcrop(R1,BBy);
                        R1(Jst==0)=0;
                        crim(crBW==0)=0;   %real subimage
                        trim=imtranslate(crim,[BBx(1) BBx(2)],'FillValues',0,'OutputView','full');
                        roim=imrotate(trim,-CurrAng,'crop');
                        R1=imfuse(R1,roim,'diff','Scaling','independent');
                    end
                    if answer(2)=="ok"
                        R1(cumulativeBinaryImage)=0;
                    end
                    if answer(3)=="ok"
                        cropy=floor(rectout(1));
                        cropx=floor(rectout(2));
                        if cropx==0
                            cropx=1;
                        end
                        if cropy==0
                            cropy=1;
                        end
                        R1=R1(floor(cropx:(cropx+rectout(4))),floor(cropy:cropy+rectout(3)));
                    end
                    if reidtx==0
                    imwrite(R1,[num2str(name),'_ed','.tif'],'WriteMode','append');
                    end
                    
                    disp('Edit Channel 1');
                    if xx>1
                        R2=imrotate(img2x,rot);
                        if answer(4)=="h"
                            R2= flip(R2 ,2);
                        elseif answer(4)=="v"
                            R2= flip(R2 ,1);
                        elseif answer(4)=="hv"||answer(4)=="vh"
                            R2= flip(R2 ,2);
                            R2= flip(R2 ,1);
                        end
                        if answer(5)=="ok"
                            [crim,~] = imcrop(R2,BBy);
                            R2(Jst==0)=0;
                            crim(crBW==0)=0;   %real subimage
                            trim=imtranslate(crim,[BBx(1) BBx(2)],'FillValues',0,'OutputView','full');
                            roim=imrotate(trim,-CurrAng,'crop');
                            R2=imfuse(R2,roim,'diff','Scaling','independent');
                        end
                        if answer(2)=="ok"
                            R2(cumulativeBinaryImage)=0;
                        end
                        if answer(3)=="ok"
                            R2=R2(floor(rectout(2):(rectout(2)+rectout(4))),floor(rectout(1):rectout(1)+rectout(3)));
                        end
                        if reidtx==0
                        imwrite(R2,[num2str(name),'_ed','.tif'],'WriteMode','append');
                        end
                        disp('Edit Channel 2');
                        
                    end
                    if xx>2
                        R3=imrotate(img3x,rot);
                        if answer(4)=="h"
                            R3= flip(R3 ,2);
                        elseif answer(4)=="v"
                            R3= flip(R3 ,1);
                        elseif answer(4)=="hv"||answer(4)=="vh"
                            R3= flip(R3 ,2);
                            R3= flip(R3 ,1);
                        end
                        if answer(5)=="ok"
                            [crim,~] = imcrop(R3,BBy);
                            R3(Jst==0)=0;
                            crim(crBW==0)=0;   %real subimage
                            trim=imtranslate(crim,[BBx(1) BBx(2)],'FillValues',0,'OutputView','full');
                            roim=imrotate(trim,-CurrAng,'crop');
                            R3=imfuse(R3,roim,'diff','Scaling','independent');
                        end
                        if answer(2)=="ok"
                            R3(cumulativeBinaryImage)=0;
                        end
                        if answer(3)=="ok"
                            R3=R3(floor(rectout(2):(rectout(2)+rectout(4))),floor(rectout(1):rectout(1)+rectout(3)));
                        end
                         if reidtx==0
                        imwrite(R3,[num2str(name),'_ed','.tif'],'WriteMode','append');
                        disp('Edit Channel 3');
                         end
                    end
                    if xx>3
                        R4=imrotate(img4x,rot);
                        if answer(4)=="h"
                            R4= flip(R4 ,2);
                        elseif answer(4)=="v"
                            R4= flip(R4 ,1);
                        elseif answer(4)=="hv"||answer(4)=="vh"
                            R4= flip(R4 ,2);
                            R4= flip(R4 ,1);
                        end
                        if answer(5)=="ok"
                            [crim,~] = imcrop(R4,BBy);
                            R4(Jst==0)=0;
                            crim(crBW==0)=0;   %real subimage
                            trim=imtranslate(crim,[BBx(1) BBx(2)],'FillValues',0,'OutputView','full');
                            roim=imrotate(trim,-CurrAng,'crop');
                            R4=imfuse(R4,roim,'diff','Scaling','independent');
                        end
                        if answer(2)=="ok"
                            R4(cumulativeBinaryImage)=0;
                        end
                        if answer(3)=="ok"
                            R4=R4(floor(rectout(2):(rectout(2)+rectout(4))),floor(rectout(1):rectout(1)+rectout(3)));
                        end
                        if reidtx==0
                        imwrite(R4,[num2str(name),'_ed','.tif'],'WriteMode','append');
                        end
                        
                        disp('Edit Channel 4');
                        
                    end
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if reidtx==1
                   while 1
                    cd(trial)  % start entering to each channel
                    ch=[name,'_ed.tif'];
                    mainy=figure('Name',[num2str(name)],'NumberTitle','off');
                    J=im2gray(R1);
                    imshow(J,[]);
                    title('Ch: Right mouse click UP and DOWN to contrast| press SPACE to save contrast and Finish');
                    temp=imcontrast;temp.Position = [100 100 300 300];
                    pause;
                    cmax = str2num(temp.Children(1).Children(3).Children.Children(2).Children.Children(2).Children(2).String);
                    cmin = str2num(temp.Children(1).Children(3).Children.Children(2).Children.Children(2).Children(5).String);
                    cont1=[cmin cmax];
                    J=imadjust(R1,[cont1/65535]);
                    imshow(J,[]);
%                     title('Ch: Right mouse click UP and DOWN to contrast| press SPACE to save contrast and Finish');
%                     temp=imcontrast;temp.Position = [100 100 300 300];
%                     pause;
%                     cmax = str2num(temp.Children(1).Children(3).Children.Children(2).Children.Children(2).Children(2).String);
%                     cmin = str2num(temp.Children(1).Children(3).Children.Children(2).Children.Children(2).Children(5).String);
%                     cont1=[cmin cmax];
%                     J=imadjust(img,[cont1/65535]);
                    prompt = {'1.Rotation(CCW):','4.Cropping Out','5.Cropping In','2.Flip(hv)','3.Stich'};
                    dlgtitle = 'Editor';
                    dims = [1 35];
                    definput = {num2str(rot),'ok','ok','hv','ok'};
                    answer = string(inputdlg(prompt,dlgtitle,dims,definput));
                    if isempty(answer)==1
                        close all
                        k=1;
                    else
                        rot=str2double(answer(1));
                        J=imrotate(J,rot);
                        if answer(4)=="h"
                            J= flip(J ,2);
                        elseif answer(4)=="v"
                            J= flip(J ,1);
                        elseif answer(4)=="hv"||answer(4)=="vh"
                            J= flip(J ,2);
                            J= flip(J ,1);
                        end
                        if answer(5)=="ok" %stich
                            %                     Jst=J;
                            close all
                            fst=figure;
                            hIm=imshow(J,[]);
                            sticher;
                            uiwait(fst);
                            try
                                CurrAng=round(CurrAng);
                            catch
                                CurrAng=0;
                            end
                            infox = regionprops(bw1,'Boundingbox') ;BBy = infox(1).BoundingBox;
                            [crim,~] = imcrop(J,BBy);
                            [crBW,~] = imcrop(bw1,BBy);
                            crim(crBW==0)=0;%real subimage
                            J(bw1==1)=0;
                            xs=CurrPos(:,1);ys=CurrPos(:,2);
                            [ms,ns]=size(J);
                            BWp=poly2mask(xs,ys,ms,ns);
                            infox = regionprops(BWp,'Boundingbox') ;BBx = infox(1).BoundingBox;
                            trim=imtranslate(crim,[BBx(1) BBx(2)],'FillValues',0,'OutputView','full');
                            roim=imrotate(trim,-CurrAng,'crop');
                            BWp=imrotate(BWp,-CurrAng,'crop');
                            J(BWp==1)=0;
                            Jst=J(:,:,1);
                            J=imfuse(roim,Jst,'diff','Scaling','independent');
                        end %stich
                        if answer(2)=="ok"
                            mask_zero;
                        end
                        if answer(3)=="ok"
                            [croppedimage,rectout] = imcrop(J);
                            J=croppedimage;
                        end
                        figure('Name',[num2str(name),'_ed'],'NumberTitle','off');
                        imshow(J)
                        disp('Writing multichannel image')
                        R1=imrotate(R1,rot);
                        if answer(4)=="h"
                            R1= flip(R1 ,2);
                        elseif answer(4)=="v"
                            R1= flip(R1 ,1);
                        elseif answer(4)=="hv"||answer(4)=="vh"
                            R1= flip(R1 ,2);
                            R1= flip(R1 ,1);
                        end
                        if answer(5)=="ok"
                            [crim,~] = imcrop(R1,BBy);
                            R1(Jst==0)=0;
                            crim(crBW==0)=0;   %real subimage
                            trim=imtranslate(crim,[BBx(1) BBx(2)],'FillValues',0,'OutputView','full');
                            roim=imrotate(trim,-CurrAng,'crop');
                            R1=imfuse(R1,roim,'diff','Scaling','independent');
                        end
                        if answer(2)=="ok"
                            R1(cumulativeBinaryImage)=0;
                        end
                        if answer(3)=="ok"
                            cropy=floor(rectout(1));
                            cropx=floor(rectout(2));
                            if cropx==0
                                cropx=1;
                            end
                            if cropy==0
                                cropy=1;
                            end
                            R1=R1(floor(cropx:(cropx+rectout(4))),floor(cropy:cropy+rectout(3)));
                        end
                        close all
                        imwrite(R1,[num2str(name),'_ed','.tif'],'WriteMode','append');
                        
                        
                        disp('Edit Channel 1');
                        if xx>1
                            R2=imrotate(R2,rot);
                            if answer(4)=="h"
                                R2= flip(R2 ,2);
                            elseif answer(4)=="v"
                                R2= flip(R2 ,1);
                            elseif answer(4)=="hv"||answer(4)=="vh"
                                R2= flip(R2 ,2);
                                R2= flip(R2 ,1);
                            end
                            if answer(5)=="ok"
                                [crim,~] = imcrop(R2,BBy);
                                R2(Jst==0)=0;
                                crim(crBW==0)=0;   %real subimage
                                trim=imtranslate(crim,[BBx(1) BBx(2)],'FillValues',0,'OutputView','full');
                                roim=imrotate(trim,-CurrAng,'crop');
                                R2=imfuse(R2,roim,'diff','Scaling','independent');
                            end
                            if answer(2)=="ok"
                                R2(cumulativeBinaryImage)=0;
                            end
                            if answer(3)=="ok"
                                R2=R2(floor(rectout(2):(rectout(2)+rectout(4))),floor(rectout(1):rectout(1)+rectout(3)));
                            end
                            imwrite(R2,[num2str(name),'_ed','.tif'],'WriteMode','append');
                            
                            disp('Edit Channel 2');
                            
                        end
                        if xx>2
                            R3=imrotate(R3,rot);
                            if answer(4)=="h"
                                R3= flip(R3 ,2);
                            elseif answer(4)=="v"
                                R3= flip(R3 ,1);
                            elseif answer(4)=="hv"||answer(4)=="vh"
                                R3= flip(R3 ,2);
                                R3= flip(R3 ,1);
                            end
                            if answer(5)=="ok"
                                [crim,~] = imcrop(R3,BBy);
                                R3(Jst==0)=0;
                                crim(crBW==0)=0;   %real subimage
                                trim=imtranslate(crim,[BBx(1) BBx(2)],'FillValues',0,'OutputView','full');
                                roim=imrotate(trim,-CurrAng,'crop');
                                R3=imfuse(R3,roim,'diff','Scaling','independent');
                            end
                            if answer(2)=="ok"
                                R3(cumulativeBinaryImage)=0;
                            end
                            if answer(3)=="ok"
                                R3=R3(floor(rectout(2):(rectout(2)+rectout(4))),floor(rectout(1):rectout(1)+rectout(3)));
                            end
                            imwrite(R3,[num2str(name),'_ed','.tif'],'WriteMode','append');
                            disp('Edit Channel 3');
                            
                        end
                        if xx>3
                            R4=imrotate(R4,rot);
                            if answer(4)=="h"
                                R4= flip(R4 ,2);
                            elseif answer(4)=="v"
                                R4= flip(R4 ,1);
                            elseif answer(4)=="hv"||answer(4)=="vh"
                                R4= flip(R4 ,2);
                                R4= flip(R4 ,1);
                            end
                            if answer(5)=="ok"
                                [crim,~] = imcrop(R4,BBy);
                                R4(Jst==0)=0;
                                crim(crBW==0)=0;   %real subimage
                                trim=imtranslate(crim,[BBx(1) BBx(2)],'FillValues',0,'OutputView','full');
                                roim=imrotate(trim,-CurrAng,'crop');
                                R4=imfuse(R4,roim,'diff','Scaling','independent');
                            end
                            if answer(2)=="ok"
                                R4(cumulativeBinaryImage)=0;
                            end
                            if answer(3)=="ok"
                                R4=R4(floor(rectout(2):(rectout(2)+rectout(4))),floor(rectout(1):rectout(1)+rectout(3)));
                            end
                            imwrite(R4,[num2str(name),'_ed','.tif'],'WriteMode','append');
                            disp('Edit Channel 4');
                            pause(0.5)
                        end
                        
                    end
                mTitle = 'Another Round';
                mQuestion = 'Do you wish to re-edit? ';
                choice = string(questdlg(mQuestion,mTitle,'re-edit','Next','Next'));
                if choice=="Next"
                break;
                end
                end
                end
            end
            clc
        end
    end
end
        disp("Finished")