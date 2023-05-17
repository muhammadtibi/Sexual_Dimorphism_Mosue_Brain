% extract roi
clear all
warning('off','all')
uuu=cd;
transcendentalDoc='';
% COA: COAa,COAp,
% sAMY: AAA,BA,CEAc,CEAl,CEAm,IA,MEAad,MEAav,MEApd,MEApv,PAA,
% CTXsp: LA,BLAa,BLAp,BLAv,BMAa,BMAp,PA,
% COAa,COAp,AAA,BA,CEAc,CEAl,CEAm,IA,MEAad,MEAav,MEApd,MEApv,LA,BLAa,BLAp,BLAv,BMAa,BMAp,PA,
% COA,AAA,BA,CEA,IA,MEA,LA,BLA,BMA,PA,NLOT,
clc
try
    % % % %     tt = readtable('/bigdata/microscope_images/roicombo.csv');
    % % % %     fig = uifigure('HandleVisibility', 'on');
    % % % %     fig.Name = 'Choose genes and drag to Settings Dialog';
    % % % %     fig.Position = [30 30 520 500];
    % % % %     uit = uitable(fig,'Data',tt);
    % % % %     uit.Position = [20 30 480 450];
    % % % %     uit.ColumnEditable = true;
    % % %     t = readtable('structure_tree_safe_2017.csv');
    % % %     vars = {'name','acronym','structure_id_path'};
    % % %     t = t(:,vars);
    % % %     t.index = (1:height(t)).';
    % % %     t = t(:, [4 1 2 3]);
    % % %     fig = uifigure('HandleVisibility', 'on');
    % % %     fig.Name = 'Choose ROI and drag to Settings Dialog';
    % % %     fig.Position = [30 30 520 500];
    % % %     uit = uitable(fig,'Data',t);
    % % %     uit.Position = [20 30 480 450];
    % % %     uit.ColumnEditable = true;
    S = {...
        struct('name','Stat. Level','type','enum','values',{{'Extract ROI','Click Check','Edit ROI','SpreadSheet','Extract ABCDE'}},'doc',transcendentalDoc);...
        struct('name','File Suffix','type','str','default','geo.');...
        struct('name','*Overlap Distance','type','int','default',5);...
        struct('name','*Threshold Lvl.','type','str','default','3,3,3');...
        struct('name','*Region of Interest','type','str','default','AAA,AOB,AON,BA,BAC,BLA,BMA,BST,CB,CEA,CLA,COA,DP,EP,HB,HPF,PVZ,PVR,MEZ,LZ,ME,IA,Isocortex,LSX,MB,MEA,MOB,NLOT,PA,PAA,PALd,PALm,PALv,PIR,STRV,STRd,TH,TR,TT,fiber tracts,');
        struct('name','EXP','type','enum','values',{{'4X','20X','AISH','VIS'}});...
        };


    Param = Settings_GUI(S);
    state=string(Param(1));
    suffix=cell2mat(Param(2));
    closestAllowableDistance=cell2mat(Param(3));
    th=cell2mat(Param(4));
    str = regexprep(th,',',' ');
    th = str2num(str);thr1=th(1);thr2=th(2);thr3=th(3);
    regions=string(Param(5));
    regions=strsplit(regions,',');
    regions=regions(1:end-1);
    expy=string(Param(6));
catch ME
    close all
    cd(uuu);
    return ;
end
close all
clc
%% switch cases
if state=="Crop ROI" || state== "Edit ROI"
else
    fprintf('Please pick a '); fprintf(2, 'Folder\n');
    listy=uipickfiles('num',[],'FilterSpec',['*',suffix,'tif'],'out','struct');
end
clc
% cd(listy(1).name)
if state=="Click Check"

    for lis=1:length(listy)  % go over multidirectiores
        trial=listy(lis).name;
        cd(trial)
        image_file_names = dir([listy(lis).name filesep ['*',suffix,'tif']]); % get the contents of the image_folder
        image_file_names = natsortfiles({image_file_names.name});
        [~, file_numCol] = size(image_file_names); %count number of image tif files in the directory
        starter = input('From which image to start?');

        for fol=starter:file_numCol
            sprintf([num2str(fol),' image out of ',num2str(file_numCol),'\nfolder ',num2str(lis),' out of ',num2str(length(listy)), ' folder/s'])
            cd(listy(lis).name);
            filename=char(image_file_names(fol));
            % disp('coloring TIFF as ND2 please wait')
            % bfr= BioformatsImage(filename);
            B=importdata([filename(1:(end-4)),'_position_counts_per_dapi.txt']);
            %                  bfr= BioformatsImage(diry);
            ch_color=["DAPI","GFP","RFP","CY5"];
            col=['b','g','r','magenta'];
            dot_b=[];
            dot_v=[];
            cent_x=B.data(:,1); cent_y=B.data(:,2);
            for j=2:3
                img=imread(filename,'index',j);
                figure('Name',num2str(filename),'NumberTitle','off'); zz=imshow(img,[]);
                axis equal;
                axis off;
                hold on;
                title([char(ch_color(j)),': Choose contrast then press <ENTER>'])
                [h, w, ~] = size(img);
                imgzoompan(gcf, 'ImgWidth', w, 'ImgHeight', h);
                %
                temp = imcontrast;temp.Position = [0 0 300 300];
                pause;
                cmax = str2num(temp.Children(1).Children(3).Children.Children(2).Children.Children(2).Children(2).String);
                cmin = str2num(temp.Children(1).Children(3).Children.Children(2).Children.Children(2).Children(5).String);
                title([char(ch_color(j)),': please select undetected points then <ENTER> image'])
                b_xfp=B.data(:,2+j)==1;
                %    b_rfp=B.data(:,5)==1;
                b_overlap=B.data(:,6)==1;
                scatter(B.data(b_xfp,1),B.data(b_xfp,2),100,'MarkerEdgeColor',char(col(j)),'LineWidth',2)
                % scatter(B.data(b_rfp,1),B.data(b_rfp,2),100,'d','r','LineWidth',2)

                scatter(B.data(b_overlap,1),B.data(b_overlap,2),100,'+','b')
                %    pl = selectdata('selectionmode','brush');
                hA = gca;
                resetplotview(hA,'InitializeCurrentView');
                [xv, yv,~] = ginputc('Color', 'y', 'LineWidth', 2,'LineStyle', ':');
                set(hA);
                [aa,bb]=size(B.data);
                B.data=[B.data(:,:);zeros(length(xv),bb)];
                B.data(((end-length(xv)+1):end),1)=xv;
                B.data(((end-length(yv)+1):end),2)=yv;
                B.data(((end-length(xv)+1):end),2+j)=ones(length(xv),1);
                title([char(ch_color(j)),': please deselect falsely detected points then CLOSE image'])
                non_x=[];non_y=[];
                if j==2
                    xy_gfp=[xv yv];
                    cont_gfp=[cmin cmax];
                else
                    xy_rfp=[xv yv];
                    cont_rfp=[cmin cmax];
                end
                impixelinfo;

                while size(findobj(zz))>0
                    try
                        [pl,xsel,ysel] = selectdata('selectionmode','closest','Action','delete');
                        non_x=[non_x,xsel];%beacuse its reverted
                        non_y=[non_y,ysel];%beacuse its reverted
                    end
                    %   scatter(cell2mat(xsel(2)),cell2mat(ysel(2)))
                end
                try
                    non_x= cell2mat(non_x(~cellfun('isempty',non_x)));
                    non_y= cell2mat(non_y(~cellfun('isempty',non_y)));
                end
                % cent_x=B.data((B.data(:,2+j)==1),1);
                % cent_y=B.data((B.data(:,2+j)==1),2);
                for xx=1:length(non_x)
                    cent_x=B.data(:,1); cent_y=B.data(:,2);

                    B.data(((round(cent_x,2))==round(non_x(xx),2)) & (round(cent_y,2)==round(non_y(xx),2)),:)=[];
                    % B.data(((floor(cent_x,2))==floor(non_x(xx),2)) & (floor(cent_y,2)==floor(non_y(xx),2)),:)=[];

                end

            end
            gf= imadjust(imread(filename,'index',2),[cont_gfp/65535]);                   rf= imadjust(imread(filename,'index',3),[cont_rfp/65535]);

            C = imfuse(rf,gf,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);
            rev=figure('Name',num2str(filename),'NumberTitle','off');
            ff=imshow((C));title('Overlap(pink): please deselect (+) falsely co-detected cells then CLOSE image')
            hold on
            b_gfp=B.data(:,4)==1;b_rfp=B.data(:,5)==1;b_overlap=B.data(:,6)==1;
            gfp_cor=[B.data(b_gfp,1),B.data(b_gfp,2)];
            %                    gfp_cor=gfp_cor();
            rfp_cor=[B.data(b_rfp,1),B.data(b_rfp,2)];
            % %                   scatter(B.data(b_gfp,1),B.data(b_gfp,2),100,'MarkerEdgeColor','g','LineWidth',2)
            % %                   scatter(B.data(b_rfp,1),B.data(b_rfp,2),100,'MarkerEdgeColor','r','LineWidth',2)
            scatter(B.data(b_overlap,1),B.data(b_overlap,2),100,'MarkerEdgeColor','m','LineWidth',2)
            %                    [distances,~] = pdist2(xy_gfp, xy_rfp,'euclidean','Smallest',1);
            [distances1,~] = pdist2(rfp_cor,gfp_cor,'euclidean','Smallest',1);
            [distances2,~] = pdist2(gfp_cor,rfp_cor,'euclidean','Smallest',1);
            closePairs1 = distances1< closestAllowableDistance;
            closePairs2 = distances2< closestAllowableDistance;

            gfp_over=gfp_cor(closePairs1,:)+1;                  rfp_over=rfp_cor(closePairs2,:)+1;
            scatter(rfp_over(:,1),rfp_over(:,2),100,'+','MarkerEdgeColor','m','LineWidth',2);
            scatter(gfp_over(:,1),gfp_over(:,2),10,'*','MarkerEdgeColor','m','LineWidth',2);
            [h, w, ~] = size(C);
            imgzoompan(gcf, 'ImgWidth', w, 'ImgHeight', h);
            non_x=[];non_y=[];
            B.data=[B.data(:,:);zeros(sum(closePairs2),bb)];
            B.data(((end-sum(closePairs2)+1):end),1)=rfp_over(:,1);
            B.data(((end-sum(closePairs2)+1):end),2)=rfp_over(:,2);
            B.data(((end-sum(closePairs2)+1):end),2+j+1)=ones(sum(closePairs2),1);
            while size(findobj(ff))>0
                try
                    [pl,xsel,ysel] = selectdata('selectionmode','closest','Action','delete');
                    non_x=[non_x,xsel];%beacuse its reverted
                    non_y=[non_y,ysel];%beacuse its reverted
                end
            end
            try
                non_x= cell2mat(non_x(~cellfun('isempty',non_x)));
                non_y= cell2mat(non_y(~cellfun('isempty',non_y)));
            end


            for xx=1:length(non_x)
                % cent_x=B.data((B.data(:,2+j+1)==1),1);
                % cent_y=B.data((B.data(:,2+j+1)==1),2);
                cent_x=B.data(:,1); cent_y=B.data(:,2);
                B.data(((round(cent_x))==round(non_x(xx))) & (round(cent_y)==round(non_y(xx))),:)=[];

                %  end
            end
            % save updated text
            rev_T = [ {'dapi_x','dapi_y','area','GFP','RFP','overlap','value'};m2c([B.data(:,1),B.data(:,2),B.data(:,3),B.data(:,4),B.data(:,5),B.data(:,6),B.data(:,7)])];
            saveCellFile(rev_T,[trial,'/',filename(1:end-4),'_position_counts_per_dapi','.txt']);
            clc
        end
    end

    disp("Finished")
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555
elseif state=="Extract ROI"
    % ROI extraction
    % my ROI that i am interested in for example
    if expy~="VIS"
        disp('Loading Libraries...')
        % Load CCF atlas
        allen_atlas_path = '/data/Aligment/AP_histology-master/AllenCCF';
        %addpath(genpath('/data/Technion_analysis/Amygdala/FISH/Aligment'))
        % tv = readNPY([allen_atlas_path filesep 'template_volume_10um.npy']);
        % av = readNPY([allen_atlas_path filesep 'annotation_volume_10um_by_index.npy']);
        st = loadStructureTree([allen_atlas_path filesep 'structure_tree_safe_2017.csv']);
        clc
        if expy=='AISH'
            cont=2;
        elseif expy=='20X'
            cont=1;
        else%4X
            cont=0;
        end
        disp('extracting ROI...')
        for lis=1:length(listy)  % go over multidirectiores
            trial=listy(lis).name;
            cd(trial)
            %       try
            % text_file_names1 = dir([trial filesep ['*counts_per_dapi.txt']]); % get the contents of the image_folder
            %       end
            text_file_names = dir([trial filesep '*counts_per_dapi.txt']); % get the contents of the image_folder
            text_file_names = natsortfiles({text_file_names.name});
            text_file_names = string(text_file_names);

            if regions(1)~="All"
                tonormal=[];
                for w=1:length(text_file_names) % for all the image in th folder
                    clear T G F

                    sprintf([num2str(w),' file out of ',num2str(length(text_file_names)),'\nfolder ',num2str(lis),' out of ',num2str(length(listy)), ' folder/s'])
                    fchar=char(text_file_names(w));
                    %                    try
                    cd(trial)

                    f_name=fchar(1:(end-33));
                    if cont~=2 % etay and fish
                        load([f_name,'_maskcrop','.mat'], 'zmask'); % get mask of rois
                    else %AISH
                        load([f_name(1:(end-1)),'_maskcrop','.mat'], 'zmask'); % get mask of rois
                    end
                    %             end
                    % %              try
                    % %                     f_name=fchar(1:(end-37));
                    % % load([f_name,'_maskcrop','.mat'], 'zmask'); % get mask of rois
                    % %             end
                    cd(trial)
                    T=readtable(text_file_names(w));
                    if height(T)>0
                        T(~T.dapi_x,:) = [];
                        T(T.dapi_x> size(zmask,2) ,:)=[];            T(T.dapi_y> size(zmask,1),:)=[];
                        dot_roi=zeros(height(T),1);
                        n_mask=zeros(size(zmask));
                        %       roi_mask.n_mask=n_mask;      roi_mask.regions=regions;
                        apendy=[];
                        allblobs=zeros(1,length(regions)+1);
                        for k=1:length(regions) %  subregions
                            notfound=0;
                            roi_area=regions(k)
                            roi_id_path = st.structure_id_path(find(strcmp(st.safe_name,roi_area)));%search full name
                            if isempty(roi_id_path)
                                roi_id_path = st.structure_id_path(find(strcmp(st.acronym,roi_area)));% search acronym
                            end
                            roi_idx = find(contains(st.structure_id_path,roi_id_path));% all subregions
                            xmask=ismember(zmask,roi_idx);
                            [bx,by]=find(xmask);
                            for nn=1:length(dot_roi)
                                dot_roi(nn)= zmask(floor(T.dapi_y(nn)),floor(T.dapi_x(nn)));%??? find zmask and dot
                            end
                            %                             blobmask=zeros(size(zmask));
                            for n=1:length(roi_idx)% find if its part of roi  and subroi
                                pos=dot_roi==roi_idx(n);
                                T.ROI(pos)=k;
                                notfound=notfound+sum(pos);
                                %                                 blobmask=blobmask+(zmask==roi_idx(n));
                            end
                            if notfound==0 % all subroiu is not found=0
                                apendy=[apendy,k];
                            else
                                % others
                            end
                            %              roi_mask(k).roi=[bx,by];

                            allblobs(k)=sum(sum(xmask));
                        end

                        T.ROI(~T.ROI,:)=length(regions)+1;
                        regionsx=[regions,'others'];
                        % T.ROI_names=zeros(height(T),1);
                        % T.ROI_names=regions(T.ROI)
                        %  T(~T.ROI,:)=[];
                        %   T.index = (1:height(T)).'
                        %   T = T(:, [7 1 2 3 4 5 6 7])

                        for aa=1: length(T.ROI)
                            T.ROI_names(aa)=regionsx(T.ROI(aa)); % get names of roi
                        end

                        if cont==1  %me tabulate
                            %F.Properties.VariableNames(4)
                            T.ch1=table2array(T(:,4))>=thr1;
                            T.ch2=table2array(T(:,5))>=thr2;
                            T.ch3=table2array(T(:,6))>=thr3;
                            T.ch1ch2=(table2array(T(:,4))>=thr1) & (table2array(T(:,5))>=thr2);
                            T.ch1ch3=(table2array(T(:,4))>=thr1) & (table2array(T(:,6))>=thr3);
                            T.ch2ch3=(table2array(T(:,5))>=thr2) & (table2array(T(:,6))>=thr3);
                            T.ch1ch2ch3=(table2array(T(:,4))>=thr1) & (table2array(T(:,5))>=thr2) & (table2array(T(:,6))>=thr3);
                            T(T.ch1+T.ch2+T.ch3==0,:)=[];
                            sum_ch1=sum(T.ch1);sum_ch2=sum(T.ch2);sum_ch3=sum(T.ch3);
                            sum_ch1ch2=sum(T.ch1ch2);sum_ch1ch3=sum(T.ch1ch3);sum_ch2ch3=sum(T.ch2ch3);
                            sum_ch1ch2ch3=sum(T.ch1ch2ch3);
                            sum_c1= sum(table2array(T((table2array(T(:,4))>=thr1),4)));
                            sum_c2= sum(table2array(T((table2array(T(:,5))>=thr2),5)));
                            sum_c3= sum(table2array(T((table2array(T(:,6))>=thr3),6)));
                            sum_of_cell=length(T.ROI);
                            ALL_val={'ALL Regions',length(regionsx) ,sum_of_cell,sum_c1,sum_c2,sum_c3,sum_ch1,sum_ch2,sum_ch3,sum_ch1ch2,sum_ch1ch3,sum_ch2ch3,sum_ch1ch2ch3};
                            %summary(T);
                            currpath=cd;
                            cd ../ % go up
                            writetable(T,[f_name,'_all_mono','.csv']); % save raw excel table with centroids
                            cd(currpath)
                            columnIndicesToDelete = [1 2 8]; % and so on
                            T(:,columnIndicesToDelete) = [];
                            G=groupsummary(T,"ROI","sum");
                            for aa=1: length(G.ROI)
                                G.ROI_names(aa)=regionsx(G.ROI(aa)); % get names of roi
                            end
                            G(:,1) = [];
                            F = G(:, [13 1 2 3 4 5 6 7 8 9 10 11 12]);

                            F=[F;ALL_val];%final summary table

                            F.sum_ROI=(F.sum_ch1+F.sum_ch2+F.sum_ch3-F.sum_ch1ch2-F.sum_ch1ch3-F.sum_ch2ch3-F.sum_ch1ch2ch3-F.sum_ch1ch2ch3);
                            cd ../ % go up
                            writetable(F,[f_name,'_sum_mono','.csv']); % save raw excel table with centroids
                            %    statit(trial,f_name,thr1,cont,1,1)
                            % save([f_name,'_roimask.mat'],'roi_mask');


                            clear T G F
                        elseif cont==2 % AISH % allen
                            sum_of_area=sum(T.area);
                            %                         sum_blob=sum(sum(zmask>0));
                            %                         ALL_val={'ALL Regions',length(regions) ,sum_of_area,sum_blob};
                            ALL_val={'ALL Regions',length(regions) ,sum_of_area};
                            %summary(T);
                            writetable(T,[f_name,'_all_mono','.csv']); % save raw excel table with centroids
                            columnIndicesToDelete = [1 2  5]; % and so on
                            try
                                T(:,columnIndicesToDelete) = [];
                                G=groupsummary(T,"ROI","sum");
                                %   G.ROI_names=all_names(G.ROI);
                                for aa=1: length(G.ROI)
                                    G.ROI_names(aa)=regionsx(G.ROI(aa)); % get names of roi
                                end
                                %                                 for ap=1:length(apendy)
                                %                                 add_row={apendy(ap),0 ,0,0,0,0,regions(apendy(ap))};
                                %                                 G=[G;add_row];
                                %                                 end
                                G=sortrows(G,1);
                                G(:,1) = [];
                                F = G(:, [3 1 2]);
                                %                                 F=[F;ALL_val];%final summary table
                            catch
                                F=T;
                            end
                            cd(trial)
                            allblobs(apendy)=[];
                            %                            allblobs(allblobs==0)=[];
                            allblobs=allblobs'/0.35;
                            allblobs(end)=[];
                            F(end,:)=[];
                            F.perc=F.GroupCount./allblobs;
                            %                             addtab=[regions(zerid)',zeros(sum(zerid),3)];
                            %                             F=[table2array(F);addtab];% add table to F with regions not found
                            %                             F((length(allblobs)+1:(sum(zerid)+length(allblobs))),2:4)=addtab;
                            %                             F=table(F);
                            F=sortrows(F,4);
                            %                             ALL_val=[ALL_val,100];
                            %                             F=[F;ALL_val];%final summary table
                            tonormal=[tonormal;F.perc];

                            writetable(F,[f_name,'_sum_mono.csv']); % save raw excel table with centroids
                            clear T G F
                        else  %lin tabulate ie 4x info image with not all regions
                            %     sum_ch1= sum(table2array(T((table2array(T(:,4))>thr1),4)));
                            % sum_ch2= sum(table2array(T((table2array(T(:,4))>thr2),5)));
                            % sum_ch3= sum(table2array(T((table2array(T(:,4))>thr3),6)));
                            %   sum_of_cell=length(T.ROI);
                            %   sum_of_GFP=sum(T.GFP);sum_of_RFP=sum(T.RFP);sum_of_both=sum(T.overlap);
                            %   ALL_val={"ALL Regions",length(regions) ,sum_of_cell,sum_ch1,sum_ch2,sum_ch3 };
                            %summary(T);
                            sum_of_area=sum(T.area);
                            sum_of_GFP=sum(T.GFP);sum_of_RFP=sum(T.RFP);sum_of_both=sum(T.overlap);mean_of_value=sum(T.value);
                            ALL_val={'ALL Regions',length(regionsx) ,sum_of_area,sum_of_GFP,sum_of_RFP,sum_of_both,mean_of_value};
                            %summary(T);
                            writetable(T,[f_name,'_all_mono','.csv']); % save raw excel table with centroids
                            columnIndicesToDelete = [1 2 8]; % and so on
                            try
                                T(:,columnIndicesToDelete) = [];
                                G=groupsummary(T,"ROI","sum");
                                %   G.ROI_names=all_names(G.ROI);
                                for aa=1: length(G.ROI)
                                    G.ROI_names(aa)=regionsx(G.ROI(aa)); % get names of roi
                                end
                                for ap=1:length(apendy)
                                    add_row={apendy(ap),0 ,0,0,0,0,regions(apendy(ap))};
                                    G=[G;add_row];
                                end
                                G=sortrows(G,1);
                                G(:,1) = [];
                                F = G(:, [6 1 2 3 4 5]);
                                F=[F;ALL_val]%final summary table
                            catch
                                F=T;
                            end
                            % %   F.sum_ROI=F;
                            %     writetable(F,[f_name,'_sum_mono','.csv']); % save raw excel table with centroids
                            cd(trial)
                            writetable(F,[f_name,'_sum_mono','.csv']); % save raw excel table with centroids
                            % save([f_name,'_roimask.mat'],'roi_mask');
                            % %   cd ../ % go up

                            % %   writetable(F,[f_name,'_sum_mono','.csv']); % save raw excel table with centroids
                            % % writetable(T,[f_name,'_all_mono','.csv']); % save raw excel table with centroids
                            % % columnIndicesToDelete = [1 2 8]; % and so on
                            % % T(:,columnIndicesToDelete) = [];
                            % %  G=groupsummary(T,"ROI","sum");
                            % % for aa=1: length(G.ROI)
                            % %   G.ROI_names(aa)=regions(G.ROI(aa)); % get names of roi
                            % % end
                            % % G(:,1) = [];
                            % % F = G(:, [6 1 2 3 4 5 ]);
                            % %
                            % % F=[F;ALL_val]%final summary table
                            % % writetable(F,[f_name,'_sum_mono','.csv']); % save raw excel table with centroids
                            % % %    statit(trial,f_name,thr1,cont,1,1)

                            clear T G F

                        end
                    end

                end      % all files
                save('tonoraml.mat','tonormal')
            else     % for all areas ==ALL true
                all_idx = 2:1327;% all subregions
                all_names= string(st.acronym);
                all_names(strcmp(all_names, 'root')) = [];
                for w=1:length(text_file_names) % for all the image in th folder
                    cd(trial)
                    w
                    T=readtable(text_file_names(w));
                    if height(T)>0
                        T.ROI= zeros(height(T),1);
                        fchar=char(text_file_names(w));
                        try
                            f_name=fchar(1:(end-33));
                            load([f_name,'_maskcrop','.mat'], 'zmask'); % get mask of rois
                        catch
                            f_name=extractBefore(fchar,'_geo');
                            load([f_name,'_maskcrop','.mat'], 'zmask'); % get mask of rois
                        end
                        %             end
                        % %              try
                        % %                     f_name=fchar(1:(end-37));
                        % % load([f_name,'_maskcrop','.mat'], 'zmask'); % get mask of rois
                        % %             end
                        for k=1:length(all_idx) %  subregions
                            roi_area=all_names(k);
                            roi_id_path = st.structure_id_path(find(strcmp(st.acronym,roi_area)));
                            roi_idx = find(contains(st.structure_id_path,roi_id_path));% all subregions

                            cent_x=floor(table2array(T(:,1)));
                            cent_y=floor(table2array(T(:,2)));
                            dot_roi=zeros(length(cent_x),1);
                            for nn=1:length(cent_x)
                                dot_roi(nn)= zmask(cent_y(nn),cent_x(nn));%??? x and y
                            end
                            for n=1:length(roi_idx)
                                pos=dot_roi==roi_idx(n);
                                T.ROI(pos)=k;
                            end
                        end
                        T(~T.ROI,:)=[];
                        %   T.index = (1:height(T)).'
                        %   T = T(:, [7 1 2 3 4 5 6 7])
                        T.ROI_names=all_names(T.ROI); % get names of roi
                        sum_of_area=sum(T.area);
                        sum_of_GFP=sum(T.GFP);sum_of_RFP=sum(T.RFP);sum_of_both=sum(T.overlap);
                        ALL_val={'ALL Regions',length(regions) ,sum_of_area,sum_of_GFP,sum_of_RFP,sum_of_both};
                        %summary(T);
                        writetable(T,[f_name,'_all_mono','.csv']); % save raw excel table with centroids
                        columnIndicesToDelete = [1 2  8]; % and so on
                        T(:,columnIndicesToDelete) = [];
                        G=groupsummary(T,"ROI","sum");
                        G.ROI_names=all_names(G.ROI);
                        G(:,1) = [];
                        F = G(:, [6 1 2 3 4 5]);
                        F=[F;ALL_val]%final summary table
                        %     writetable(F,[f_name,'_sum_mono','.csv']); % save raw excel table with centroids

                        cd(trial)
                        writetable(F,[f_name,'_sum_mono','.csv']); % save raw excel table with centroids

                        % %   cd ../ % go up

                        % %   writetable(F,[f_name,'_sum_mono','.csv']); % save raw excel table with centroids
                        %pause(3)

                        % statit(trial,f_name,thr1,cont,1,auto)
                        clear T
                        clear G
                        clear F
                    end
                end
            end
        end
    else %visy==1
        % ROI extraction
        disp('Loading Libraries...')
        allen_atlas_path = '/data/Aligment/AP_histology-master/AllenCCF';
        st = loadStructureTree([allen_atlas_path filesep 'structure_tree_safe_2017.csv']);
        disp('extracting ROI...')
        for lis=1:length(listy)  % go over multidirectiores
            clear T G F
            disp('Loading barcodes file..')
            path=listy(lis).name;
            cd(path)
            outs_files = dir([path filesep 'out*']); % get the contents of the image_folder
            outs_files = natsortfiles({outs_files.name});
            outs_files = char(outs_files(end))
            try
                load('10xMtxFile.mat', 'barcodesx')
                disp('Loaded previoius barcodes')
            catch
                mtx_file = [path,'/',outs_files,'/outs/filtered_feature_bc_matrix/matrix.mtx'];
                bc_file = [path,'/',outs_files,'/outs/filtered_feature_bc_matrix/barcodes.tsv'];
                gene_file =[path,'/',outs_files,'/outs/filtered_feature_bc_matrix/features.tsv'];
                [~, ~, barcodesx] = load10xMtxFile(mtx_file,bc_file,gene_file);
            end
            % lm=length(geneid_8_1);
            ln=length(barcodesx);
            str_barcodes=string(barcodesx);
            tissue_positions_list=readtable([path,'/',outs_files,'/outs/spatial/tissue_positions_list.csv']);
            xst=tissue_positions_list.Var5;yst=tissue_positions_list.Var6;
            bar_org=string(tissue_positions_list.Var1);
            bar_xy=cell(ln,3);
            for ib=1:ln
                curr_barcode=str_barcodes(ib);
                bar_ind=find(curr_barcode==bar_org);
                bar_xy(ib,:)={curr_barcode,xst(bar_ind),yst(bar_ind)};
            end
            bar_ar=cell2mat(bar_xy(:,(2:3)));
            mapped_xy=bar_ar(:,1:2);
            T = array2table(mapped_xy,...
                'VariableNames',{'cent_x','cent_y'});
            im_tif = dir([path filesep 'Vis*gray.tif']);%future ; take gray out % get the contents of the image_folder
            im_tif = natsortfiles({im_tif.name});
            im_tif = char(im_tif(end));
            load([im_tif(1:(end-4)),'_maskcrop','.mat'], 'zmask'); % get mask of rois
            dot_roi=zeros(length(mapped_xy),1);
            n_mask=zeros(size(zmask));
            for nn=1:length(dot_roi)
                %get value of zmask given x and y
                dot_roi(nn)= zmask(floor(mapped_xy((nn),1)),floor(mapped_xy((nn),2)));
            end
            if regions(1)~="All"
                %       roi_mask.n_mask=n_mask;      roi_mask.regions=regions;

                for k=1:length(regions) %  subregions
                    roi_area=regions(k)
                    roi_id_path = st.structure_id_path(find(strcmp(st.safe_name,roi_area)));%search full name
                    if isempty(roi_id_path)
                        roi_id_path = st.structure_id_path(find(strcmp(st.acronym,roi_area)));% search acronym
                    end
                    roi_idx = find(contains(st.structure_id_path,roi_id_path));% all subregions
                    [bx,by]=find(ismember(zmask,roi_idx));   % for submask n_mask

                    for n=1:length(roi_idx)
                        pos=dot_roi==roi_idx(n);
                        T.ROI(pos)=k;
                    end
                    %              roi_mask(k).roi=[bx,by];
                end
                T(~T.ROI,:)=[];

                for aa=1: length(T.ROI)
                    T.ROI_names(aa)=regions(T.ROI(aa)); % get names of roi
                end
                writetable(T,[im_tif(1:(end-4)),'_all_mono','.csv']); % save raw excel table with centroids
                try
                    T(:,4) = [];
                    G=groupsummary(T,"ROI","sum");
                    %   G.ROI_names=all_names(G.ROI);
                    for aa=1: length(G.ROI)
                        G.ROI_names(aa)=regions(G.ROI(aa)); % get names of roi
                    end
                    G(:,[1 3 4]) = [];
                catch
                    F=T;
                end

                writetable(G,[im_tif(1:(end-4)),'_sum_mono','.csv']); % save raw excel table with centroids
                % save([im_tif(1:(end-4)),'_roimask.mat'],'roi_mask');

            else     % for all areas ==ALL true
                all_idx = 2:1327;% all subregions
                all_names= string(st.acronym);
                % all_names(strcmp(all_names, 'root')) = [];

                for k=1:length(all_idx) %  subregions
                    roi_area=all_names(k);
                    roi_id_path = st.structure_id_path(find(strcmp(st.acronym,roi_area)));
                    roi_idx = find(contains(st.structure_id_path,roi_id_path));% all subregions
                    for nn=1:length(dot_roi)
                        %get value of zmask given x and y
                        dot_roi(nn)= zmask(floor(mapped_xy((nn),1)),floor(mapped_xy((nn),2)));
                    end
                    for n=1:length(roi_idx)
                        pos=dot_roi==roi_idx(n);
                        T.ROI(pos)=k;
                    end
                end
                T(~T.ROI,:)=[];
                T.ROI_names=all_names(T.ROI); % get names of roi
                writetable(T,[im_tif(1:(end-4)),'_all_mono','.csv']); % save raw excel table with centroids
                T(:,4) = [];
                G=groupsummary(T,"ROI","sum");
                G.ROI_names=all_names(G.ROI);
                G(:,[1 3 4]) = []
                writetable(G,[im_tif(1:(end-4)),'_sum_mono','.csv']); % save raw excel table with centroids

            end
        end
    end
elseif state=="Adjust ROI"
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cd('/bigdata/microscope_images/presentation_demo/rabies/LIN')
    load('B4 002_4x_maskline.mat', 'maskline')
    J=imadjust(imread('B4 002_4x_geo.tif','index',1));
    load('B4 002_4x_maskcrop.mat', 'zmask')
    %% mask2poly
    nmask=zmask>1;
    im=zmask>1;
    CC=imfuse(J,maskline,'blend');

    % Get the contour
    im=imfill(im>0,4,'holes');
    C=bwboundaries(im);
    C=C{1};
    % Visualize binary image and the contour
    fprintf('Decimate by a factor of 500\n')
    disp('please wait')
    my_vertices=flip(DecimatePoly(C,[0.002 2]),2);
    fpoly=figure;
    imshow(CC);
    polyx = drawpolygon('Position',my_vertices);
    axis equal;  axis off;
    addlistener(polyx,'MovingROI',@allevents);
    addlistener(polyx,'ROIMoved',@allevents);
    uiwait(fpoly);
    %transform
    [mp,fp] = cpselect(maskline,J,CurrPos,my_vertices,'Wait',true);
    t = fitgeotrans(mp,fp,'polynomial',4);
    Rfixed = imref2d(size(J));
    masklinex = imwarp(maskline,t,'nearest','OutputView',Rfixed);
    disp('Fusing mask and subimage..')
    imshowpair(J,masklinex,'blend');


elseif state=="Edit ROI"
    editor2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55

elseif state=="SpreadSheet" % show files agains gfp rfp starter per roi
    for lis=1:length(listy)
        dirx=listy.name;
        cd(dirx);
        diry=strsplit(dirx,'/');
        diry=char(diry(end));
        batchy=dir(fullfile([diry,'*all_batch.csv']));

        for rg=1:length(regions)
            T = readtable(batchy.name);
            roinames=string(T.ROI_names);

            % for kk=1:height(T) % change names
            %  yvalues{1,1}{kk,1}=[num2str(kk),yvalues{1,1}{kk,1}];
            % end
            currt=T(roinames==num2str(regions(rg)),:);
            if ~isempty(currt)
                % read the data into  a variable
                % % % % V = T.Properties.VariableNames; % take out irrilavant
                % % % % for i = [1:width(T)]
                % % % %     v_is_cell(i) = iscell(T.(V{i}));
                % % % % end
                %use logical indexing to delete the required columns
                % T(:,v_is_cell) = [];
                tconvert=table2array(currt(:,[5 6 7]));
                [files,ndx]=natsortfiles(currt.file_name);% sort them naturally
                files=string(files);
                files=  extractBefore(files,"_sum"); % extract in any name before _sum!!!
                tconvert=tconvert(ndx,:)';
                % sort them naturally
                % figure('Name',num2str(dirx),'NumberTitle','off');
                figure;
                h = heatmap(files,{'GFP','RFP','Both'},tconvert);
                h.ColorScaling = 'scaledrows';
                title(num2str(regions(rg)))
                clear T
                clear currt
            else
                regions(rg)
                disp('region was not found')
            end
        end
    end
    %%%%%%
elseif  state=="Crop ROI"
    if visy==0
        fprintf('Please Select '); fprintf(2, 'Annotated Image/s \n');
        listy=uipickfiles('num',[],'FilterSpec',['*',suffix,'tif'],'out','struct'); % loop images
        clc
        for i=1:length(listy)
            cd(listy(i).folder);
            im_tif = string(listy(i).name);
            diry=strsplit(im_tif,'/');
            diry=char(diry(end));
            f_name=diry(1:(end-8));
            load([f_name,'_maskcrop','.mat']); % get mask of rois
            im1=imread(diry,'index',1);C=imadjust(im1);
            bfr= BioformatsImage(diry);
            im2=imread(diry,'index',2); im3=imread(diry,'index',3);
            if bfr.sizeT>3  || bfr.sizeZ>3 || bfr.sizeC>3
                im4=imread(diry,'index',4);
            end
            bfr= BioformatsImage(diry);
            xmask=imadjust(zmask);
            J = imfuse(C,xmask,'blend');
            cropin = input('CROP IN(y)?','s');
            clc
            if cropin=="y"
                f1=figure('Name',num2str(diry),'NumberTitle','off');
                imshow(J);imcontrast; [h, w, ~] = size(J);
                imgzoompan(gcf, 'ImgWidth', w, 'ImgHeight', h);
                h1 = drawassisted('label',"Select"); bw1= createMask(h1);title('Please Select accurate ROI, then CLOSE image');
                uiwait(f1); im1(~bw1)=0;  im2(~bw1)=0; im3(~bw1)=0;J(~bw1)=0;
                if bfr.sizeT>3  || bfr.sizeZ>3 || bfr.sizeC>3
                    im4(~bw1)=1;
                end
                zmask(~bw1)=1;
            end

            cropout = input('CROP OUT(y)?','s');
            clc
            if cropout=="y"
                f2=figure('Name',num2str(diry),'NumberTitle','off');
                imshow(J); imcontrast; [h, w, ~] = size(J);
                imgzoompan(gcf, 'ImgWidth', w, 'ImgHeight', h);
                title('Please de-Select accurate ROI, then CLOSE image');
                h2= drawassisted('Color','r','label',"deSelect");  bw2= createMask(h2);uiwait(f2);
                im1(bw2)=0; im2(bw2)=0;  im3(bw2)=0;
                if bfr.sizeT>3  || bfr.sizeZ>3 || bfr.sizeC>3
                    im4(bw2)=0;
                end
                zmask(bw2)=1;
            end
            if cropout=="y" || cropin=="y"
                prompt = {'Slice Name (keep format: geo.tif)'};
                dlgtitle = 'Name';
                dims = [1 35];
                definput = {diry};
                answer = char(inputdlg(prompt,dlgtitle,dims,definput));
                if isempty(answer)==1
                    close all
                else
                    imwrite(im1,answer);
                    disp('Edit Channel 1');

                    disp('Edit Channel 2');
                    imwrite(im2,answer,'WriteMode','append');

                    imwrite(im3,answer,'WriteMode','append');
                    disp('Edit Channel 3');
                    if bfr.sizeT>3  || bfr.sizeZ>3 || bfr.sizeC>3
                        im4(bw2)=0;
                        imwrite(im4,answer,'WriteMode','append');
                        disp('Edit Channel 4');
                    end
                    disp('saving mask')
                    save([answer(1:(end-4)),'_maskcrop.mat'],'zmask','masklined','RGBmask');
                end
            end
        end
    else % visy==1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5555
        fprintf('Please Select '); fprintf(2, 'Annotated Image/s \n');
        listy=uipickfiles('num',[],'FilterSpec',['*',suffix,'tif'],'out','struct'); % loop images
        clc
        for i=1:length(listy)
            path=listy(i).name;cd(path);
            im_tif = dir([path filesep 'Vis*gray.tif']);%future ; take gray out % get the contents of the image_folder
            im_tif = natsortfiles({im_tif.name});
            im_tif = char(im_tif(end));
            load([im_tif(1:(end-4)),'_maskcrop','.mat'], 'zmask'); % get mask of rois
            imgv=imread([im_tif(1:(end-4)),'.tif']);
            f1=figure('Name',im_tif(1:(end-4)),'NumberTitle','off');
            imshow(imgv);
            h1 = drawassisted('label',"Select"); bw1= createMask(h1);title('Please Select accurate ROI, then CLOSE image');
            uiwait(f1); imgv(~bw1)=0;
            zmask(~bw1)=1;

            prompt = {'Slice Name (keep format: geo.tif)'}; % break point
            dlgtitle = 'Name';
            dims = [1 35];
            definput = {[im_tif(1:(end-4)),'.tif']};
            answer = char(inputdlg(prompt,dlgtitle,dims,definput));
            if isempty(answer)==1
                close all
            else
                imwrite(imgv,[answer]);
                disp('Edit Channel 1');
                disp('saving mask')
                save([answer(1:(end-4)),'_maskcrop','.mat'], 'zmask','masklined','RGBmask'); % get mask of rois
                %  make plot mask lines and save them;
                mask_val=unique(zmask);
                total_areas=length(mask_val);
                % figure; imagesc(zmask);hold on
                pb = CmdLineProgressBar('Finding updated cooredniates for each region...');
                for maskr=2:total_areas
                    pb.print(maskr,total_areas)
                    BW =zmask==mask_val(maskr);
                    [B,~] = bwboundaries(BW,'noholes');
                    for kkk = 1:length(B)
                        boundaryx = B{kkk};
                        hold on
                        plot(boundaryx(:,2),boundaryx(:,1), 'r', 'LineWidth', 2)
                        hold on
                        k_bound(kkk)={boundaryx};
                    end

                    maskline(maskr-1)={k_bound};
                    clear k_boundend
                end
                clear maskline
            end
            set ( gca, 'xdir', 'reverse' )
            set ( gca, 'ydir', 'reverse' )
            axis equal;axis tight;axis off;
            exportgraphics(gca,'masklines.png');
        end % list folders
    end     % visium
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Extract ABCD
elseif state=="Extract ABCDE"
readlines
end % swithes
warning('on','all')
disp('Finished');

%% function for updating polygon movement
% function allevents(~,evt)
% evname = evt.EventName;
% switch(evname)
%     case{'MovingROI'}
%         % %             disp(['ROI moving previous position: ' mat2str(evt.PreviousPosition)]);
%         % %             disp(['ROI moving current position: ' mat2str(evt.CurrentPosition)]);
%     case{'ROIMoved'}
%         %             disp(['ROI moved previous position: ' mat2str(evt.PreviousPosition)]);
%         %             disp(['ROI moved current position: ' mat2str(evt.CurrentPosition)]);
%         assignin('base','CurrPos',evt.CurrentPosition);
% end
% end
