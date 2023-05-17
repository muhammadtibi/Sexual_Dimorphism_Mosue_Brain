%get statistics for all
%%
% cd('/bigdata/microscope_images/Etay/Amygdala/Brains/Males/EA10_1_BMA/Tiff files/10_1_L')
clear all
warning('off','all')
uuu=cd;
transcendentalDoc = '';
clc
construct=0;
try
    S = {...
        struct('name','Stat. Level','type','enum','values',{{'Batch1','Batch2','Versus','Mono'}},'doc',transcendentalDoc);...
        struct('name','EXP type Name','type','enum','values',{{'4X','20x','AISH','VIS'}});...
        struct('name','suffix (w|o _sum/all.csv)','type','str','default','.');...%_L_20x
        struct('name','Magnification','type','int','default',4);...
        struct('name','preffix','type','str','default','.');...%_L_20x
        
        % % struct('name','Condition name/s','type','buttonGroup','groupCount',[0 2],...
        % %     'groupFields',...
        % %     {{
        % %         struct('name','Condition name','type','str')...
        % %     }},...
        % %      'default',...
        % %     {{...
        % %         {".",3.141592},...
        % %     }}...
        % %     );...
        % % % struct('name','Condition 1 color','type','colorPicker','default',[0 0 1]);...
        % % % struct('name','Condition 2 color','type','colorPicker','default',[1 0 0]);...
        struct('name','*Counter Threshold','type','str','default','3,3,3');...
        % % % struct('name','Repeats?','type','int','default',1);...
        % % % struct('name','Use DROI','type','checkbox');...
        };
    
    Param = Settings_GUI(S);
    state=string(cell2mat(Param(1)));
    exp=string(cell2mat(Param(2)));
    if exp=="20x"
        cont=1;
    elseif exp=="AISH"
        cont=2;
    elseif exp=="VIS"
        cont=3;
    else% 4X
        cont=0;
    end
    suffix=char(Param(3));
    if suffix=="."
        suffix='';
    end
    mag=cell2mat(Param(4));
    if mag==10
        mag=0.65;
    elseif mag==20
        mag=0.32;
    elseif mag==40
        mag=0.16;
    else  %4
        mag=1.63;
    end
    preffix=char(Param(5));
    if preffix=="."
        preffix='';
    end
    % % % condy=Param(5);condy=string(condy{1,1});preffix=condy(1);
    % % % if preffix=="."
    % % %     preffix='';
    % % % end
    % % % preffix=char(preffix);
    % % % try
    % % %     cond2=condy(2);
    % % % if cond2=="."
    % % %     cond2="";
    % % % end
    % % % cond2=char(cond2);
    % % % end
    % % col1=cell2mat(Param(6));
    % % col2=cell2mat(Param(7));
    thr=Param(6);
    str = regexprep(thr,',',' ');
    thr = str2num(cell2mat(str));thr1=thr(1);thr2=thr(2);thr3=thr(3);
    
    
    % % rep=Param(8);
    % % draw=cell2mat(Param(9));
    % % if draw==1
    % %     auto=0;
    % % else
    auto=1;
    % %
    % % end
    
catch ME
    close all
    cd(uuu);
    return;
end


%% switch state

if state=="Mono"
    fprintf('Please Select '); fprintf(2, 'Iner *_sum_mono.csv file/s \n');
    listy=uipickfiles('num',[],'FilterSpec',[preffix,'*',suffix,'_sum_mono.csv'],'out','struct'); % loop images
    for lis=1:length(listy)
        %       wheres= input('aligment Folders name: ','s')
        trial=listy(lis).folder;
        cd(trial)
        
        im_tif = string(listy(lis).name);%update for every image
        diry=strsplit(im_tif,'/');
        diry=char(diry(end));
        diry=char(diry(1:(end-13)));
        statit%%%%%%%%%%%%%%here stat
    end
elseif  state=="Batch1"
    fprintf('Please Select '); fprintf(2, ' a MOTHER Folder/s \n');
    listy=uipickfiles('num',[],'FilterSpec',[preffix,'*',suffix,'_sum_mono.csv'],'out','struct'); % loop images
    for lis=1:length(listy)
        trial=listy(lis).name;
        cd(trial)
        x_file = string(listy(lis).name);
        diry=strsplit(x_file,'/');
        diry=char(diry(end));
        Allfile=dir(fullfile([preffix,'*',suffix,'_sum_mono.csv']));
        StrAllfile=struct2cell(Allfile);
        Filename=StrAllfile(1,:);
        [~,nFile]=size(Filename);
        Allfilex=dir(fullfile([preffix,'*',suffix,'_all_mono.csv']));
        StrAllfilex=struct2cell(Allfilex);
        Filenamex=StrAllfilex(1,:);
        all={};allx={};
        for icsv=1:1:nFile
            T_curr=readtable(Filename{icsv});
            T_currx=readtable(Filenamex{icsv});
                       if height(T_curr)>0
                T_curr.file_name(:)=Filename(icsv);
                all=[all;T_curr];
                T_currx.file_name(:)=Filenamex(icsv);
                allx=[allx;T_currx];
            end
        end
        currpath=cd;
        cd ../
        writetable(allx,[diry,preffix,suffix,'_all_cells.csv']);
        cd(currpath);
        if cont~=1
            all = all(:, [7 1 2 3 4 5 6]);
            writetable(all,[diry,preffix,suffix,'_all_batch.csv']);
            all(:,1) = [];
            % except_all=gcsv(table2array(gcsv.ROI_names)=='All Regions');
            sumy=grpstats(all,'ROI_names','sum');
            %  fcsv=groupsummary(fcsv,"ROI","sum");
            % r=zeros(length(fcsv.ROI_names),1);
            for q=1:length(sumy.ROI_names)
                if   (string(sumy.ROI_names{q}) == "ALL Regions")% find alls and oreder them
                    break;
                end
            end
            ALL=sumy(q,:);
            sumy(q,:)=[];
            x_areas=height(sumy);
            xreport= summary(sumy);
            mean_gfp=mean(sumy.sum_GFP);mean_rfp=mean(sumy.sum_RFP);
            mean_overlap=mean(sumy.sum_overlap);mean_value=mean(sumy.sum_value);
            sumy=[sumy;ALL];
            batchit%%%%%%%%%%%%%% here stat
            cd ..
            sumy.prc_gfp=prcg;sumy.prc_rfp=prcr;sumy.prc_both=prcs;
            writetable(sumy,[diry,'_',preffix,'_',suffix,'_sum_batch.csv']); % save raw excel table with centroids
        else %cont==1
            all = all(:, [15 1 2 3 4 5 6 7 8 9 10 11 12 13 14]);
            allxy=all;
            writetable(all,[diry,'_',preffix,'_',suffix,'_all_batch.csv']);
            all(:,[1 3 ]) = [];
            % except_all=gcsv(table2array(gcsv.ROI_names)=='All Regions');
            sumy=grpstats(all,'ROI_names','sum');
            %  fcsv=groupsummary(fcsv,"ROI","sum");
            % r=zeros(length(fcsv.ROI_names),1);
            for q=1:length(sumy.ROI_names)%loop till you find all regions
                if   (string(sumy.ROI_names{q}) == "ALL Regions")
                    break;
                end
            end
            ALL=sumy(q,:);
            sumy(q,:)=[];%exclude all regions
            x_areas=height(sumy);
            xreport= summary(sumy);
            mean_gfp=mean(table2array(sumy(:,4)));mean_rfp=mean(table2array(sumy(:,5)));mean_cy5=mean(table2array(sumy(:,6)));
            sumy=[sumy;ALL];
            batchit%%%%%%%%%%%%%% here stat
            cd ..
            writetable(sumy,[diry,'_',preffix,'_',suffix,'_sum_batch.csv']); % save raw excel table with centroids
            end
       
        
      
        clear sumy
        clear gcsv
        clear T_curr
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Batch2
elseif state=="Batch2"
    fprintf('Please Select '); fprintf(2, ' a GrandMother Folder/s \n');
    listy=uipickfiles('num',[],'FilterSpec',[preffix,'*',suffix,'_sum_batch.csv'],'out','struct'); % loop images
    clc
        cd(listy(1).name)
        x_file = string(listy(1).name);
        diry=strsplit(x_file,'/');
        diry=char(diry(end));
        Allfile=dir(fullfile([preffix,'*',suffix,'_sum_batch.csv']));
        StrAllfile=struct2cell(Allfile);
        Filename=StrAllfile(1,:);
        [~,nFile]=size(Filename);
        all={};
        for icsv=1:1:nFile
            T_curr=readtable(Filename{icsv});
            if height(T_curr)>0
                T_curr.file_name(:)=Filename(icsv);
                all=[all;T_curr];
            end
        end
        if cont~=1
            all = all(:, [10 1 2 3 4 5 6 7 8 9]);
            writetable(all,[diry,'_',preffix,'_',suffix,'_all_batch2.csv']);
            allx=all;
            brainames=unique(table2array(allx(:,1)));
            all(:,[1 3 ]) = [];
            % except_all=gcsv(table2array(gcsv.ROI_names)=='All Regions');
            sumy=grpstats(all,'ROI_names','sum');
            meany=grpstats(all,'ROI_names','mean');
%             miny=grpstats(all,'ROI_names','min');
%             maxy=grpstats(all,'ROI_names','max');
%             stdy=grpstats(all,'ROI_names','std');
%             mediany=grpstats(all,'ROI_names','median');
            %  fcsv=groupsummary(fcsv,"ROI","sum");
            % r=zeros(length(fcsv.ROI_names),1);
            for q=1:length(sumy.ROI_names)
                if   (string(sumy.ROI_names{q}) == "ALL Regions")% find alls and oreder them
                    break;
                end
            end
            ALL=sumy(q,:);
            sumy(q,:)=[];
            namcell=sumy.ROI_names';
            gfp_all=zeros(nFile,length(namcell));
            rfp_all=zeros(nFile,length(namcell));
            overlap_all=zeros(nFile,length(namcell));
            all_prcg=zeros(nFile,length(namcell));
            all_prcr=zeros(nFile,length(namcell));
            all_prcs=zeros(nFile,length(namcell));
            for qq=1:length(namcell)
%             gfpx=[0,0,0];rfpx=[0,0,0];
%             overlapx=[0,0,0];
            roid=(string(all.ROI_names) == string(namcell(qq)));% find scatters data for each region
            gfpx=all.sum_sum_GFP(roid);
            rfpx=all.sum_sum_RFP(roid);
            overlapx=all.sum_sum_overlap(roid);
            prcgx=all.prc_gfp(roid);
            prcrx=all.prc_rfp(roid);
            prcsx=all.prc_both(roid);
            gfp_all(1:sum(roid),qq)=gfpx;
            rfp_all(1:sum(roid),qq)=rfpx;
            overlap_all(1:sum(roid),qq)=overlapx;
            all_prcg(1:sum(roid),qq)=prcgx;
            all_prcr(1:sum(roid),qq)=prcrx;
            all_prcs(1:sum(roid),qq)=prcsx;
            end
            datax={gfp_all,rfp_all,overlap_all};
            x_areas=height(sumy);
            xreport= summary(sumy);
            batch2;
%             x_areas=height(sumy);
%             mean_gfp=mean(sumy.sum_sum_GFP);mean_rfp=mean(sumy.sum_sum_RFP);mean_overlap=mean(sumy.sum_sum_overlap);
%             sumy=[sumy;ALL];
            
        else %cont==1
            all = all(:, [15 1 2 3 4 5 6 7 8 9 10 11 12 13 14]);
            writetable(all,[diry,'_',preffix,'_',suffix,'_all_batch2.csv']);
            all(:,[1 3 ]) = [];
            % except_all=gcsv(table2array(gcsv.ROI_names)=='All Regions');
            sumy=grpstats(all,'ROI_names','sum');
            %  fcsv=groupsummary(fcsv,"ROI","sum");
            % r=zeros(length(fcsv.ROI_names),1);
            for q=1:length(sumy.ROI_names)%loop till you find all regions
                if   (string(sumy.ROI_names{q}) == "ALL Regions")
                    break;
                end
            end
            ALL=sumy(q,:);
            sumy(q,:)=[];%exclude all regions
            x_areas=height(sumy);
            xreport= summary(sumy);
            mean_gfp=mean(table2array(sumy(:,4)));mean_rfp=mean(table2array(sumy(:,5)));mean_cy5=mean(table2array(sumy(:,6)));
            sumy=[sumy;ALL];
        end
%         cd ..
           sumy(:,[7 8 9])=[];
           sumy.sum_prc_gfp=table2array(meany(:,7));sumy.sum_prc_rfp=table2array(meany(:,8));sumy.sum_prc_both=table2array(meany(:,9));
           sumy=[sumy;ALL];
           writetable(sumy,[diry,'_',preffix,'_',suffix,'_sum_batch2.csv']); % save raw excel table with centroids
           writetable(meany,[diry,'_',preffix,'_',suffix,'_mean_batch2.csv']); % save raw excel table with centroids

%         batchit%%%%%%%%%%%%%% here stat
        clear meany 
        clear sumy
        clear gcsv
        clear T_curr
      
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% VS
else     %"versus"
    siglvl=input('type significance level:');
    fprintf('Please Select '); fprintf(2, ' TWO MOTHER Folders \n');
    listy=uipickfiles('num',[],'FilterSpec',[preffix,'*',suffix,'_sum_mono.csv'],'out','struct'); % loop images
    clc
    % 1
    versusit
end % state switch 
warning('on','all')

clc
