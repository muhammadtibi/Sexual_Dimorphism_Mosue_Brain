clear all
close all
addpath('/data/matlab_functions/');
 %%
 % read single sample and show ditribution of starter in the different
 % areas of the brain

 

foldername_M = '//bigdata/microscope_images/Etay/Amygdala/Analysis/Starter_distribution/Males';
foldername_F = '//bigdata/microscope_images/Etay/Amygdala/Analysis/Starter_distribution/Females';
fname_M = dir([foldername_M,'/*__*']);
fname_M = {fname_M.name}'
fname_F = dir([foldername_F,'/*__*']);
fname_F = {fname_F.name}'

cd(foldername_M);
file_cell_M = cell(length(fname_M),1);

for i=1:length(fname_M)
    t = readtable(fname_M{i});
    hline = t.Properties.VariableNames;
    file_cell_M{i} = table2cell(t);
end


cd(foldername_F);
file_cell_F = cell(length(fname_F),1);

for i=1:length(fname_F)
    t = readtable(fname_F{i});
    hline = t.Properties.VariableNames;
    file_cell_F{i} = table2cell(t);
end

% Brain region color from allen
foldername_structure_tree = '/bigdata/microscope_images/Etay/Amygdala/Analysis/Anotation';
cd (foldername_structure_tree);
t = readtable('structure_tree_safe_2017.csv');
hline_structure_tree = t.Properties.VariableNames;
file_cell_structure_tree = table2cell(t);

%
brain_regions = {'others'};
for i=1:length(file_cell_M)-1
tmp = unique([file_cell_M{i}(:,1); file_cell_M{i+1}(:,1)],'stable');
brain_regions= unique([brain_regions; tmp],'stable');
end




%% delete "fiber tracts" "others" and "all region from data"

exlude_regions = [{'fiber tracts'} {'others'} {'ALL Regions'}];

loc = ismember(brain_regions,exlude_regions);
brain_regions(loc)=[];

for i=1:length(file_cell_M)
    loc = ismember(file_cell_M{i}(:,1),exlude_regions);
    file_cell_M{i}(loc,:)=[];
end

for i=1:length(file_cell_F)
    loc = ismember(file_cell_F{i}(:,1),exlude_regions);
    file_cell_F{i}(loc,:)=[];
end

% set brain region in allen order

[~,loc]=ismember(brain_regions,file_cell_structure_tree(:,4));
brain_regions= file_cell_structure_tree(sort(loc),4);

%% adjust areas for all amygdala
AMY_sub_regions = [{'sAMY'} {'BMA'} {'PAA'} {'BLA'} {'LA'} {'COA'} {'PA'}];
AMY_meta_region = [{'Cortical Amy'};{'Subplate Amy'};{'Striatal Amy'}];
% starter distribution within the Amy

for i = 1:length(file_cell_M)
    loc1=(find(ismember(file_cell_M{i}(:,1),AMY_sub_regions([3,6]))));
    loc2=(find(ismember(file_cell_M{i}(:,1),AMY_sub_regions([2,4,5,7]))));
    loc3=(find(ismember(file_cell_M{i}(:,1),AMY_sub_regions(1))));
    AMY_file_cell_M{i,1}(1,1) = sum(cell2mat(file_cell_M{i}(loc1,6)));
    AMY_file_cell_M{i,1}(2,1) = sum(cell2mat(file_cell_M{i}(loc2,6)));
    AMY_file_cell_M{i,1}(3,1) = sum(cell2mat(file_cell_M{i}(loc3,6)));
    tmp1= sum(AMY_file_cell_M{i});
    AMY_file_cell_M{i}(:,2) = AMY_file_cell_M{i}(:,1)./tmp1;
end
for i = 1:length(file_cell_F)
    loc1=(find(ismember(file_cell_F{i}(:,1),AMY_sub_regions([3,6]))));
    loc2=(find(ismember(file_cell_F{i}(:,1),AMY_sub_regions([2,4,5,7]))));
    loc3=(find(ismember(file_cell_F{i}(:,1),AMY_sub_regions(1))));
    AMY_file_cell_F{i,1}(1,1) = sum(cell2mat(file_cell_F{i}(loc1,6)));
    AMY_file_cell_F{i,1}(2,1) = sum(cell2mat(file_cell_F{i}(loc2,6)));
    AMY_file_cell_F{i,1}(3,1) = sum(cell2mat(file_cell_F{i}(loc3,6)));
    tmp1= sum(AMY_file_cell_F{i});
    AMY_file_cell_F{i}(:,2) = AMY_file_cell_F{i}(:,1)./tmp1;

end
%% set regions order and add % in column 11
loc=(find(ismember(brain_regions,AMY_sub_regions)));
brain_regions(loc,:)=[];
brain_regions(2:end+1,1)=brain_regions;
brain_regions(1,1)={'MEA'};
    
for i=1:length(file_cell_M)
    loc=(find(ismember(file_cell_M{i}(:,1),AMY_sub_regions)));
    AMY_starter = sum(cell2mat(file_cell_M{i}(loc,6)));
    file_cell_M{i}(loc,:)=[];
    file_cell_M{i}(2:end+1,:)=file_cell_M{i};
    file_cell_M{i}(1,1)={'MEA'};
    file_cell_M{i}(1,6)={AMY_starter};
    tmp= sum(cell2mat(file_cell_M{i}(:,6)));
    file_cell_M{i}(:,11) = m2c(cell2mat(file_cell_M{i}(:,6))./tmp);
end
for i=1:length(file_cell_F)
    loc=(find(ismember(file_cell_F{i}(:,1),AMY_sub_regions)));
    AMY_starter = sum(cell2mat(file_cell_F{i}(loc,6)));
    file_cell_F{i}(loc,:)=[];
    file_cell_F{i}(2:end+1,:)=file_cell_F{i};
    file_cell_F{i}(1,1)= {'MEA'};
    file_cell_F{i}(1,6)={AMY_starter};
    tmp= sum(cell2mat(file_cell_F{i}(:,6)));
    file_cell_F{i}(:,11) = m2c(cell2mat(file_cell_F{i}(:,6))./tmp);
end

%% make starter % matrix

starter_perc_M = zeros(length(brain_regions),length(file_cell_M));
for i=1:length(file_cell_M)
    for k=1:length(brain_regions)
        loc=find(ismember(file_cell_M{i}(:,1),brain_regions(k)));
        if loc~=0
            starter_perc_M(k,i) = cell2mat(file_cell_M{i}(loc,11));
        end
    end
end


starter_perc_F = zeros(length(brain_regions),length(file_cell_F));
for i=1:length (file_cell_F)
    for k=1:length(brain_regions)
        loc=find(ismember(file_cell_F{i}(:,1),brain_regions(k)));
        if loc~=0
            starter_perc_F(k,i) = cell2mat(file_cell_F{i}(loc,11));
        end
    end
end

%% Make amygdala matrix % divide AMY to CTXpl CTXsp and Striatum



AMY_starter_perc_M = zeros(length(AMY_meta_region),length(AMY_file_cell_M));
for i=1:length(AMY_file_cell_M)
    AMY_starter_perc_M(:,i) = AMY_file_cell_M{i}(:,2);
end
AMY_starter_perc_F = zeros(length(AMY_meta_region),length(AMY_file_cell_F));
for i=1:length(AMY_file_cell_F)
    AMY_starter_perc_F(:,i) = AMY_file_cell_F{i}(:,2);
end

%%

brain_region_color = {};

for k=1:length(brain_regions)
    for i=1:length(file_cell_structure_tree)
        if strcmp(brain_regions(k),file_cell_structure_tree(i,4))==1
            brain_region_color(k)=file_cell_structure_tree(i,14);
        end
    end
end
for i=1:length(brain_region_color)
    if isempty(brain_region_color{1,i})
        brain_region_color(i)= {'AAAAAA'};
    end
end
str='#';
brain_region_color= append(str,brain_region_color);
brain_region_color = string(brain_region_color);
w=strcmp("#19399",brain_region_color(1,:));
for i=1:length(w)
    if w(i)==1
        brain_region_color(i)= '#AAAAAA';
    end
end

brain_region_color(1) = string({'#497689'});
%%
AMY_meta_color (1,1:3) = [{'#70FF70'};{'#8ADA87'};{'#98D6F9'}];
AMY_meta_color = string(AMY_meta_color);



%% 
% AVG of each sex
AVG_starter_perc_M = mean(starter_perc_M');
AVG_starter_perc_F = mean(starter_perc_F');
AVG_AMY_starter_perc_M = mean(AMY_starter_perc_M');
AVG_AMY_starter_perc_F = mean(AMY_starter_perc_F');

%% figures
figure ('position',[100,400,400,800],'color','w');
subplot(2,2,1);
h_m = bar(1:length(file_cell_M), starter_perc_M(1:end,:),1,'stack');
for i = 1 : length(brain_region_color)
    h_m(i).FaceColor = brain_region_color(i);
end
ylim([0 1])
xtickangle(90)
axis tight
set(gca,'xticklabel',fname_M)
title('Starter distribution Males');
xlabel('sampels');
ylabel('starter %');
legend('Amy','Isocortex','OLF', 'HPF','CTXsp','STR','PAL','TH','HY','MB','HB');
subplot(2,2,2);
h_f = bar(1:length(file_cell_F), starter_perc_F(1:end,:), 1, 'stack');
for i = 1 : length(brain_region_color)
    h_f(i).FaceColor = brain_region_color(i);
end
ylim([0 1])
xtickangle(90)
axis tight
set(gca,'xticklabel',fname_F)
title('Starter distribution Females');
xlabel('sampels');
ylabel('starter %');

%
rgb=hex2rgb(cell2mat(brain_region_color'));

% Pie of AVG
subplot(2,2,3);
explode = [1 0 0 0 0 0 0 0 0 0 0];
ax = gca();
p_M=pie(ax,AVG_starter_perc_M,explode);
colormap(rgb)
% ax.Colormap=brain_region_color;
subplot(2,2,4);
ax = gca();
P_F=pie(ax,AVG_starter_perc_F,explode);
colormap(rgb)

figure ('position',[500,800,600,600],'color','w');
AMY_rgb = hex2rgb(cell2mat(AMY_meta_color'));
%Pie of AVG AMY
subplot(3,2,5);
ax = gca();
p_M=pie(ax,AVG_AMY_starter_perc_M);
colormap(AMY_rgb)
% ax.Colormap=brain_region_color;
subplot(3,2,6);
ax = gca();
P_F=pie(ax,AVG_AMY_starter_perc_F);
colormap(AMY_rgb)
legend(AMY_meta_region);

