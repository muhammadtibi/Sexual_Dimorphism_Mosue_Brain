function fillOutUIControls( paramarray, paramValues,paramEditHandles )
% AUTHOR:   Abel Szkalisity
% DATE:     Dec 09, 2016
% NAME:     fillOutUIControls
%
% In cooperation with generateUIControls this function fills up UI fields
% with parameters. (e.g. this can be used to see the current settings of
% your method). This functionality can be replaced during generation with
% the default fields of generateUIControls.
%
% INPUT:
%   paramarray      see generateUIControls function same input.
%   paramValues     A cellarray (exactly as long as paramarray and the
%                   matching positions refers to the same parameter) The
%                   type of the parameter value in this array is determined
%                   by the paramarray matching position. See
%                   fetchUIControlValues paramcell output.
%   paramEditHandles See the output of generateUIControls.
%
% COPYRIGHT
% Settings Template Toolbox. All Rights Reversed. 
% Copyright (C) 2019 Abel Szkalisity
% BIOMAG, Synthetic and System Biology Unit, Institute of Biochemsitry,
% Biological Research Center, Szeged, Hungary
% Ikonen group Department of Anatomy, Faculty of Medicine, University of
% Helsinki, Helsinki, Finland.

for i=1:length(paramarray)
    switch paramarray{i}.type
        case 'int'
            set(paramEditHandles{i},'string',num2str(paramValues{i}));
        case 'str'
            set(paramEditHandles{i},'string',paramValues{i});
        case 'enum'                        
            for j=1:length(paramarray{i}.values)
                %TODO: handle cases when the enum is not string.
                if strcmp(paramValues{i},paramarray{i}.values{j})
                    break;
                end                    
            end
            set(paramEditHandles{i},'value',j);
        case 'checkbox'
            set(paramEditHandles{i},'value',paramValues{i});
        case 'slider'
            set(paramEditHandles{i},'value',paramValues{i});
        case 'colorPicker'
            set(paramEditHandles{i},'BackgroundColor',paramValues{i});
        case 'multiple-enum'
            %TODO needs to be tested
            selectedVals = false(1,length(paramarray{i}.values));
            for j=1:length(paramValues{i})
                %TODO: handle cases when the enum is not string.
                if ismember(paramValues{i}{j},paramarray{i}.values{j})
                    selectedVals(j) = true;
                end                    
            end
            set(paramEditHandles{i},'value',find(selectedVals))
        case 'dir'
            set(paramEditHandles{i}{1},'string',paramValues{i});
        case 'file'
            set(paramEditHandles{i}{1},'string',paramValues{i});
    end
end

