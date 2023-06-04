function [ paramcell ] = fetchUIControlValues( paramarray , paramEditHandles, checkHandle)
% AUTHOR:   Abel Szkalisity
% DATE:     June 30, 2017
% NAME:     createClassImage
%
% This function extracts the parameters from uicontrol elements and returns
% them in a cellarray. If the given parameter in the control field does not
% fulfil the condition implied by the corresponding paramarray type
% (currently we only check that the 'int' typed input field must be a
% number and that directories are valid existing directories) then we show
% a warning message to the user. In that case the paramcell output is
% unreliable. NOTE: In case of default behaviour (i.e. no checkHandle
% provided), then only the currently visible items are checked. Otherwise
% it is all left to your own check function.
%
% INPUT:
%   paramarray      See the documentation of generateUIControls.
%   paramEditHandles A cellarray of handles to the uicontrol objects. We
%                   suppose that the GUI was generated with the
%                   generateUIControl function from the paramarray
%                   specification therefore paramarray and paramEditHandles
%                   has equal length. The type of the ith entry of the
%                   paramEditHandles is determined by the ith entry of the
%                   paramarray. (e.g. for enum type parameter there is a
%                   popuplist graphical object)
%   checkHandle     A handle to a function (which can be the string with
%                   the function name), that is going to be called with one
%                   parameter: the potential paramcell output. This
%                   parameter is OPTIONAL. If it is provided then the only
%                   parameter check that is run is the given function. The
%                   function must return a boolean value and a message
%                   (string) to display that describes the problem if the
%                   boolean is not 1.
%
% OUTPUT:
%   paramcell       A cellarray with the fetched parameters. The type of
%                   the ith entry of this is determined by the type of the
%                   corresponding paramarray entry. The length
%                   is equal to paramarray's length. Types given back:
%                       - int -> numeric value (can be matrix)
%                       - slider -> numberic value,
%                       - checkbox -> 0/1 (numeric)
%                       - colorPicker -> array (RGB)
%                       - enum -> string
%                       - multiple-enum -> cellarray of strings
%                       - string -> string
%                       - dir -> string
%                       - buttonGroup -> cellarray with the inner
%                       paramcells (recursively)
%
% COPYRIGHT
% Settings Template Toolbox. All Rights Reversed. 
% Copyright (C) 2019 Abel Szkalisity
% BIOMAG, Synthetic and System Biology Unit, Institute of Biochemsitry,
% Biological Research Center, Szeged, Hungary
% Ikonen group Department of Anatomy, Faculty of Medicine, University of
% Helsinki, Helsinki, Finland.

nofParams = length(paramarray);
paramcell = cell(1,nofParams);
ok = 1;
checkFuncNotProvided = nargin<3; %bool

for i=1:length(paramEditHandles)
    if strcmp(paramarray{i}.type,'str')
        paramcell{i} = get(paramEditHandles{i},'string');
        if ~isfield(paramarray{i},'strSpecs')
            strSpecs = struct('forbidden', 0, 'capital', [], 'number', [], 'glyph', [], 'space', []);
        else
            strSpecs = paramarray{i}.strSpecs;
        end                
        editHandle = paramEditHandles{i};
        if checkFuncNotProvided && strcmp(editHandle.Visible,'on')
            [ok,msg] = checkString(paramcell{i},strSpecs.forbidden,strSpecs.capital,strSpecs.number,strSpecs.glyph,strSpecs.space,paramarray{i}.name);
            if ~ok               
                break;
            end
        end
    elseif strcmp(paramarray{i}.type,'int')
        paramcell{i} = str2num(get(paramEditHandles{i},'string')); %#ok<ST2NM> It can happen that we need convert into array.        
        %check for number, not for int.
        if ~isfield(paramarray{i},'numSpecs')
            numSpecs = struct('integer',0,'scalar',0,'limits',[-Inf,Inf]);
        else
            numSpecs = paramarray{i}.numSpecs;
        end              
        editHandle = paramEditHandles{i};
        if checkFuncNotProvided && strcmp(editHandle.Visible,'on')
            [ok,msg] = checkNumber(paramcell{i},numSpecs.integer,numSpecs.scalar,numSpecs.limits,paramarray{i}.name);
            if ~ok               
                break;
            end
        end
    elseif strcmp(paramarray{i}.type,'enum') 
        paramcell{i} = paramarray{i}.values{(get(paramEditHandles{i},'value'))};
    elseif strcmp(paramarray{i}.type,'multiple-enum')
        paramcell{i} = paramarray{i}.values(get(paramEditHandles{i},'value'));
    elseif strcmp(paramarray{i}.type,'slider') || strcmp(paramarray{i}.type,'checkbox')        
        paramcell{i} = get(paramEditHandles{i},'Value');
    elseif strcmp(paramarray{i}.type,'colorPicker')
        paramcell{i} = get(paramEditHandles{i},'BackgroundColor');
    elseif strcmp(paramarray{i}.type,'dir')
        paramcell{i} = get(paramEditHandles{i}{1},'String');
        editHandle = paramEditHandles{i}{1};
        if checkFuncNotProvided && strcmp(editHandle.Visible,'on') && ~exist(paramcell{i},'dir')
            msg = ['Parameter ''' paramarray{i}.name ''' has to be an existing directory!'];
            ok = 0;
            break;
        end
    elseif strcmp(paramarray{i}.type,'file')
        paramcell{i} = get(paramEditHandles{i}{1},'String');        
        editHandle = paramEditHandles{i}{1};
        if strcmp(paramarray{i}.fileOpenType,'get') && checkFuncNotProvided && strcmp(editHandle.Visible,'on') && ~exist(paramcell{i},'file')
            msg = ['Parameter ''' paramarray{i}.name ''' does not exist on the file system!'];
            ok = 0;
            break;
        end
    elseif strcmp(paramarray{i}.type,'buttonGroup')
        nofSubGroups = length(paramEditHandles{i})-1;
        paramcell{i} = cell(nofSubGroups,1);
        for j=1:nofSubGroups
            parEditHandles = paramEditHandles{i}{j+1}; %Do like this as currently only level 2 indexing works for the myHandleCellArrays class.
            if checkFuncNotProvided               
                paramcell{i}{j} = fetchUIControlValues(paramarray{i}.groupFields,parEditHandles(1:end-1));
            else
                %If there was a check function provided then on lower levels
                %it's just a fake one
                fakeCheckFcn = @(x)deal(1,'Fake check function');            
                paramcell{i}{j} = fetchUIControlValues(paramarray{i}.groupFields,parEditHandles(1:end-1),fakeCheckFcn);
            end
        end
    end
end

if nargin>=3
    [ok,msg] = feval(checkHandle,paramcell);
end

if ~ok
    errordlg(msg,'Parameter error');
    %throw error unreliable output
    error('Settings_GUI:errorConstraintFault',msg);
end

end