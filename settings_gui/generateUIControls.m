function [ paramTextHandles, paramEditHandles] = generateUIControls( paramarray, parentHandle, heightOfItems, margin, varargin)
% AUTHOR:   Abel Szkalisity
% DATE:     Dec 09, 2016
% NAME:     generateUIControls
%
% To create the 'form' of any method on graphical user interface. The
% margins affect the appearance as described in placeUIControls.
%
%   The most important outcome of the function is that the uicontrols are
%   placed on the parent control in a dinamic way according to paramarray.
%   However for further usage the handles to the uicontrols are also
%   returned.
%
% INPUT:
%   paramarray      a cellarrray where each entry is a structure. Each
%                   structure describes one input parameter. This structure
%                   has the following fields:
%                   .name: a string with the name of the parameter (this
%                   will be printed to the user when asked for the
%                   parameter value on the GUI)
%                   .type: This determines the type of the parameter. It is
%                   a string and should be one of the followings:
%                       - 'int'** for numbers [numeric] See the optional
%                       parameters for other specifications (numSpecs)
%                       - 'enum'** Possible values listed in .values
%                       cellarray,
%                       .default describes the index in the list. Displays
%                       as a drop-down list [index of the entry] 
%                       - 'str' standing for a text input [string] See the
%                       optional parameters for other specifications (strSpecs)
%                       - 'slider'** Limits are listed in .lim [double]
%                       - 'checkbox' 0 or 1 and displays a checkbox [int: 0 / 1]
%                       - 'colorPicker' RGB triplet to store [RGB triplet
%                       each entry between 0-1]
%                       - 'multiple-enum'** Displays a list from which
%                       multiple entries can be selected [vector of indices
%                       to be selected]
%                       - 'dir' Displays a browse button or the user may
%                       copy-paste the directory into the field [string]
%                       - 'file'** Displays a browse button or the user may
%                       copy-paste the file path into the field It can be
%                       used to query file name for a new file OR also for
%                       selecting an existing file on the system (see **
%                       for the fileOpenType field)[string]
%
%                       Meta-fields:
%                       - 'buttonGroup'** Creates a set of buttons (with
%                       any number of buttons and any type described above)
%                       which may be dynamically extended or removed by the
%                       user. The basic structure to be repeated must be
%                       described in the .groupFields field. The number of
%                       group of possible appearances must be described in
%                       the groupCount field. [cellarray where each entry
%                       is a cellarray itself. The outer cellarray
%                       describes the group multiplications, the inner one
%                       the fields within that].
%                       
%                   .default: Optional. A default value to put to the input
%                   in advance. The type of the default field is indicated
%                   in brackets "[]" above.
%                   .doc: Optional. Documentation/extra information to the
%                   given parameter as string.
%                   .dynamicVis: Optional. A boolean value which indicates
%                   whether the change in answer to the given input field
%                   should initiate the call of the changeVisibility
%                   function. Default is false for all fields. Note: for
%                   proper functionality the changeVisibility function in
%                   the name-valie arguments has to be specified.
%
%                   ** There can be additional fields according to the
%                   specified type of the input.
%                   .values: This is a cellarray of strings and required
%                   for 'enum' and 'multiple-enum'. It should contain all
%                   the options to be listed.
%                   .lim: A vector with length of 2 specifying the minimum
%                   and maximum (in this order) values to be allowed for
%                   the slider.
%                   .fileOpenType: This a string value for file selection.
%                   Use either 'get' or 'put'. Use get if you want to query
%                   for an existing file and put for new files
%                   .groupFields: A cellarray of structures (just like
%                   paramarray) describing the basic structure of the
%                   gruop. NOTE: default values for these are specified
%                   separately in the .default field of the buttonGroup
%                   typed structure.
%                   .groupCount: An array with length 2 containing the
%                   minimum and maximum number of group multiplications.
%                   .filter (Optional) A filter specifying which filetypes
%                   to be able to select via the browser button. Check the
%                   documentation of uigetfile and uiputfile. However this
%                   does not  prevent the user from copying an inproper
%                   file into the text field.
%                   .numSpecs (Optional) A structure with the following
%                   fields:
%                       .integer A boolean indicating whether the input
%                       should be integer
%                       .scalar A boolean indicating whether the input has
%                       to be scalar
%                       .limits A 2 element vector of numbers indicating
%                       the acceptable limits of the number (inclusive)
%                   .strSpecs (Optional) TODO.
%                       
%   parentHandle    handle to the parent object to which we're going to add
%                   the uicontrols. It must have a position property.
%   heightOfItems   The height of each ui item and the distance between
%                   them as well in pixels. (in the drawing at the top this
%                   is the height of a row). The multiple-enum is an
%                   exemption as it resizes its height according to the
%                   number of options listed there.
%   margins         Array with two values: HM and VM (horizontal and
%                   vertical margin) in this order, and in the unit of
%                   pixels.
%   
%   NAME,VALUE pairs, optional arguments
%   resizeParent    A bool value indicating whether to resize the parent ui
%                   object so that the data fits into it (only height).
%                   Default is false
%   lineSpacing     Display parameter, see placeUIControls for details.
%   divisionRatio   Display parameter, see placeUIControls for details.
%   checkVisibility Handle to a function that defines the elements'
%                   visibility based on the current answers. It is executed
%                   upon load (hence the default value modifications take
%                   effect), and upon modifying an answer in a field with
%                   true dynamicVis property. The function gets the
%                   paramarray, the current set of answers to them, and
%                   must return a boolean array, in the same format as the
%                   current answers (the answers follow the description in
%                   the output of fetchUIControlValues)
%   checkHandle     It is possible to provide already the handle
%                   checkFunction for the fetchUIControlValues function
%                   (see description there). By doing so, the
%                   checkVisibility function is called ONLY if the
%                   parameters are faithful to the checkFunction. Note: the
%                   empty checkhandle indicates to run the default checks
%                   (such as numSpecs, strSpecs, directory exists etc.)
%   
%
% OUTPUT:
%   paramTextHandles handles to the parameter names.
%   paramEditHandles handles to the input fields.
%
% COPYRIGHT
% Settings Template Toolbox. All Rights Reversed. 
% Copyright (C) 2020 Abel Szkalisity
% BIOMAG, Synthetic and System Biology Unit, Institute of Biochemsitry,
% Biological Research Center, Szeged, Hungary
% Ikonen group Department of Anatomy, Faculty of Medicine, University of
% Helsinki, Helsinki, Finland.

% TODO: make dynamicVis working also for files, directories AND
% buttonGroups

p = inputParser;
addParameter(p,'resizeParent',false,@islogical);
addParameter(p,'lineSpacing',0.5,@isnumeric);
addParameter(p,'divisionRatio',0.5,@isnumeric);
addParameter(p,'checkVisibility',@(a,b)(defaultCheckVisibility(a,b,true)));
addParameter(p,'checkHandle',@(x)deal(1,'Fake check function'));
parse(p,varargin{:});
resizeParent = p.Results.resizeParent;
lineSpacing = p.Results.lineSpacing;
divisionRatio = p.Results.divisionRatio;

nofParams = length(paramarray);

paramTextHandles = myHandleCellArray(1,nofParams);
paramEditHandles = myHandleCellArray(1,nofParams);
dataCollection.pth = paramTextHandles;
dataCollection.peh = paramEditHandles;
dataCollection.hoi = heightOfItems;
dataCollection.m = margin;
dataCollection.ls = lineSpacing;
dataCollection.dr = divisionRatio;
dataCollection.fullParamArray = paramarray;
dataCollection.checkHandle = p.Results.checkHandle;
dataCollection.checkVisibility = p.Results.checkVisibility;
dataCollection.resizeParent = resizeParent;
dataCollection.parentHandle = parentHandle;

sizeOfPanel = get(parentHandle,'Position');

%Create handles
for i=1:nofParams
    [paramTextHandles{i},paramEditHandles{i}] = generateItem(paramarray{i},parentHandle,dataCollection);
end

%Place them dynamically
[sumHeightItems] = placeUIControls(...
    paramarray, paramTextHandles, paramEditHandles, sizeOfPanel, heightOfItems, margin,...
    'lineSpacing',lineSpacing,...
    'divisionRatio',divisionRatio,...
    'checkHandle',p.Results.checkHandle,...
    'checkVisibility',p.Results.checkVisibility);

if resizeParent
    parentHandle.Position(4) = sumHeightItems;
    sizeOfPanel(4) = sumHeightItems;
    placeUIControls(...
        paramarray, paramTextHandles, paramEditHandles, sizeOfPanel, heightOfItems, margin,...
        'lineSpacing',lineSpacing,...
        'divisionRatio',divisionRatio,...
        'checkHandle',p.Results.checkHandle,...
        'checkVisibility',p.Results.checkVisibility);
end
   
end

%% Internal functions
function [textHandle,editHandle] = generateItem(paramStruct,parentHandle,dataCollection)
    %Creating text handles
    if strcmp(paramStruct.type,'buttonGroup')
        actTextControl = uicontrol('Parent',parentHandle,'Units','pixels','Style','Text','String',['--- Specify ' paramStruct.name ' parameters ---']);
    else
        actTextControl = uicontrol('Parent',parentHandle,'Units','pixels','Style','Text','String',paramStruct.name);
    end
    if isfield(paramStruct,'doc')
        infoIcon = im2double(imread('infoIcon.png'));
        bgPix = all(infoIcon == 1,3);        
        nofAllBgPix = sum(sum(bgPix));
        infoIcon(repmat(bgPix,1,1,3)) = repmat(reshape(parentHandle.BackgroundColor,1,1,3),nofAllBgPix,1);
        textHandle{1} = uicontrol('Parent',parentHandle,'Units','pixels','Style','pushbutton','CData',infoIcon,'tooltip','Click for documentation',...
            'UserData',struct('docText',paramStruct.doc,'paramName',paramStruct.name,'parentHandle',parentHandle),...
            'Callback',@(hObject,eventdata)docCallback(hObject,eventdata,guidata(hObject)));
        textHandle{2} = actTextControl;
    else
        textHandle = actTextControl;
    end      
    
    %Creating the different edit handles
    switch paramStruct.type
        case 'int'
            editHandle = uicontrol('Parent',parentHandle,'Units','pixels','Style','edit','BackgroundColor','white');
            if isfield(paramStruct,'default')
                set(editHandle,'string',num2str(paramStruct.default));
            end            
            if isfield(paramStruct,'dynamicVis') && paramStruct.dynamicVis, set(editHandle,'Callback',@(hO,ed)dynamicCallback(hO,ed,dataCollection,[])); end
        case 'str'
            editHandle = uicontrol('Parent',parentHandle,'Units','pixels','Style','edit','BackgroundColor','white');
            if isfield(paramStruct,'default')
                set(editHandle,'string',paramStruct.default);
            end            
            if isfield(paramStruct,'dynamicVis') && paramStruct.dynamicVis, set(editHandle,'Callback',@(hO,ed)dynamicCallback(hO,ed,dataCollection,[])); end
        case 'enum'        
            editHandle = uicontrol('Parent',parentHandle,'Units','pixels','Style','popupmenu','BackgroundColor','white');
            set(editHandle,'String',paramStruct.values);
            if isfield(paramStruct,'default')
                set(editHandle,'value',paramStruct.default);
            end            
            if isfield(paramStruct,'dynamicVis') && paramStruct.dynamicVis, set(editHandle,'Callback',@(hO,ed)dynamicCallback(hO,ed,dataCollection,[])); end
        case 'slider'
            editHandle = uicontrol('Parent',parentHandle,'Units','pixels','Style','slider','BackgroundColor','white');
            set(editHandle,'Min',paramStruct.lim(1));
            set(editHandle,'Max',paramStruct.lim(2));
            if isfield(paramStruct,'default')
                set(editHandle,'value',paramStruct.default);
            end            
            if isfield(paramStruct,'dynamicVis') && paramStruct.dynamicVis, set(editHandle,'Callback',@(hO,ed)dynamicCallback(hO,ed,dataCollection,[])); end
        case 'checkbox'
            editHandle = uicontrol('Parent',parentHandle,'Units','pixels','Style','checkbox');
            if isfield(paramStruct,'default')
                set(editHandle,'value',paramStruct.default);
            end            
            if isfield(paramStruct,'dynamicVis') && paramStruct.dynamicVis, set(editHandle,'Callback',@(hO,ed)dynamicCallback(hO,ed,dataCollection,[])); end
        case 'colorPicker'
            editHandle = uicontrol('Parent',parentHandle,'Units','pixels','Style','pushbutton','BackgroundColor','white','Callback',@(hObject,eventdata)colorPickerCallback(hObject,eventdata,guidata(hObject)));
            if isfield(paramStruct,'default')
                set(editHandle,'BackgroundColor',paramStruct.default);
            end            
            if isfield(paramStruct,'dynamicVis') && paramStruct.dynamicVis
                origCB = editHandle.Callback;
                set(editHandle,'Callback',@(hO,ed)dynamicCallback(hO,ed,dataCollection,origCB)); 
            end
        case 'multiple-enum'
            editHandle = uicontrol('Parent',parentHandle,'Units','pixels','Style','listbox','BackgroundColor','white','String',paramStruct.values,'min',0,'max',length(paramStruct.values));            
            if isfield(paramStruct,'default')
                set(editHandle,'value',paramStruct.default);
            end
            if isfield(paramStruct,'dynamicVis') && paramStruct.dynamicVis, set(editHandle,'Callback',@(hO,ed)dynamicCallback(hO,ed,dataCollection,[])); end
        case 'dir'
            editHandle{1} = uicontrol('Parent',parentHandle,'Units','pixels','Style','edit','BackgroundColor','white');
            editHandle{2} = uicontrol('Parent',parentHandle,'Units','pixels','Style','pushbutton','String','browse',...
                'UserData',struct('paramStruct',paramStruct,'editHandle',editHandle{1}),...
                'Callback',@(hObject,eventdata)browseCallback(hObject,eventdata,guidata(hObject),'dir'));
            if isfield(paramStruct,'default')
                set(editHandle{1},'string',paramStruct.default);
            end            
        case 'file'
            editHandle{1} = uicontrol('Parent',parentHandle,'Units','pixels','Style','edit','BackgroundColor','white');
            editHandle{2} = uicontrol('Parent',parentHandle,'Units','pixels','Style','pushbutton','String','browse',...
                'UserData',struct('paramStruct',paramStruct,'editHandle',editHandle{1}),...
                'Callback',@(hObject,eventdata)browseCallback(hObject,eventdata,guidata(hObject),'file'));
            if isfield(paramStruct,'default')
                set(editHandle{1},'string',paramStruct.default);
            end
        case 'buttonGroup'            
            if isfield(paramStruct,'default')
                nofDefaultGroups = length(paramStruct.default); 
            else
                nofDefaultGroups = paramStruct.groupCount(1); %Set it to minimum
            end
            if paramStruct.groupCount(1)<0, error('Settings_GUI:definitionFault','The minimum number of groups must be at least 0.'); end
            if nofDefaultGroups>paramStruct.groupCount(2) %Maximize the number of groups.
                nofDefaultGroups = paramStruct.groupCount(2);
            end
            
            %Note the n+1 st elemetns are used internally for other
            %reasons.
            infoLine = textHandle;
            textHandle = myHandleCellArray(nofDefaultGroups+1,1);
            editHandle = myHandleCellArray(nofDefaultGroups+1,1);
            for j=(1:nofDefaultGroups)+1
                [textHandle{j},editHandle{j}] = createSingleButtonGroup(paramStruct,parentHandle,dataCollection,j-1,textHandle,editHandle);
            end
            %Put it to the beginning
            editHandle{1} = uicontrol('Parent',parentHandle,'Units','pixels','Style','pushbutton',...
                'String',['Add another ' paramStruct.name],...
                'UserData',struct(...
                            'paramStruct',paramStruct,...
                            'ownEditHandles',editHandle,...
                            'ownTextHandles',textHandle,...                            
                            'parentHandle',parentHandle,...
                            'dataCollection',dataCollection),...
                'Callback',@(hObject,eventdata)changeButtonGroupNumberCallback(hObject,eventdata,guidata(hObject),'add'));
            textHandle{1} = infoLine;
    end            
end

function colorPickerCallback(hObject,~,~)
    pickedColor = uisetcolor();
    if length(pickedColor)==3
        set(hObject,'BackgroundColor',pickedColor);
    end      
end

function docCallback(hObject,~,~)
    uData = get(hObject,'UserData');
    M = msgbox(uData.docText,['Documentation for: ' uData.paramName ],'help');
    parentUData = uData.parentHandle.UserData;
    if isempty(parentUData), parentUData.docFigs = {M};
    else, parentUData.docFigs{end+1} = M;
    end
    uData.parentHandle.UserData = parentUData;
    origDelFunc = uData.parentHandle.DeleteFcn;
    uData.parentHandle.DeleteFcn = @(hO,ed)(settingsPanelCloseRequestFunc(hO,ed,origDelFunc));
end

function settingsPanelCloseRequestFunc(hObject,eventdata,origPanelFunc)
    uData = get(hObject,'UserData');
    if isfield(uData,'docFigs')
        for i=1:length(uData.docFigs)
            if ishandle(uData.docFigs{i}), delete(uData.docFigs{i}); end
        end
    end        
    if ~isempty(origPanelFunc)
        origPanelFunc(hObject,eventdata);
    end
    % If it is still a handle
    if ishandle(hObject), delete(hObject); end
end

function browseCallback(hObject,~,~,fileType)
%Filytype should be dir or file
    uData = get(hObject,'UserData');
    defPath = get(uData.editHandle,'String');
    if isempty(defPath)
        defPath = pwd;           
    end    
    switch fileType
        case 'dir'
            fileToOpen = uigetdir(defPath,'Please select a folder');
        case 'file'
            if isfield(uData.paramStruct,'filter')
                filter = uData.paramStruct.filter;
            else
                filter = '*';
            end
            switch uData.paramStruct.fileOpenType
                case 'get'                    
                    [fileName,filePath] = uigetfile(filter,'Please select a file',defPath);                
                case 'put'
                    [fileName,filePath] = uiputfile(filter,'Please select a file',defPath);
            end
            if ~isnumeric(fileName) && ~isnumeric(filePath)
                fileToOpen = fullfile(filePath,fileName);
            else
                fileToOpen = 0;
            end                
    end        
    if ~isnumeric(fileToOpen)
        set(uData.editHandle,'String',fileToOpen);
    end
end

function dynamicCallback(hObject,eventdata,dataCollection,origCallback) 
    if ~isempty(origCallback)
        origCallback(hObject,eventdata);
    end
    parentHandle = dataCollection.parentHandle;
    sizeOfPanel = get(parentHandle,'Position');
    [sumHeightItems] = placeUIControls(...
        dataCollection.fullParamArray, dataCollection.pth, dataCollection.peh, sizeOfPanel , dataCollection.hoi, dataCollection.m,...
        'lineSpacing',dataCollection.ls ,...
        'divisionRatio',dataCollection.dr,...
        'checkHandle',dataCollection.checkHandle,...
        'checkVisibility',dataCollection.checkVisibility);
    if dataCollection.resizeParent
        parentHandle.Position(4) = sumHeightItems;
        sizeOfPanel(4) = sumHeightItems;
        placeUIControls(...
            dataCollection.fullParamArray, dataCollection.pth, dataCollection.peh, sizeOfPanel , dataCollection.hoi, dataCollection.m,...
            'lineSpacing',dataCollection.ls ,...
            'divisionRatio',dataCollection.dr,...
            'checkHandle',dataCollection.checkHandle,...
            'checkVisibility',dataCollection.checkVisibility);
    end
end


function changeButtonGroupNumberCallback(hObject,~,~,changeType)
    uData = get(hObject,'UserData');        
    oth = uData.ownTextHandles;
    oeh = uData.ownEditHandles;    
    DC = uData.dataCollection;    
    parentPanel = uData.parentHandle;        
    switch changeType
        case 'add'
            paramStruct = uData.paramStruct;
            %Always add to the end. Always check for the greatest ID.
            maxID = 0;
            for i=2:length(oth), if oth{i}{length(oth{i})}>maxID, maxID = oth{i}{length(oth{i})}; end; end
            maxID = maxID + 1;
            [oth{length(oth)+1},oeh{length(oeh)+1}] = createSingleButtonGroup(paramStruct,parentPanel,DC,maxID,oth,oeh); %#ok<NASGU> Fine, it is supposed to be a handle
            sizeOfPanel = get(parentPanel,'Position');
                        
            sumHeightItems = placeUIControls( DC.fullParamArray, DC.pth, DC.peh, sizeOfPanel, DC.hoi, DC.m,...
                'lineSpacing',DC.ls,...
                'divisionRatio',DC.dr,...
                'checkHandle',DC.checkHandle,...
                'checkVisibility',DC.checkVisibility);
            
            parentPanel.Position(4) = sumHeightItems;
            sizeOfPanel(4) = sumHeightItems;
            placeUIControls( DC.fullParamArray, DC.pth, DC.peh, sizeOfPanel, DC.hoi, DC.m,...
                'lineSpacing',DC.ls,...
                'divisionRatio',DC.dr,...
                'checkHandle',DC.checkHandle,...
                'checkVisibility',DC.checkVisibility);
            
        case 'remove'
            myID = uData.myID;
            for i=2:length(oth)
                if myID == oth{i}{length(oth{i})}
                    idx = i;
                    break;
                end
            end
            deleteHandlesRecursively(oth{idx});
            deleteHandlesRecursively(oeh{idx});
            oth(idx) = []; %#ok<NASGU> It is a handle
            oeh(idx) = []; %#ok<NASGU> It is a handle
            placeUIControls( DC.fullParamArray, DC.pth, DC.peh, get(parentPanel,'Position'), DC.hoi, DC.m,...
                'lineSpacing',DC.ls,...
                'divisionRatio',DC.dr,...
                'checkHandle',DC.checkHandle,...
                'checkVisibility',DC.checkVisibility);
            % In case of deletion don't resize the parent
    end
end

function [textHandle,editHandle] = createSingleButtonGroup(paramStruct,parentHandle,dataCollection,j,textHandles,editHandles)

    nofBasicFields = length(paramStruct.groupFields);
    textHandle = cell(nofBasicFields,1);
    editHandle = cell(nofBasicFields,1);
    for k=1:nofBasicFields
        if isfield(paramStruct,'default')
            currField = paramStruct.groupFields{k};
            jLocal = mod(j-1,length(paramStruct.default))+1;
            currField.default = paramStruct.default{jLocal}{k};
            [textHandle{k},editHandle{k}] =  generateItem(currField,parentHandle,dataCollection);
        else
            [textHandle{k},editHandle{k}] =  generateItem(paramStruct.groupFields{k},parentHandle,dataCollection);
        end
    end
    textHandle{k+1} = j;
    editHandle{k+1} = uicontrol('Parent',parentHandle,'Units','pixels','Style','pushbutton',...
        'String',['Remove this ' paramStruct.name],...
        'UserData',struct(...            
            'ownEditHandles',editHandles,...
            'ownTextHandles',textHandles,...
            'myID',j,...
            'parentHandle',parentHandle,...
            'dataCollection',dataCollection),...
        'Callback',@(hObject,eventdata)changeButtonGroupNumberCallback(hObject,eventdata,guidata(hObject),'remove'));
end

function deleteHandlesRecursively(cellA)
    if iscell(cellA)
        for i=1:length(cellA)
            if length(cellA{i})==1 && ishandle(cellA{i}) && ~isa(cellA{i},'myHandleCellArray')
                delete(cellA{i});
            elseif iscell(cellA{i}) || isa(cellA{i},'myHandleCellArray')
                deleteHandlesRecursively(cellA{i});
            end
        end
    end
end

