function [sumHeightItems] = placeUIControls( paramarray, paramTextHandles, paramEditHandles, sizeOfPanel, heightOfItems, margin, varargin)
% AUTHOR:   Abel Szkalisity
% DATE:     Oct 20, 2019
% NAME:     placeUIControls
%
% This function calculates how to place the already created graphical
% objects onto the parent panel.
%
%   | ^VM
%   | <- HM -> <- info -> PARAMETER1 <- HM ->/2 || <- HM ->/2 INPUT_FIELD1 <- HM -> |
%   | lineSpacing
%   | <- HM -> <- info -> PARAMETER2 <- HM ->/2 || <- HM ->/2 INPUT_FIELD2 <- HM -> |
%   | ?VM
%   
%   | < ----------- text size --------------- > || < ------- edit size ---------- > |
%   
%
%   VM and HM stands for Vertical and Horizontal Margin and come as input.
%   
%   In case of an info icon, it is placed into the left panel right before
%
% INPUTS:
%   paramarray      See in generateUIControls. This is required as some
%                   inputs height is dependent on the inner data.
%   paramTextHandles,paramEditHandles Cellarrays of handle structures as
%                   created by generateUIControls.
%   sizeOfPanel     The size property of the parent handle.
%   heightOfItems   The height of each ui item and the double of distance
%                   between them as well in pixels. (in the drawing at the
%                   top this is the height of a row). The multiple-enum is
%                   an exemption as it resizes its height according to the
%                   number of options listed there.
%   margins         Array with two values: HM and VM (horizontal and
%                   vertical margin) in this order, and in the unit of
%                   pixels.
%   NAME,VALUE pairs, optional arguments
%   lineSpacing     A double value for the size of line spacing
%                   (multiplication factor for the height of items)
%   groupSpacing    A double value for the size of the space before and
%                   after grouped itmes (multiplication factor for the
%                   height of items)
%   divisionRatio   The horizontal size ratio of the left (text, parameter
%                   description) and right (edit boxes) panels in the GUI.
%                   The value refers to the size of the left panel compared
%                   to the value 1. I.e. the default 0.5 value halves the
%                   GUI, but a 0.3 value would make the texts cover 1/3 of
%                   the total width and leave 2/3 space for the edit boxes.
%                   (textSize/(textSize+editSize) == divisionRatio)
%   checkVisibility See the doc of generateUIControls
%   checkHandle     See the doc of the fetchUIControlValues
%
% OUTPUT:
%   sumHeightItems  The total sum of the height of the placed items in
%                   pixels
%
% COPYRIGHT
% Settings Template Toolbox. All Rights Reversed. 
% Copyright (C) 2019 Abel Szkalisity
% BIOMAG, Synthetic and System Biology Unit, Institute of Biochemsitry,
% Biological Research Center, Szeged, Hungary
% Ikonen group Department of Anatomy, Faculty of Medicine, University of
% Helsinki, Helsinki, Finland.

%TODO: do checks for negative widths

p = inputParser;
addParameter(p,'lineSpacing',0.5,@isnumeric);
addParameter(p,'groupSpacing',0.5,@isnumeric);
addParameter(p,'divisionRatio',0.5,@(x)(isnumeric(x) && x>0 && x<1))
addParameter(p,'checkVisibility',@(a,b)(defaultCheckVisibility(a,b,true)));
addParameter(p,'checkHandle',@(x)deal(1,'Fake check function'));
parse(p,varargin{:});
lineSpacing = p.Results.lineSpacing;
groupSpacing = p.Results.groupSpacing;
textRatio = p.Results.divisionRatio;

nofParams = length(paramarray);

sumHeightItems = 0;
maxMultipleEnumHeight = 200;

%Try to fetch here the current values
changeVis = false;
try 
    if isempty(p.Results.checkHandle)
        [ paramcell ] = fetchUIControlValues( paramarray , paramEditHandles);
    else
        [ paramcell ] = fetchUIControlValues( paramarray , paramEditHandles, p.Results.checkHandle);
    end    
    changeVis = true;
catch e
    if ~strcmp(e.identifier,'Settings_GUI:errorConstraintFault')           
        rethrow(e);
    end
end
if changeVis
    displayItems = p.Results.checkVisibility(paramarray,paramcell);
else
    displayItems = checkCurrentVisibility(paramarray,paramEditHandles,paramTextHandles);
end

%Set here if some field needs to be greater than standard (i.e. multiple list input)
for i=1:nofParams
    [sumHeightItems] = placeSingleItem(paramarray{i},paramTextHandles{i},paramEditHandles{i},sumHeightItems,displayItems{i});
end

sumHeightItems = sumHeightItems + 2*margin(2);

%Local function for recursion
function [sumHeightItems] = placeSingleItem(paramStruct,paramTextHandle,paramEditHandle,sumHeightItems,displayItem)
    currentHeight = heightOfItems;
    editColWidth = sizeOfPanel(3)*(1-textRatio)-margin(1)*1.5;
    if strcmp(paramStruct.type,'buttonGroup')
               
        if displayItem{1}
            %Set the header for the button group. This one covers the whole
            %width.
            sumHeightItems = sumHeightItems + currentHeight + heightOfItems*groupSpacing;
            placeTextHandle(paramTextHandle{1},margin,sizeOfPanel,sumHeightItems, textRatio, currentHeight, 'full');        
            %Add regular margin also below this.
            sumHeightItems = sumHeightItems + heightOfItems*lineSpacing;

            for j=2:length(paramTextHandle)                
                for k=1:length(paramStruct.groupFields)
                    [sumHeightItems] = placeSingleItem(paramStruct.groupFields{k},paramTextHandle{j}{k},paramEditHandle{j}{k},sumHeightItems,displayItem{j}{k});
                end

                % After each item we need to place the remove button
                sumHeightItems = sumHeightItems + currentHeight;
                set(paramEditHandle{j}{k+1},'Position',[sizeOfPanel(3)*textRatio+margin(1)/2, sizeOfPanel(4)-margin(2)-sumHeightItems, editColWidth, currentHeight],'Visible','on');
                %Add regular margin also below this.
                sumHeightItems = sumHeightItems + heightOfItems*lineSpacing;

                %Check here for the minimality constraint to disable the
                %button.
                if length(paramTextHandle)-1==paramStruct.groupCount(1)
                    set(paramEditHandle{j}{k+1},'Enable','off');
                else
                    set(paramEditHandle{j}{k+1},'Enable','on');
                end            
            end
            if isempty(j), sumHeightItems = sumHeightItems + heightOfItems*(1+lineSpacing); end

            % After each item broup we need to place the add button. NOTE: this
            % one is in ONE LINE WITH THE LAST REMOVE BUTTON
            % Hence we need to remove the regular line spacing:
            sumHeightItems = sumHeightItems - heightOfItems*lineSpacing;
            set(paramEditHandle{1},'Position',[margin(1), sizeOfPanel(4)-margin(2)-sumHeightItems, sizeOfPanel(3)*textRatio-margin(1)*1.5, currentHeight],'Visible','on');
            %And then add it back.
            sumHeightItems = sumHeightItems + heightOfItems*lineSpacing;

            %Check also for this button's availability:
            if length(paramTextHandle)-1==paramStruct.groupCount(2)
                set(paramEditHandle{1},'Enable','off');
            else
                set(paramEditHandle{1},'Enable','on');
            end

            %Aand the group margin.
            sumHeightItems = sumHeightItems + currentHeight + heightOfItems*groupSpacing;        
        else
            hideTextHandle(paramTextHandle{1});
            % Just to make sure everything is false from down to this.    
            displayItem = assignToCell(displayItem,false);
            % Hide everything inside recursively
            for j=2:length(paramTextHandle)
                for k=1:length(paramStruct.groupFields)                    
                    [sumHeightItems] = placeSingleItem(paramStruct.groupFields{k},paramTextHandle{j}{k},paramEditHandle{j}{k},sumHeightItems,displayItem{j}{k});
                end                
                % Hide the remove button
                set(paramEditHandle{j}{k+1},'Visible','off');                               
            end            
            % Hide the addition button
            set(paramEditHandle{1},'Visible','off');
        end
        
    else        
        if displayItem
            %Modify item height for multiple enums.
            if strcmp(paramStruct.type,'multiple-enum')
                currentHeight = min(maxMultipleEnumHeight,heightOfItems*length(paramStruct.values));                            
            end

            sumHeightItems = sumHeightItems + currentHeight;        

            if strcmp(paramStruct.type,'dir') || strcmp(paramStruct.type,'file')        
                editPositionLeft =  [sizeOfPanel(3)*textRatio+margin(1)/2, sizeOfPanel(4)-margin(2)-sumHeightItems, editColWidth*2/3 currentHeight];
                editPositionRight = [sizeOfPanel(3)*textRatio+margin(1)/2+editColWidth*2/3, sizeOfPanel(4)-margin(2)-sumHeightItems, editColWidth*1/3 currentHeight];
                set(paramEditHandle{1},'Position',editPositionLeft,'Visible','on');
                set(paramEditHandle{2},'Position',editPositionRight,'Visible','on');
            else
                editPosition = [sizeOfPanel(3)*textRatio+margin(1)/2, sizeOfPanel(4)-margin(2)-sumHeightItems, editColWidth currentHeight];    
                set(paramEditHandle,'Position',editPosition,'Visible','on');
            end

            placeTextHandle(paramTextHandle,margin,sizeOfPanel,sumHeightItems, textRatio, currentHeight, 'default');

            % Add margin
            sumHeightItems = sumHeightItems + heightOfItems*lineSpacing;
        else
            hideTextHandle(paramTextHandle);
            if strcmp(paramStruct.type,'dir') || strcmp(paramStruct.type,'file')                        
                set(paramEditHandle{1},'Visible','off');
                set(paramEditHandle{2},'Visible','off');
            else                
                set(paramEditHandle,'Visible','off');
            end                
        end
    end
end


end

function placeTextHandle(myTextHandles,margin,sizeOfPanel,sumHeightItems, textRatio, currentHeight, width)
%Width refers to the fact if the text label is 'full' or anything else (default)
    infoMargin = 2;
    if strcmp(width,'full')
        w = sizeOfPanel(3)-margin(1)*2;
    else %text ratio
        w = sizeOfPanel(3)*textRatio-margin(1)*1.5;
    end
    if iscell(myTextHandles) && length(myTextHandles) == 2
        set(myTextHandles{1},'Position',[margin(1), sizeOfPanel(4)-margin(2)-sumHeightItems currentHeight currentHeight],'Visible','on');
        resizedImg = imresize(myTextHandles{1}.CData,[currentHeight currentHeight]);
        resizedImg(resizedImg>1) = 1;
        set(myTextHandles{1},'CData',resizedImg);
        set(myTextHandles{2},'Position',[margin(1)+currentHeight+infoMargin, sizeOfPanel(4)-margin(2)-sumHeightItems w-currentHeight-infoMargin  currentHeight],'Visible','on');
    else
        set(myTextHandles,'Position',[margin(1)+currentHeight+infoMargin, sizeOfPanel(4)-margin(2)-sumHeightItems w-currentHeight-infoMargin currentHeight],'Visible','on');
    end    
end

function hideTextHandle(myTextHandles)
%Width refers to the fact if the text label is 'full' or anything else (default)        
    if iscell(myTextHandles) && length(myTextHandles) == 2
        set(myTextHandles{1},'Visible','off');
        set(myTextHandles{2},'Visible','off');
    else
        set(myTextHandles,'Visible','off');
    end    
end


function result = anyForCell(cellWithLogicals)
%This function is a generalized any function for cellarrays. It returns
%true of there is at least one true value somewhere down in the array.
result = false;
for i=1:length(cellWithLogicals)
    if iscell( cellWithLogicals{i} )
        result = result || anyForCell(cellWithLogicals{i});
    elseif islogical(cellWithLogicals{i})
        result = result || cellWithLogicals{i};
    end
end
    
end

function cellArray = assignToCell(cellArray,valueToAssign)
%This function recursively assigns a value to each entry in the cellarray
%which is not a cellarray itself.
    if iscell( cellArray )
        for i=1:length(cellArray)
            cellArray{i} = assignToCell(cellArray{i},valueToAssign);
        end
    else
        cellArray = valueToAssign;
    end
    
end


function visibilityBool = checkCurrentVisibility(paramarray,paramEditHandles,paramTextHandles)
%Checks the current visibility ststus of the settings. Output should follow
%similar conventions as the checkVisibility functions (described in generateUIControls)
    visibilityBool = cell(size(paramarray));
    for i=1:length(paramarray)
        if strcmp(paramarray{i}.type,'buttonGroup')
            innerVisibilities = cell(length(paramEditHandles)-1,1);
            for j=2:length(paramEditHandles{i}) 
               innerVisibilities{j} = checkCurrentVisibility(paramarray{i}.groupFields,paramEditHandles{i}{j},paramTextHandles{i}{j});
            end
            visibilityBool{i} = innerVisibilities;
        else %It is not recursive here any more
            %Text handle is either a handle to the text itself OR a 2
            %length cellarray to the info icon and the text itself
            
            %Edit handle is either an edit handle, OR a length of 2 (for
            %dir and file) cell array to the handles.
            visibilityBool{i} = checkCellHandleVisibility(paramEditHandles{i}) && checkCellHandleVisibility(paramTextHandles{i});
        end
    end
end

function visible = checkCellHandleVisibility(cellOrHandle)
%Visible only if everything is visible.
    visible = true;
    if iscell(cellOrHandle)
        for i=1:length(cellOrHandle)
            visible = visible & checkCellHandleVisibility(cellOrHandle{i});
        end
    elseif ishandle(cellOrHandle)
        visible = visible & strcmp(cellOrHandle.Visible,'on');
    end
end


