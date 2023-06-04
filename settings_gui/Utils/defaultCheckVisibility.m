function visibilityBool = defaultCheckVisibility(paramarray,paramcell,bool)
%Must return a cellarray describing the visibility of each element in
%paramcell with a boolean value. The overall structure must follow the
%structure of paramcell, except for the case when paramarray is a button
%group. In this case the corresponding cellarray entry is one longer the
%first one referring to the general visibility of the button group.
%NOTE: the whole button group's visibility always depends on the main
%visibility option, regardless of what is visibilities are specified
%internally within the group. However the structure still must follow tha
%paramcell's structure!
    visibilityBool = cell(size(paramcell));
    for i=1:length(paramarray)
        if strcmp(paramarray{i}.type,'buttonGroup')
            nofSubGroups = length(paramcell{i});
            visibilityBool{i} = cell(1,nofSubGroups+1);
            visibilityBool{i}{1} = bool;
            for j=1:nofSubGroups
                visibilityBool{i}{j+1} = defaultCheckVisibility(paramarray{i}.groupFields,paramcell{i}{j},bool);
            end
        else
            visibilityBool{i} = bool;
        end
    end
end