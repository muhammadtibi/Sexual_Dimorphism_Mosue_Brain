function [bool,msg] = checkString(toBeChecked, forbidden, capital, number, glyph, space, text)
% AUTHOR:   Attila Beleon
% DATE:     July 1, 2020
% NAME:     checkString
%
% Checks for specific string properties of an input. If the input is not a
% valid string then it must be empty.
% The different constraing inputs can be left empty in which case that
% condition is not examined.
%
% INPUT:
%   toBeChecked     The string to be checked.
%   forbidden       Cell array contains forbidden strings. If it is set
%                   then the toBeCHecked must be different.
%   capital         Bool - Requirement of capital letter
%   number          Bool - Requirement of number
%   glyph           Bool - Requirement of glyph
%   space           Bool - Requirement of space
%
% OUTPUT:
%   bool            Boolean indicating if the process was successful or not
%   msg             If bool is 0 so that at least one constraint is not met
%                   then in the msg there is description about the problem
%
% COPYRIGHT
% Settings Template Toolbox. All Rights Reversed. 
% Copyright (C) 2019 Abel Szkalisity
% BIOMAG, Synthetic and System Biology Unit, Institute of Biochemsitry,
% Biological Research Center, Szeged, Hungary
% Ikonen group Department of Anatomy, Faculty of Medicine, University of
% Helsinki, Helsinki, Finland.


msg = '';
bool = 1;
if isempty(toBeChecked)
    bool = 0;
    msg = [ text ' can not be empty!'];
    return;
end

if ~isempty(forbidden)
    idx = any(strcmp(forbidden, toBeChecked));
    if idx
        bool = 0;
        msg = [ text ' already exist!'];
        return;
    end
end

if ~isempty(capital)
    if capital == 1 && ~any(isstrprop(toBeChecked,'upper'))
        bool = 0;
        msg = [text ' must contain capital letter!'];
        return;
    elseif capital == 0 && any(isstrprop(toBeChecked,'upper'))
        bool = 0;
        msg = [text ' must not contain capital letter!'];
        return;
    end
end

if ~isempty(number)
    if number == 1 && isempty(regexp(toBeChecked,'\d','once'))
        bool = 0;
        msg = [text ' must contain number!'];
        return;
    elseif number == 0 && ~isempty(regexp(toBeChecked,'\d','once'))
        bool = 0;
        msg = [text ' must not contain number!'];
        return;
    end
end

if ~isempty(glyph)
    if glyph == 1 && isempty(regexp(toBeChecked,'\W','once'))
        bool = 0;
        msg = [text ' must contain glyph!'];
        return;
    elseif glyph == 0 && ~isempty(regexp(toBeChecked,'\W','once'))
        bool = 0;
        msg = [text ' must not contain glyph!'];
        return;
    end
end

if ~isempty(space)
    if space == 1 && isempty(regexp(toBeChecked,'\s','once'))
        bool = 0;
        msg = [text ' must contain space!'];
        return;
    elseif space == 0 && ~isempty(regexp(toBeChecked,'\s', 'once'))
        bool = 0;
        msg = [text ' must not contain space!'];
        return;
    end
end

end