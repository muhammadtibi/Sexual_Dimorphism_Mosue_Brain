function [bool,msg] = checkNumber( toBeChecked, integer, scalar, limits, text)
% AUTHOR:   Abel Szkalisity
% DATE:     March 13, 2017
% NAME:     checkNumber
%
% Checks for specific number properties of an input. If the input is not a
% valid number then it must be empty.
% The different constraing inputs can be left empty in which case that
% condition is not examined
%
% INPUT:
%   toBeChecked     The number to be checked, empty or scalar or matrix.
%   integer         Boolean checks for integer numbers. If it is set then
%                   all elements of toBeChecked must be integer.
%   scalar          The cardinality of the matrix must be 1.
%   limits          2 element vector, specifying inclusive limits for
%                   values of toBeChecked. Use -inf and inf if the limit is
%                   needed only from one way
%   text            The text to be included in the output msg before the
%                   output of the function.
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
    msg = [ text ' has to be a number!'];
    return;
end

if ~isempty(integer) && integer && any(reshape(rem(toBeChecked,1)~=0,1,[]))
    bool = 0;
    msg = [ text ' has to be integer!'];
    return;
end

if ~isempty(scalar) && scalar && numel(toBeChecked)~=1
    bool = 0;
    msg = [text ' has to be scalar!'];
    return;
end

if ~isempty(limits)
    if min(toBeChecked(:))<limits(1)
        bool = 0;
        msg = [ text ' must be at least ' num2str(limits(1)) ' !'];
    end
    if max(toBeChecked(:))>limits(2)
        bool = 0;
        msg = [ text ' must be at most ' num2str(limits(2)) ' !'];
    end
end


end