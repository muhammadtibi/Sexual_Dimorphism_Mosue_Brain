classdef myHandleCellArray < handle
    %To facilitate the dynamic change of cell array values.
    
    properties
        traditionalCellArray;
    end
    
    methods
        function obj = myHandleCellArray(varargin)
            obj.traditionalCellArray = cell(varargin{:});
        end
        
        function varargout = subsref(obj,s)
            switch s(1).type
                case '.'
                    if length(s)==1
                        varargout{1} = obj.(s.subs);
                    elseif length(s)==2 
                         [varargout{1:nargout}] = subsref(obj.(s(1).subs),s(2));
                    else
                         [varargout{1:nargout}] = subsref(obj,s);
                    end                        
                case '{}' 
                    if length(s)==1
                        [varargout{1:nargout}] = subsref(obj.traditionalCellArray,s);
                    elseif length(s)==2
                        %Get the first level: it is curly brace so it must
                        %have length 1.
                        firstLevel = subsref(obj.traditionalCellArray,s(1));
                        [varargout{1:nargout}] = subsref(firstLevel,s(2));
                    else
                        builtin(subsref,obj,s);
                    end
                case '()'
                    varargout = {subsref(obj.traditionalCellArray,s)};
                otherwise
                    error('Not a valid indexing expression')    
            end            
        end
        
        function obj = subsasgn(obj,s,varargin)
            switch s(1).type
                case '.'
                    if length(s)==1
                        obj.(s.subs) = varargin{:};
                    elseif length(s)==2 && strcmp(s(1).subs,'traditionalCellArray')
                        obj.traditionalCellArray = builtin('subsasgn',obj.traditionalCellArray,s(2),varargin{:});
                    else
                        
                    end                        
                case '{}'
                    if length(s) == 1
                        obj.traditionalCellArray = builtin('subsasgn',obj.traditionalCellArray,s,varargin{:});
                    elseif length(s) > 1
                        innerLevel = subsref(obj.traditionalCellArray,s(1));
                        innerLevel = subsasgn(innerLevel,s(2:end),varargin{:});
                        obj.traditionalCellArray = builtin('subsasgn',obj.traditionalCellArray,s(1),innerLevel);
                    end
                case '()'
                    obj.traditionalCellArray = builtin('subsasgn',obj.traditionalCellArray,s,varargin{:});
                otherwise
                    error('Not a valid indexing expression')    
            end
        end
        
        function n = length(obj)
            n = length(obj.traditionalCellArray);
        end
        
        function n = numArgumentsFromSubscript(~,s,~)            
            n = 1;
            for i=1:length(s)
                if strcmp(s(i).type,'{}') || strcmp(s(i).type,'()')
                    n = n * prod(cellfun(@length,s(i).subs));
                end
            end            
        end
                        
    end
    
    methods (Static)
        function locSubsAsgn()
        end
    end
end

