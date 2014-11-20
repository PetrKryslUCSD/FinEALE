% stores mhash object as a persistent variable and manages all its methods
% *************************************************************************
% function ArgOut = HashtableStore(action, ind, varargin)
%
% Purpose: this function stores mhash object as a persistent variable and
%          manages all its methods - maps keys to values, store and 
%          retrieve objects from a hashtable, ecc. The methods are
%          equivalent to the Java methods.
%
%          For the hashtable methods description see:
%          http://java.sun.com/j2se/1.4.2/docs/api/java/util/Hashtable.html
%
% Author:
%     Mikhail Poda-Razumtsev
%
% Initial Version:
%     25.11.2005
%
% Change History:
%
% *************************************************************************
function ArgOut = HashtableStore(action, ind, varargin)

ArgOut = [];

persistent HTable
if isempty(HTable), HTable = local_newElement; end

switch lower(action)
    
% Java methods
    case 'clear'
        HTable = local_Clear(HTable, ind);
    case 'clone'
        [HTable, ArgOut] = local_Clone(HTable, ind);
    case 'contains'
        ArgOut = local_ContainsValue(HTable, ind, varargin{:});
    case 'containskey'
        ArgOut = local_ContainsKey(HTable, ind, varargin{:});
    case 'containsvalue'
        ArgOut = local_ContainsValue(HTable, ind, varargin{:});
    case 'elements'
        ArgOut = local_Elements(HTable, ind);
    case 'entryset'
        % do nothing
    case 'equals'
        ArgOut = local_Equals(HTable, ind, varargin{:});
    case 'get'
        [HTable, ArgOut] = local_Get(HTable, ind, varargin{:});
    case 'hashcode'
        ArgOut = ind;
    case 'isempty'
        ArgOut = local_IsEmpty(HTable, ind);
    case 'keys'
        ArgOut = local_Keys(HTable, ind);
    case 'keyset'
        ArgOut = local_Keys(HTable, ind);
    case 'put'
        [HTable, ArgOut] = local_Put(HTable, ind, varargin{:});
    case 'putall'
        HTable = local_PutAll(HTable, ind, varargin{:});
    case 'rehash'
        % do nothing
    case 'remove'
        HTable = local_Remove(HTable, ind, varargin{:});
    case 'size'
        ArgOut = local_Size(HTable, ind);
    case 'tostring'
        ArgOut = local_ToString(HTable, ind);
    case 'values'
        ArgOut = local_Elements(HTable, ind);
                
% Non-Java methods
    case 'new'
        [HTable, ArgOut] = local_New(HTable);
    case 'delete'
        [HTable, ArgOut] = local_removeElement(HTable, ind);
    case 'createcopy'
        [HTable, ArgOut] = local_CreateCopy(HTable, ind);
    case 'display'
        ArgOut = local_Display(HTable);
    case 'isinternid'
        ArgOut = local_IsInternID(ind);
    case 'ismhashtable'
        ArgOut = local_isKey(HTable, ind);
    case 'getelement'
        ArgOut = local_getElement(HTable, ind);
        
    otherwise
        % do nothing
end


% *************************************************************************
% function HTable = local_Clear(HTable, ind)
%
% Purpose: Clears this hashtable so that it contains no keys.
%
% *************************************************************************
function HTable = local_Clear(HTable, ind)

element = local_newElement;
HTable  = local_putElement(HTable, ind, element);


% *************************************************************************
% function [HTable, HClone] = local_Clone(HTable, ind)
%
% Purpose: Creates a shallow copy of this hashtable. Uses intern ID to make
%           sure that no outside call could make the same thing
%
% Status: DIRTY
%
% *************************************************************************
function [HTable, HClone] = local_Clone(HTable, ind)

[HTable, ind]   = local_CreateCopy(HTable, ind);
HClone          = mhashtable(ind, 849784.87374/exp(pi));


% *************************************************************************
% function isContained = local_ContainsValue(HTable, ind, value)
%
% Purpose: Tests if some key maps into the specified value in this hashtable.
%
% *************************************************************************
function isContained = local_ContainsValue(HTable, ind, value)

isContained = 0;
element     = local_getElement(HTable, ind);
if isempty(element), return; end

data        = {element.num_data{:}, element.char_data{:}};

if ischar(value)
    for i = 1:length(data)
        val = data{i};
        if ischar(val)
            if strcmp(val, value)
                isContained = 1;
                return
            end
        end
    end
elseif isnumber(value)
    for i = 1:length(data)
        val = data{i};
        if isnumber(val)
            if val == value
                isContained = 1;
                return
            end
        end
    end
else
    % not yet supported match for non char/num values
    isContained = [];
    return;
end


% *************************************************************************
% function isKey = local_ContainsKey(HTable, ind, value)
%
% Purpose: Tests if the specified object is a key in this hashtable.
%
% *************************************************************************
function isKey = local_ContainsKey(HTable, ind, key)

isKey       = 0;
element     = local_getElement(HTable, ind);
if isempty(element), return; end

isKey       = local_isKey(element, key);


% *************************************************************************
% function isKey = local_isKey(element, key)
%
% Purpose: Tests if the specified key exists in this hashtable.
%
% *************************************************************************
function isKey = local_isKey(element, key)

isKey       = 0;
index       = [];

if ischar(key)
    index = find(strcmp(element.char_keys, key));
elseif isnumeric(key)
    index = find(element.num_keys == key(1));
else
    return
end

if ~isempty(index), isKey = 1; end


% *************************************************************************
% function data = local_Elements(HTable, ind)
%
% Purpose: Returns an enumeration of the values in this hashtable.
%
% *************************************************************************
function data = local_Elements(HTable, ind)

data    = [];
element	= local_getElement(HTable, ind);
if isempty(element), return; end

data  	= {element.num_data{:}, element.char_data{:}};


% *************************************************************************
% function isEqual = local_Equals(HTable, ind, obj)
%
% Purpose: Compares the specified Object with this hashtable for equality
%
% *************************************************************************
function isEqual = local_Equals(HTable, ind, obj)

isEqual = 0;
if ~isa(obj, 'mhashtable'), return; end
obj = struct(obj);
if ind == obj.ind, isEqual = 1; end


% *************************************************************************
% function [HTable, val] = local_Get(HTable, ind, key)
%
% Purpose: Returns the value to which the specified key is mapped in this 
%          hashtable. If no such key is existent the return an empty array
%
% *************************************************************************
function [HTable, val] = local_Get(HTable, ind, key)

val     = [];
element = local_getElement(HTable, ind);
if isempty(element), return; end

val     = local_getElement(element, key);


% *************************************************************************
% function isHEmpty = local_IsEmpty(HTable, ind)
%
% Purpose: Tests if this hashtable maps no keys to values.
%
% *************************************************************************
function isHEmpty = local_IsEmpty(HTable, ind)

isHEmpty = 1;
element  = local_getElement(HTable, ind);
if isempty(element), return; end
if ~(isempty(element.num_data) && isempty(element.char_data)), isHEmpty = 0; end


% *************************************************************************
% function keys = local_Keys(HTable, ind)
%
% Purpose: Returns an enumeration of the keys in this hashtable.
%
% *************************************************************************
function keys = local_Keys(HTable, ind)

keys        = [];
element     = local_getElement(HTable, ind);
if isempty(element), return; end

num_keys    = num2cell(element.num_keys);
keys        = {num_keys{:}, element.char_keys{:}};


% *************************************************************************
% function [HTable, success] = local_Put(HTable, ind, key, val)
%
% Purpose: Maps the specified key to the specified value in this hashtable.
%
% *************************************************************************
function [HTable, success] = local_Put(HTable, ind, key, val)

success             = 0;
element             = local_getElement(HTable, ind);
if isempty(element), return; end

[element, success]  = local_putElement(element, key, val);
HTable              = local_putElement(HTable, ind, element);


% *************************************************************************
% function HTable = local_PutAll(HTable, ind, obj)
%
% Purpose: Copies all of the mappings from the specified obj to this
%           Hashtable These mappings will replace any mappings that this 
%           Hashtable had for any of the currently specified keys
%
% *************************************************************************
function HTable = local_PutAll(HTable, ind, obj)

if ~isa(obj, 'hashtable'), return; end

element     = local_getElement(HTable, ind);
if isempty(element), return; end

obj_data  	= obj.getElement;
if isempty(obj_data), return; end

num_data    = obj_data.num_data;
char_data   = obj_data.char_data;

num_keys    = obj_data.num_keys;
char_keys   = obj_data.char_keys;

% loop for all numeric keys
for i = 1:length(num_keys)
    element = local_putElement(element, num_keys(i), num_data{i});
end

% loop for all char keys
for i = 1:length(char_keys)
    element = local_putElement(element, char_keys{i}, num_data{i});
end

HTable  = local_putElement(HTable, ind, element);


% *************************************************************************
% function HTable = local_Remove(HTable, ind, key)
%
% Purpose: Removes the key (and its corresponding value) from this hashtable.
%
% *************************************************************************
function HTable = local_Remove(HTable, ind, key)

element = local_getElement(HTable, ind);
if isempty(element), return; end

element = local_removeElement(element, key);
HTable  = local_putElement(HTable, ind, element);


% *************************************************************************
% function HSize = local_Size(HTable, ind)
%
% Purpose: Returns the number of keys in this hashtable.
%
% *************************************************************************
function HSize = local_Size(HTable, ind)

HSize   = [];
element = local_getElement(HTable, ind);
if isempty(element), return; end

HSize   = length(element.char_keys) + length(element.num_keys);


% *************************************************************************
% function data = local_ToString(HTable, ind)
%
% Purpose: Returns a string representation of this Hashtable object in the
%           form of a set of entries
%
% *************************************************************************
function dataString = local_ToString(HTable, ind)

try
dataString  = '[]';
element     = local_getElement(HTable, ind);
if isempty(element), return; end

num_keys    = num2cell(element.num_keys);
char_keys   = element.char_keys;

char_data   = element.char_data;
num_data    = element.num_data;

data        = {num_data{:}, char_data{:}};
keys        = {num_keys{:}, char_keys{:}};

data_length = length(data);

if data_length > 0
    
    char_ones   = ones(data_length, 1);
    dpoints     = char(58*char_ones);
    empty       = char(32*char_ones);
    
    data_cell   = cell(1, data_length);
    key_cell    = data_cell;
    
    % loop for all data and key elements
    for i = 1:data_length
        data_cell{i}    = local_DisplayCell(data(i));
        cell_str        = local_DisplayCell(keys(i));
        key_cell{i}     = cell_str(end:-1:1);
    end
    
    data_str    = char(data_cell);
    key_str     = char(key_cell);
        
    dataString  = [key_str(:,end:-1:1), dpoints, empty, data_str];
end
catch
dataString = 'an error has occured while producing string rapresentation of hashtable contents';    
end


% *************************************************************************
% function cell_cell = local_DisplayCell(cell_cell, cell_data, ind, isInverse)
%
% Purpose: returns string view of the cell contents
%
% *************************************************************************
function cell_str = local_DisplayCell(cell_cell)

iCell       = cell_cell;
iElement    = cell_cell{1};

if ischar(iElement)
    if length(iElement) > 25
        cell_str = sprintf('[%0.0fx%0.0f char]', size(iElement));
    else
        cell_str = sprintf('''%s''', iElement);
    end
elseif isnumeric(iElement)
    if length(iElement) == 1
        cell_str = sprintf('%g', iElement);
    elseif length(iElement) > 5
        cell_str = sprintf('[%0.0fx%0.0f %s]', size(iElement), class(iElement(1)));
    else
        cell_str = hsprintf(iElement);
    end
else
    cell_str = evalc('disp(iCell)');
    cell_str([1, 2, 3, 4, end - 1, end]) = [];
end


% *************************************************************************
% function [HTable, ind] = local_New(HTable)
%
% Purpose: creates a new hashtable object
%
% *************************************************************************
function [HTable, ind] = local_New(HTable)

element = local_newElement;
ind     = max(HTable.num_keys) + 1;
if isempty(ind), ind = 1; end
HTable  = local_putElement(HTable, ind, element);


% *************************************************************************
% function [element, success] = local_removeElement(element, key)
%
% Purpose: removes hashtable element
%
% *************************************************************************
function [element, success] = local_removeElement(element, key)

success = 0;

if ischar(key)
    index = find(strcmp(element.char_keys, key));
    if ~isempty(index)
        element.char_data(index)    = []; 
        element.char_keys(index)    = []; 
        success                     = 1;
    end
elseif isnumeric(key)
    index = find(element.num_keys == key(1));
    if ~isempty(index)
        element.num_data(index)     = [];
        element.num_keys(index)     = []; 
        success                     = 1;
    end
else
    return
end


% *************************************************************************
% function [HTable, ind] = local_CreateCopy(HTable, ind)
%
% Purpose: creates a copy of the element and returns it index
%
% *************************************************************************
function [HTable, ind] = local_CreateCopy(HTable, ind)

ind     = [];
element = local_getElement(HTable, ind);
if isempty(element), return; end

ind     = max(HTable.num_keys) + 1;
if isempty(ind), ind = 1; end
HTable  = local_putElement(HTable, ind, element);


% *************************************************************************
% function isIntern = local_IsInternID(internID)
%
% Purpose: verifies whether the argument matches the intern ID constant
%
% Status: DIRTY
%
% *************************************************************************
function isIntern = local_IsInternID(internID)

isIntern = 0;
if ~isnumeric(internID) || length(internID) > 1, return; end
if internID == 849784.87374/exp(pi), isIntern = 1; end


% *************************************************************************
% function element = local_newElement
%
% Purpose: returns a struct corresponding to an empty element
%
% *************************************************************************
function element = local_newElement

element.num_data    = {};
element.char_data   = {};
element.num_keys    = [];
element.char_keys   = [];


% *************************************************************************
% function [val, success] = local_getElement(element, key)
%
% Purpose: returns the value corresponding to the given key
%
% *************************************************************************
function [val, success] = local_getElement(element, key)

success = 0;
val     = [];
index   = [];

if isempty(element) || isempty(key), return; end

if ischar(key)
    index = find(strcmp(element.char_keys, key));
    if ~isempty(index), val = element.char_data{index}; end
elseif isnumeric(key)
    index = find(element.num_keys == key(1));
    if ~isempty(index), val = element.num_data{index}; end
else
    return
end

if ~isempty(index), success = 1; end


% *************************************************************************
% function [element, success] = local_putElement(element, key, val)
%
% Purpose: assigns value to the given key
%
% *************************************************************************
function [element, success] = local_putElement(element, key, val)

success	= 0;
if isempty(element) || isempty(key), return; end

if ischar(key)
    index = find(strcmp(element.char_keys, key));
    if isempty(index)
        index = length(element.char_keys) + 1;
        element.char_keys{index} = key;
        element.char_data{index} = val;
    else
        element.char_data{index} = val;
    end
    success	= 1;
elseif isnumeric(key)
    index = find(element.num_keys == key(1));
    if isempty(index)
        index = length(element.num_keys) + 1;
        element.num_keys(index) = key(1);
        element.num_data{index} = val;
    else
        element.num_data{index} = val;         
    end
    success	= 1;
end