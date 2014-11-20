% constructor for the mhashtable class
% *************************************************************************
% function mhash = mhashtable(varargin)
%
% Purpose: returns an object of the mhashtable class. The object is a
%          REFERENCE!!!
%
%          When called with no arguments the reference to a new hashtable is
%          returned. When called with the mhashtable object as argument, the
%          copy of this obejct is made and its reference is returned.
%
%          For the hashtable methods description see:
%          http://java.sun.com/j2se/1.4.2/docs/api/java/util/Hashtable.html
%
% Author:
%     Mikhail Poda-Razumtsev
%
% Initial Version:
%     26.11.2005
%
% Change History:
%
% *************************************************************************
function mhash = mhashtable(varargin)

try
h = [];

if nargin == 0
    h.ind = HashtableStore('new');
else
    obj = varargin{1};
    if isa(obj, 'mhashtable')
        h.ind = HashtableStore('createCopy', obj.ind);
        
    % private call, tried to minimize the probability of an accidental call from outside
    elseif HashtableStore('isInternID', varargin{2})
        h.ind = obj;
    end
end
    
if ~isempty(h), mhash = class(h, 'mhashtable'); end

catch
error('HashTable:mhashtable', 'Invalid number of arguments.');
end