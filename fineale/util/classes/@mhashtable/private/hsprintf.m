% gets the string from a number the way you would like to see it
% *************************************************************************
% function val_string = hsprintf(varargin)
%
% Purpose: gets the string from a number for the best human viewing, there
%					 is no unnecessary information (digits), but also not more information
%					 can be lost lost then specified by the SignificantDigitsNo attribute.
%					 It extends sprintf to the functionality that any very large or very
%				   small number could be displayed the way you expect and would like.
%
% Arguments:
%       varargin : see example below. Attribute-value pairs in arbitrary order.
%									 All attributes have a default value and can be not explicitly
%									 specified
%
% Return values
%     val_string : string
%
% Author:
%     Mikhail Poda-Razumtsev
%
% Initial Version:
%     26.11.2005
%
% Examples:
%
%   1. With default settings:
%           str = hsprintf(859)
%
%   2. Array input with default settings
%           str = hsprintf([859, 798])
%           >> str = 859, 798
%
%   3. MATLAB - sprintf behaviour
%           str = hsprintf('%5.6f', 798)
%
%   4. Set all parameters
%           str = hsprintf('Value', [798, pi], 'CharsNo', 8, 'SignificantDigitsNo', 4, 'Separator', '; ')
%           >> str = 798; 3.142
%
% *************************************************************************
function val_string = hsprintf(varargin)

val_string  = [];
ArgIn       = {};

% resolve different cases: with one input, or C conversion mode
if length(varargin) == 1
    ArgIn = {'Value', varargin{1}};
else
    ind_format = strfind(varargin{1}, '%');

    % special case, sprintf syntax with empty format
    if isempty(varargin{1})
        ArgIn = {'Value', varargin{2}};

    % sprintf syntax
    elseif ~isempty(ind_format)
        val_string = sprintf(varargin{:});
        return;

    % dafault hsprintf syntax
    else
        ArgIn = varargin;
    end
end

% set properties
Value               = [];
CharsNo             = 8;
SignificantDigitsNo	= 4;
Separator           = ', ';

for i = 1:2:length(ArgIn)
    switch lower(ArgIn{i})
        case 'value'
            Value = ArgIn{i+1};
        case 'CharsNo'
            CharsNo = ArgIn{i+1};
        case 'SignificantDigitsNo'
            SignificantDigitsNo = ArgIn{i+1};
        case 'Separator'
            Separator = ArgIn{i+1};
    end
end

if (isempty(Value) || ~isnumeric(Value)), return; end

CharsNo                 = max(4, round(CharsNo));
SignificantDigitsNo     = min(CharsNo - 1, round(SignificantDigitsNo));

% skalar
if length(Value) == 1
    val_string = local_getValString(Value, CharsNo, SignificantDigitsNo);

% vector
else
    val_cell = cell(1, length(Value));
    for iValue = 1:length(Value)
        val = Value(iValue);
        val_cell{iValue} = local_getValString(val, CharsNo, SignificantDigitsNo);
    end
    val_string = local_cell2string(val_cell, Separator);
end


% *************************************************************************
% function val_string = local_getValString(Value, CharsNo, SignificantDigitsNo)
%
% Purpose:
%     function to
%
% Arguments:
%     Value:
%   CharsNo:
%   SignificantDigitsNo:
%
% Return values
%     val_string
%
% *************************************************************************
function val_string = local_getValString(Value, CharsNo, SignificantDigitsNo)

val_string  = [];
realVal     = real(Value);
imagVal     = imag(Value);

if realVal ~= 0
    val_string  = local_resolve(realVal, CharsNo, SignificantDigitsNo);
end

if imagVal ~= 0
    imag_string = local_resolve(imagVal, CharsNo, SignificantDigitsNo);
    if imagVal < 0
        val_string  = [val_string, imag_string, 'i'];
    else
        val_string  = [val_string, '+', imag_string, 'i'];
    end
end

if isempty(val_string), val_string = '[]'; end


% *************************************************************************
% function String = local_cell2string(Cell, Delimiter)
%
% Purpose:
%     function to
%
% Arguments:
%       Cell:
%  Delimiter:
%
% Return values
%     String:
%
% *************************************************************************
function String = local_cell2string(Cell, Delimiter)

String = Cell{1};

for iCell = 2:length(Cell)
    String = [String, Delimiter, Cell{iCell}];
end


% *************************************************************************
% function val_string = local_resolve(Value, CharsNo, SignificantDigitsNo)
%
% Purpose:
%     function to
%
% Arguments:
%     Value:
%   CharsNo:
%   SignificantDigitsNo:
%
% Return values
%     val_string
%
% *************************************************************************
function val_string = local_resolve(Value, CharsNo, SignificantDigitsNo)

Value = double(Value);

if Value < 0
    isSign	= 1;
    Value 	= -Value;
elseif Value > 0
    isSign	= 0;
else
    val_string = '0';
    return
end

isSmall     = Value < 1;

val_log10  	= log10(Value);
decNo       = floor(val_log10);
val         = Value;
val_norm    = val/10^(decNo);


isExp       = (isSmall && (abs(decNo) >= (CharsNo - SignificantDigitsNo))) || (~isSmall && decNo > CharsNo - 3);


% cut all values beyond the accuracy given by CharsNo
if ~isExp
    if isSmall
        digNo = SignificantDigitsNo - 1;
    else
        if SignificantDigitsNo > decNo
            digNo = SignificantDigitsNo - 1;
        else
            digNo = decNo;
        end
    end
else
    digNo = SignificantDigitsNo;
end

val_norm = round(val_norm*10^(digNo))/(10^(digNo));

% pre-allocate array
val_array   = zeros(1, CharsNo);

% loop for all digits
for iPlace = 1:CharsNo
    val_floor           = floor(val_norm);
    val_rem             = val_norm - val_floor;
    val_norm            = val_rem*10;

    isNeg               =  (1 - val_rem) < 1e-5;
    isPos               =  val_rem < 1e-5;

    if isNeg
        val_array(iPlace) = val_floor + 1;
    else
        val_array(iPlace) 	= val_floor;
    end

    if isNeg || isPos
        val_array(iPlace+1:end) = [];
        break
    end
end

val_string    = sprintf('%d', val_array);

if ~isExp
    if decNo >= 0
        if length(val_string) > (decNo + 1)
            val_string  = [val_string(1:decNo + 1), '.', val_string(decNo + 2:end)];
        else
            val_string  = [val_string, sprintf('%d', zeros((decNo + 1 - length(val_string)), 1))];
        end
    else
        val_string  = ['0.', sprintf('%d', zeros(-decNo-1, 1)), val_string];
    end
else
    if length(val_string) > 1
        val_string = [val_string(1), '.', val_string(2:end), 'e', sprintf('%+d', decNo)];
    else
        val_string = [val_string(1), 'e', sprintf('%+d', decNo)];
    end
end

if isSign, val_string = ['-', val_string]; end
