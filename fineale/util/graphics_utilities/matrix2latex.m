%  Print a matrix in LaTeX tabular format.
%
% function s = matrix2latex(M,options)
%
%   matrix2latex(M) prints out the numeric matrix M in a LaTeX tabular
%   format. The '&' character appears between entries in a row, '\\'
%   is appended to the ends of rows.
%
% precision,separator,symbol,Zerotolerance,displaystyle,medskip
%   Copyright 2008,2010 by Petr Krysl
%
function s = matrix2latex(M,options)
    if (~exist('options'))
        options= [];
    end
    if (~isfield(options,'precision'))
        precision = 2;
    else
        precision = options.precision;
    end
    if (~isfield(options,'separator'))
        separator = '';
    else
        separator = options.separator;
    end
    if (~isfield(options,'symbol'))
        symbol = [];
    else
        symbol = options.symbol;
    end
    if (~isfield(options,'Zerotolerance'))
        Zerotolerance = eps;
    else
        Zerotolerance = options.Zerotolerance;
    end
    if (~isfield(options,'displaystyle'))
        displaystyle = '';
    else
        displaystyle = '\displaystyle{}';
    end
    if (~isfield(options,'medskip'))
        medskip = '';
    else
        medskip = '\medskip';
    end
    
    
    if ndims(M) > 2
        error('Works only for 2D arrays.')
    end
    
    
    
    % Extend the format specifiers.
    [m,n] = size(M);
    
    S=['\begin{array}{' sprintf('%c',repmat('c',1,n)) '}'   ];
    
    for i=1:m
        for j=1:n-1
            S=[S displaystyle str(M(i,j))  separator ' & ' ];
        end
        S=[S  str(M(i,n)) medskip ' \\' ];
    end
    
    S=[S '\end{array}' sprintf('\n')];
    
    % Display or output?
    if nargout==0
        disp(S)
    else
        s = S;
    end
    
    function String =str(Element)
        if (isempty(symbol))
            if isa(Element,'double')
                String =num2str(Element, precision);
            else
                String =texlabel(char(Element));
            end
        else
            if (abs(Element)>Zerotolerance)
                String =symbol;
            else
                String ='';
            end
        end
    end
    
    
end

function T = texlabel(varargin)
    %TEXLABEL Produces the TeX format from a character string.
    %   TEXLABEL(f) converts the expression f into the TeX equivalent
    %   for title/label application.  It processes transliterated Greek
    %   variables to print (in titles/labels) as actual Greek letters.
    %
    %   TEXLABEL(f,'literal') prints the literal label.
    %
    %   If the title/label is too long to fit into a plot window, then
    %   the center of the expression is removed and an ellipsis ...
    %   is inserted.
    %
    %   TEXLABEL is used in EZSURF, EZMESH, etc. to generate TeX format
    %   for the title, x-, y-, and z-labels for these plots.
    %
    %   Examples:
    %    texlabel('sin(sqrt(x^2 + y^2))/sqrt(x^2 + y^2)')
    %          returns
    %    {sin}({sqrt}({x}^{2} + {y}^{2}))/{sqrt}({x}^{2} + {y}^{2})
    %
    %    texlabel(['3*(1-x)^2*exp(-(x^2) - (y+1)^2) - 10*(x/5 - x^3 - y^5)*' ...
    %         'exp(-x^2-y^2) - 1/3*exp(-(x+1)^2 - y^2)'])
    %          returns
    %    {3} ({1}-{x})^{2} {exp}(-({x}^{2}) - ({y}+{1})^{2}) -...- {1}/{3} {exp}(-({x}+{1})^{2} - {y}^{2})
    %
    %    texlabel('lambda12^(3/2)/pi - pi*delta^(2/3)')
    %          returns
    %    {\lambda_{12}}^{{3}/{2}}/{\pi} - {\pi} {\delta}^{{2}/{3}}
    %
    %    texlabel('lambda12^(3/2)/pi - pi*delta^(2/3)','literal')
    %          returns
    %    {lambda12}^{{3}/{2}}/{pi} - {pi} {delta}^{{2}/{3}}
    
    %   Copyright 1984-2008 The MathWorks, Inc.
    %   $Revision: 1.12.4.8 $  $Date: 2008/05/01 20:14:24 $
    
    narginchk(1,2);% Required for 2016
    switch nargin
        % texlabel(f)
        case 1
            flag = 1;
            % texlabel(f,'literal')
        case 2
            flag = varargin{2};
            % texlabel(f,1,'literal') or texlabel(f,'literal',1)
    end
    
    % Change the title to TeX format.
    T = char(varargin{1});
    if ~isequal(flag,'literal') && (flag ~= 1)
        len = flag;
        flag = 1;
    else
        len = 49;
    end
    T = titlelen(T,len);
    Lp = findstr(T,'^(');
    % For each instance of the pattern "^(", find the corresponding ")" for it
    % (assuming it exists).
    rightParenInds = [];
    len = numel(T);
    for i=Lp
        balance = 1;
        currInd = i+1;
        while balance~=0 && currInd <= len
            currInd = currInd+1;
            if T(currInd) == ')'
                balance = balance-1;
            elseif T(currInd) == '('
                balance = balance+1;
            end
        end
        if currInd <= len
            rightParenInds(end+1) = currInd; %#ok<NASGU>
        end
    end
    % If there are no right parentheses, don't do any replacement
    if ~isempty(rightParenInds)
        T = strrep(T,'^(','^{');
        T(rightParenInds) = '}';
    end
    % If the string is of the form x^10 or 2^alpha, then
    % change the string to {x}^{10}, etc.
    % Find the indices where T is alpha-numeric.
    N = isletter(T) | (T >= '0' & T <= '9');
    % Vector of right and left brackets.
    bracket = '{}';
    % Flag for 1 or 0.
    lookfor = 1;
    % Parse the string to include right and left brackets around
    % all alpha-numeric characters.
    
    % Special case: If the value is numeric, but the value directly to its left
    % is a ".", the treat the "." as part of the number:
    for k=2:length(N)
        if N(k) == lookfor && T(k-1) == '.'
            N(k-1) = true;
        end
    end
    for k=length(N):-1:1
        if N(k)==lookfor
            T = [T(1:k), bracket(lookfor+1), T(k+1:end)];
            lookfor = ~lookfor;
        end
    end
    % Place a right bracket at the beginning of the string
    % if necessary.
    if (lookfor == 0)
        T = ['{',T];
    end
    % Remove multiplication sign .* and *:  x*y -> x y.
    T = strrep(T,'.*',' ');
    T = strrep(T,'*',' ');
    % Replace .^ with ^ and ./ with /
    T = strrep(T,'.^','^');
    T = strrep(T,'./','/');
    if flag == 1
        % Change to Greek and subscripted variables.
        T = greeks(T);
    else
        % Allow for literal x_1 versus TeX version.
        T = strrep(T,'_','\_');
    end
end
%---------------------------------

function L = titlelen(f,len)
    %TITLELEN Reduces the title length to fit into the title space
    %   of length len.
    %   TITLELEN(f,len) takes the char f and returns a sym of length
    %   len to fit into the title space.
    
%     if ~isempty(get(0,'CurrentFigure')) && ~isempty(get(gcf,'CurrentAxes'))
%         ax = gca;
%         fig = ancestor(ax,'figure');
%         p = hgconvertunits(fig,get(ax,'Position'),get(ax,'Units'),...
%             'characters',get(ax,'Parent'));
%         tlen = max(p(3)*.65,len);
%     else
%         tlen = 71;
%     end
%     
%     if length(f) > tlen
%         L = slen(f,tlen);
%     else
        L = f;
        % end 
end

%---------------------------------

function L = greeks(L)
    %GREEKS Changes English words to Greek letters.
    %  GREEKS('alpha')  returns '\alpha' so that it will
    %  appear as Greek letter in a title.
    
    G = {'alpha','beta','gamma','delta','epsilon','zeta', ...
        'eta','theta','iota','kappa','lambda','mu','nu', ...
        'xi','pi','rho','sigma','tau','upsilon', ...
        'phi','chi','psi','omega',...
        'Gamma','Delta', ...
        'Theta','Lambda', ...
        'Xi','Pi','Sigma','Upsilon', ...
        'Phi','Psi','Omega'};
    % Replace _ with + to isolate all possible variable names to
    % return {C}_{\alpha} in the case of texlabel('C_alpha').
    V = symvar(strrep(L,'_','+'));
    % Find the common elements of V and G;
    [VG,IV] = intersect(V,G);
    if ~isequal(VG,{''})
        % Replace gamma with \gamma, etc. in L
        for j = 1:length(VG)
            VG{j} = ['\' VG{j}];
            L = strrep(L,V{IV(j)},VG{j});
        end
    end
    % Now replace gamma1 with \gamma_1, etc.
    for j = 1:length(V)
        for k = 1:length(G)
            % Determine whether V{j} contains a Greek symbol
            C{j,k} = strmatch(G{k},V{j});
            VG{j,k} = isequal(V{j},G{k});
            % If V{j} contains a Greek symbol and is NOT identical
            % to a Greek symbol, then place a \ before V{j} and
            % _{N} about its trailing section.
            if ~isempty(C{j,k}) && ~VG{j,k}
                num = find( (V{j} >= '0' & V{j} <= '9') == 1);
                if ~isempty(num)
                    L = strrep(L,V{j},['\' V{j}(1:num(1)-1) ...
                        '_{' V{j}(num(1):end) '}']);
                end
            end
        end
        % Now make x1 look like x_1, etc.
        if isempty(intersect(V,G)) && ~isequal(V{j},'atan2')
            num = find( (V{j} >= '0' & V{j} <= '9') == 1);
            if ~isempty(num)
                L = strrep(L,V{j},[V{j}(1:num(1)-1) '_{' V{j}(num(1):end) '}']);
            end
        end
    end
    
    % Replace '{pi}' by '{\pi}'
    L = strrep(L,'{pi}','{\pi}');
end

%----------------------------------%
function S = slen(S,len)
    % SLEN Reduces the length of a mathematical expression by removing
    %   terms until the string is less than LEN characters long.
    %
    %   Example:
    %     S = ['1+x+1/2*x^2-1/8*x^4-1/15*x^5-1/240*x^6+1/90*x^7+' ...
    %          '31/5760*x^8+1/5670*x^9-2951/3628800*x^10-1/3150*x^11'];
    %     slen(S)  returns
    %         1+x+1/2*x^2-1/8*x^4-1/15*x^5-1/240*x^6+1/90*x^7+...-1/3150*x^11
    
    B = plmin(S); % Determine where S has + or -.
    B1 = find(B == 1);
    S1 = S;
    for j = 1:length(B1)-1
        if length(S) < len
            return;
        else
            S = [S(1:B1(end-j)) '...' S1(B1(end):end)];
        end
    end
end

%-------------------------

function c = plmin(s)
    %PLMIN  PLMIN(s) is true for + or - not inside parentheses.
    p = cumsum((s == '(') - (s == ')'));
    c = (s == '+' | s == '-') & (p == 0);
end
