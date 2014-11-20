classdef  physical_units_machine
    % Class for evaluation of measurement units.
    %
    %
    % Examples: 
    %
    %     m=physical_units_machine; %
    %     m('kg/mm^3')
    %     mu_mud=30*(1/100)*mu('GM/CM/SEC');% dynamic viscosity, centiPoise
    %
    % See also: physical_units_struct
        
        
    properties (SetAccess= private)
        system_of_units,base_time_units
    end
    
    properties ( Hidden)
        u % structure with definitions of the units of measurement
    end
    
    methods
        
        function  self=physical_units_machine(Parameters)
        %  Constructor
        % Parameters:
        %    system_of_units='SI';
        %    base_time_units='SEC';
        % See discussion of constructors in <a href="matlab:helpwin 'fineale/classes/Contents'">classes/Contents</a>.
        
        self.system_of_units='SI';
        self.base_time_units='SEC';
        self.u = physical_units_struct(self.system_of_units,self.base_time_units);
        self.u.S=self.u.SEC;% Alias for a second
        self.u.N=self.u.NT;% The symbol for the Newton
        self.u.MPA=self.u.MEGA*self.u.PA;% The symbol for the mega Pascal
        self.u.GPA=self.u.GIGA*self.u.PA;% The symbol for the giga Pascal
        self.u.KIPS=self.u.KILO*self.u.LBF;% The symbol for the thousands of pounds force
        self.u.KSI=self.u.KILO*self.u.PSI;% The symbol for the thousands of pounds per square inch
        self.u.CP=1/1000*self.u.PA*self.u.SEC;% The symbol for the centiPoise
        if (nargin<1),return,end
        if (isfield(Parameters,'system_of_units'))
            self.system_of_units=Parameters.system_of_units;
        end
        if (isfield(Parameters,'base_time_units'))
            self.base_time_units=Parameters.base_time_units;
        end
        end
        
        function fact = subsref(self,s)
        % Implement a special subscripted assignment
        switch s(1).type
            case '()'
                fact=factor(self,s(1).subs{1});
            case '.'
                switch s(1).subs
                    case 'system_of_units'
                        fact = self.system_of_units;
                    case 'base_time_units'
                        fact = self.base_time_units;
                    otherwise
                        
                end
            otherwise
                error('Specify value for x as obj(x)')
        end
        end
        
        function  fact=factor(self,expr)
        % Compute the conversion factor for the units in the expression.
        outstr=self.replace_symbols(expr,'self.u');
        try
          fact=eval(outstr);
        catch
          error (['Invalid units specification: ' expr]);
        end
        end
        
    end
    
    methods( Access = private )
        
        function outstr=replace_symbols(self,str,objectname)
        str=upper(str);
        outstr='';
        i=1;
        while i<=length(str)
            if isletter(str(i))
                k=i+1;
                while (k<=length(str))&&(isletter(str(k)))
                    k=k+1;
                end
                outstr=[outstr,'(',objectname,'.',str(i:k-1),')'];
                i=k-1;
            else
                outstr=[outstr,str(i)];
            end
            i=i+1;
        end
        
        end
        
    end
    
end
