classdef  Abaqus_exporter  < handle
    %     Simple class to export an abaqus INP file.
    %
    % %     Example:
    % % In this example no part, instance, assembly are defined.
    % AE=Abaqus_exporter;
    % AE.open(['f-' num2str(ref) '.inp']);
    % AE.HEADING([mfilename ' ' 'ElType=' ElType ' ' 'ref=' num2str(ref)]);
    % AE.NODE(model_data.geom.values);
    % AE.ELEMENT(ElType,'All',1,fes.conn);
    % AE.ELEMENT('SFM3D4','Traction',count(fes)+1,model_data.boundary_conditions.traction{1}.fes.conn);
    % AE.NSET_NSET('Corner',nc);
    % AE.ORIENTATION('Orientation', axis, aangle);
    % AE.MATERIAL('Material');
    % AE.ELASTIC(E1,E2,E3,nu12,nu13,nu23,G12,G13,G23);
    % AE.SOLID_SECTION('Material','Orientation','All');
    % AE.SURFACE_SECTION('Traction');
    % AE.STEP_PERTURBATION('Linear solution');
    % AE.BOUNDARY_direct(model_data.u.is_fixed,model_data.u.fixed_values);
    % AE.DLOAD('Traction',model_data.boundary_conditions.traction{1}.traction);
    % AE.NODE_PRINT('Corner');
    % AE.END_STEP();
    % AE.close();
    % system(['abaqus job=' ['f-' num2str(ref)]]);
    %
    %     Example:
    %      In this example we define a part, and an instance within an assembly.
    %     AE=Abaqus_exporter;
    %     AE.open(['f-' num2str(ref) '.inp']);
    %     AE.HEADING([mfilename ' ' 'ElType=' ElType ' ' 'ref=' num2str(ref)]);
    %     AE.PART('part1');
    %     AE.END_PART();
    %     AE.ASSEMBLY('ASSEM1');
    %     AE.INSTANCE('INSTNC1','PART1');
    %     AE.NODE(model_data.geom.values);
    %     AE.ELEMENT(ElType,'All',1,fes.conn);
    %     AE.ELEMENT(SurfElType,'Traction',count(fes)+1,model_data.boundary_conditions.traction{1}.fes.conn);
    %     AE.NSET_NSET('Corner',nc);
    %     AE.ORIENTATION('Orientation', Rm(:,1), Rm(:,2));
    %     AE.SOLID_SECTION('Material','Orientation','All');
    %     AE.SURFACE_SECTION('Traction');
    %     AE.NSET_NSET('xfix',find(model_data.u.is_fixed(:,1)));
    %     AE.NSET_NSET('yfix',find(model_data.u.is_fixed(:,2)));
    %     AE.NSET_NSET('zfix',find(model_data.u.is_fixed(:,3)));
    %     AE.END_INSTANCE();
    %     AE.END_ASSEMBLY();
    %     AE.MATERIAL('Material');
    %     AE.ELASTIC(E1,E2,E3,nu12,nu13,nu23,G12,G13,G23);
    %     AE.STEP_PERTURBATION('Linear solution');
    %     AE.DLOAD('ASSEM1.INSTNC1.Traction',model_data.boundary_conditions.traction{1}.traction);
    %     AE.BOUNDARY('ASSEM1.INSTNC1.xfix',1);
    %     AE.BOUNDARY('ASSEM1.INSTNC1.yfix',2);
    %     AE.BOUNDARY('ASSEM1.INSTNC1.zfix',3);
    %     AE.NODE_PRINT('ASSEM1.INSTNC1.Corner');
    %     AE.ENERGY_PRINT();
    %     AE.END_STEP();
    %     AE.close();
    
    
    properties
        filename='f';
        fid=[];
        status = [];
    end
    
    methods % constructor
        
        function self = Abaqus_exporter(Parameters)
            if nargin <1
                Parameters = struct([]);
            end
            
        end
        
    end
    
    
    methods % concrete methods
        
        function open(self, File)
            % Open INP file File.
            [pathstr,name,ext]= fileparts(File);
            self.filename = name;
            self.fid=fopen([self.filename '.inp'],'w+');
        end
        
        function close(self)
            % Close current file.
            fclose(self.fid);
        end
        
        function HEADING(self, Text)
            % Write out the *HEADING option.
            % Text is an arbitrary string.
            self.status =fprintf(self.fid,'%s\n', ['*HEADING']);
            self.status =fprintf(self.fid,'%s\n', Text);
        end
        
        function PART(self, NAME)
            % Write out the *PART option.
            % NAME is an arbitrary string.
            self.status =fprintf(self.fid,'%s\n', ['*PART, NAME=' NAME]);
        end
        
        function END_PART(self)
            % Write out the *END PART option.
            self.status =fprintf(self.fid,'%s\n', ['*END PART']);
        end
        
        function ASSEMBLY(self, NAME)
            % Write out the *ASSEMBLY option.
            % NAME is an arbitrary string.
            self.status =fprintf(self.fid,'%s\n', ['*ASSEMBLY, NAME=' NAME]);
        end
        
        function END_ASSEMBLY(self)
            % Write out the *END ASSEMBLY option.
            self.status =fprintf(self.fid,'%s\n', ['*END ASSEMBLY']);
        end
        
        function INSTANCE(self, NAME, PART)
            % Write out the *INSTANCE option.
            % Text is an arbitrary string.
            self.status =fprintf(self.fid,'%s\n', ['*INSTANCE, NAME=' NAME ', PART=' PART]);
        end
        
        function END_INSTANCE(self)
            % Write out the *END INSTANCE option.
            self.status =fprintf(self.fid,'%s\n', ['*END INSTANCE']);
        end
        
        function NODE(self,xyz)
            % Write out the *NODE option.
            %    xyz=array of node coordinates
            self.status =fprintf(self.fid,'%s\n', ['*NODE']);
            for j=1:size(xyz,1)
                self.status =fprintf(self.fid,'%g,%g,%g,%g\n', j,xyz(j,:));
            end
        end
        
        function ELEMENT(self,TYPE,ELSET,start,conn)
            % Write out the *ELEMENT option.
            %    TYPE= element type code,
            %    ELSET= element set to which the elements belong,
            %    start= start the element numbering at this integer,
            %    conn= connectivity array for the elements, one row per element
            self.status =fprintf(self.fid,'%s\n', ['*ELEMENT, TYPE =' TYPE ', ELSET=' ELSET]);
            nn9=size(conn,2); nn8=min  ([nn9,15]);
            for j=1:size(conn,1)
                self.status =fprintf(self.fid,'%g', j+start-1);
                self.status =fprintf(self.fid,',%g', conn(j,1:nn8));
                if (nn9>15)
                    self.status =fprintf(self.fid,',\n');
                    self.status =fprintf(self.fid,'%g', conn(j,15+1));
                    self.status =fprintf(self.fid,',%g', conn(j,15+1+1:end));
                end
                self.status =fprintf(self.fid,'\n');
            end
        end
        
        function NSET_NSET(self,NSET,n)
            % Write out the *NSET option.
            self.status =fprintf(self.fid,'%s\n', ['*NSET, NSET=' NSET]);
            for j=1:length(n)
                self.status =fprintf(self.fid,'%g,\n', n(j));
            end
        end
        
        function ELSET_ELSET(self,ELSET,n)
            % Write out the *ELSET option.
            self.status =fprintf(self.fid,'%s\n', ['*ELSET, ELSET=' ELSET]);
            for j=1:length(n)
                self.status =fprintf(self.fid,'%g,\n', n(j));
            end
        end
        
        function ORIENTATION(self, ORIENTATION, a, b)
            % Write out the *ORIENTATION option.
            self.status =fprintf(self.fid,'%s\n', ['*ORIENTATION,NAME=' ORIENTATION]);
            self.status =fprintf(self.fid,'%g,%g,%g,%g,%g,%g\n',a,b);
            self.status =fprintf(self.fid,'%g,%g\n', 1, 0);
        end
        
        function MATERIAL(self, MATERIAL)
            % Write out the *MATERIAL option.
            self.status =fprintf(self.fid,'%s\n', ['*MATERIAL,NAME=' MATERIAL]);
        end
        
        function ELASTIC(self, E1,E2,E3,nu12,nu13,nu23,G12,G13,G23)
            % Write out the *ELASTIC option.
            self.status =fprintf(self.fid,'%s\n', ['*ELASTIC,TYPE=ENGINEERING CONSTANTS ']);
            self.status =fprintf(self.fid,'%g,', E1,E2,E3,nu12,nu13,nu23,G12,G13);
            self.status =fprintf(self.fid,'\n%g,\n', G23);
        end
        
        function ELASTIC_ISOTROPIC(self, E,nu)
            % Write out the *ELASTIC option.
            self.status =fprintf(self.fid,'%s\n', ['*ELASTIC,TYPE=ISOTROPIC ']);
            self.status =fprintf(self.fid,'%g,%g\n', E,nu);
        end
        
        function DENSITY(self, rho)
            % Write out the *DENSITY option.
            self.status =fprintf(self.fid,'%s\n', ['*DENSITY ']);
            self.status =fprintf(self.fid,'%g\n', rho);
        end
        
        function SECTION_CONTROLS(self, NAME,WHAT)
            % Write out the *SECTION CONTROLS option.
            %    MATERIAL= material name,
            %    ORIENTATION= orientation name
            %    ELSET=  element set to which  the section applies
            self.status =fprintf(self.fid,'%s\n', ['*SECTION CONTROLS, NAME=' NAME ',' WHAT]);
        end
        
        function SOLID_SECTION(self, MATERIAL,ORIENTATION,ELSET,CONTROLS)
            % Write out the *SOLID SECTION option.
            %    MATERIAL= material name,
            %    ORIENTATION= orientation name
            %    ELSET=  element set to which  the section applies
            if (~exist ('CONTROLS', 'var'))
                self.status =fprintf(self.fid,'%s\n', ['*SOLID SECTION,MATERIAL=' MATERIAL ',ORIENTATION =' ORIENTATION ',ELSET=' ELSET]);
            else
                self.status =fprintf(self.fid,'%s\n', ['*SOLID SECTION,MATERIAL=' MATERIAL ',ORIENTATION =' ORIENTATION ',ELSET=' ELSET ', CONTROLS  =' CONTROLS]);
            end
        end
        
        function HOURGLASS(self, WHICH,VALUE)
            % Write out the *HOURGLASS option.
            %             Example:
            % *SOLID SECTION,ELSET=SOLID3,MATERIAL=MAT,CONTROL=A
            % *HOURGLASS STIFFNESS
            % 5.E8
            
            self.status =fprintf(self.fid,'%s\n', ['*HOURGLASS ' WHICH]);
            self.status =fprintf(self.fid,'%g\n', VALUE);
        end
        
        function SURFACE_SECTION(self, ELSET)
            % Write out the *SURFACE SECTION option.
            %    ELSET=  element set to which  the section applies
            self.status =fprintf(self.fid,'%s\n', ['*SURFACE SECTION, ELSET=' ELSET]);
        end
        
        
        function STEP_PERTURBATION(self, Text)
            % Write out the *STEP,PERTURBATION option.
            self.status =fprintf(self.fid,'%s\n', '*STEP,PERTURBATION');
            self.status =fprintf(self.fid,'%s\n', '*STATIC');
        end
        
        function STEP_PERTURBATION_STATIC(self, Text)
            % Write out the *STEP,PERTURBATION option for linear static analysis.
            self.status =fprintf(self.fid,'%s\n', '*STEP,PERTURBATION');
            self.status =fprintf(self.fid,'%s\n', '*STATIC');
        end
        
        function STEP_PERTURBATION_BUCKLE(self, Text, neigv)
            % Write out the *STEP,PERTURBATION option for linear buckling analysis.
            self.status =fprintf(self.fid,'%s\n', '*STEP, name=Buckling, nlgeom=NO, perturbation');
            self.status =fprintf(self.fid,'%s\n', '*BUCKLE');
            self.status =fprintf(self.fid,'%d, , , , \n', neigv);
        end
        
        function STEP_FREQUENCY(self, Nmodes)
            % Write out the *STEP,FREQUENCY option.
            self.status =fprintf(self.fid,'%s\n', '*STEP');
            self.status =fprintf(self.fid,'%s\n', '*FREQUENCY, EIGENSOLVER=LANCZOS');
            self.status =fprintf(self.fid,'%g, , ,-1.E6 \n', Nmodes);
        end
        
        function BOUNDARY_direct(self,is_fixed,fixed_value)
            % Write out the *BOUNDARY option.
            %    is_fixed= array of Boolean flags (true means fixed, or prescribed),  one row per node,
            %    fixed_value=array of displacements to which the corresponding displacement components is fixed
            self.status =fprintf(self.fid,'%s\n', '*BOUNDARY');
            for j=1:size(is_fixed,1)
                for k=1:size(is_fixed,2)
                    %<node number>, <first dof>, <last dof>, <magnitude of displacement>
                    if (is_fixed(j,k))
                        self.status =fprintf(self.fid,'%g,%g,%g,%g\n', j,k,k,fixed_value(j,k));
                    end
                end
            end
        end
        
        function BOUNDARY(self,NSET,dof)
            % Write out the *BOUNDARY option to fix displacements.
            %    NSET= node set,
            %    dof=Degree of freedom, 1, 2, 3
            self.status =fprintf(self.fid,'%s\n', '*BOUNDARY');
            self.status =fprintf(self.fid,'%s,%g\n', NSET,dof);
        end
        
        
        function BOUNDARY_type(self,NSET,dof,value)
            % Write out the *BOUNDARY option to prescribed nonzero displacements.
            %    NSET= node set,
            %    dof=Degree of freedom, 1, 2, 3
            self.status =fprintf(self.fid,'%s\n', '*BOUNDARY,TYPE=DISPLACEMENT');
            self.status =fprintf(self.fid,'%s,%g,%g,%g\n', NSET,dof,dof,value);
        end
        
        function DLOAD(self,ELSET,traction)
            % Write out the *DLOAD option.
            %    ELSET= element set to which the traction is applied,
            %    traction=traction vector
            self.status =fprintf(self.fid,'%s\n', '*DLOAD, follower=NO');
            self.status =fprintf(self.fid,'%s,%s,%g,%g,%g,%g\n', ELSET,'TRVEC',norm(traction),traction/norm(traction));
        end
        
        
        function CLOAD(self,NSET,dof, magnitude)
            % Write out the *CLOAD option.
            %    NSET=Number of node
            %    dof= 1, 2, 3,
            %    magnitude= signed multiplier
            self.status =fprintf(self.fid,'%s\n', '*CLOAD');
            self.status =fprintf(self.fid,'%s,%g,%g\n', NSET,dof, magnitude);
        end
        
        function NODE_PRINT(self,NSET)
            % Write out the *NODE PRINT option.
            self.status =fprintf(self.fid,'%s\n', ['*NODE PRINT, NSET=' NSET]);
            self.status =fprintf(self.fid,'U\n');
        end
        
        function EL_PRINT(self,ELSET, KEYS)
            % Write out the *EL PRINT option.
            self.status =fprintf(self.fid,'%s\n', ['*EL PRINT, ELSET=' ELSET ', POSITION=INTEGRATION POINTS, SUMMARY= YES']);
            self.status =fprintf(self.fid,'\n');
            self.status =fprintf(self.fid,'%s\n', KEYS);
        end
        
        function ENERGY_PRINT(self)
            % Write out the *ENERGY PRINT option.
            self.status =fprintf(self.fid,'%s\n', ['*ENERGY PRINT']);
        end
        
        function END_STEP(self)
            % Write out the *END STEP option.
            self.status =fprintf(self.fid,'%s\n', '*END STEP');
        end
        
        
        function A= extract_energy_from_abaqus_dat(File,Tagline)
            A=[];
            try
                fid=fopen(File,'r');
                while true
                    temp = fgetl(fid);   if ~ischar(temp), break, end
                    temp=deblank(fliplr(deblank(temp(end:-1:1))));
                    if (strcmpi(temp,Tagline))
                        for j=1:3
                            temp = fgetl(fid);
                        end
                        temp = fgetl(fid);
                        A = sscanf(temp(40:end), '%g');
                        fclose(fid);
                        return;
                    end
                end
            catch,fclose(fid);end
        end
        
        function d= extract_displacement_from_abaqus_dat(File,Tagline,nnodes)
            if (~exist('nnodes','var') ),nnodes=1;end
            d=zeros(nnodes,3);
            try
                fid=fopen(File,'r');
                while true
                    temp = fgetl(fid);   if ~ischar(temp), break, end
                    temp=deblank(fliplr(deblank(temp(end:-1:1))));
                    if (strcmpi(temp,Tagline))
                        for j=1:4
                            temp = fgetl(fid);
                        end
                        for j=1:nnodes
                            temp = fgetl(fid);
                            A = sscanf(temp, '%g %g %g %g');
                            d(j,:)=A(2:end);
                        end
                        fclose(fid);
                        return;
                    end
                end
            catch,fclose(fid);end
        end
        
        function d= extract_buckling_from_abaqus_dat(File,Tagline,neigv)
            if (~exist('neigv','var') ),neigv=1;end
            d=zeros(neigv,1);
            try
                fid=fopen(File,'r');
                while true
                    temp = fgetl(fid);   if ~ischar(temp), break, end
                    temp=deblank(fliplr(deblank(temp(end:-1:1))));
                    if (strcmpi(temp,Tagline))
                        for j=1:2
                            temp = fgetl(fid);
                        end
                        for j=1:neigv
                            temp = fgetl(fid);
                            A = sscanf(temp, '%g %g');
                            d(j,:)=A(2:end);
                        end
                        fclose(fid);
                        return;
                    end
                end
            catch,fclose(fid);end
        end
        
        
        
    end
    
end


