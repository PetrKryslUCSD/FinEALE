classdef  fe_set_H64 < fe_set_3_manifold
% H64 (hexahedron with 64 nodes)  finite element (FE) set class.
%
%

    
    properties (Constant, GetAccess = public)
    end
    
    methods % constructor
        
        function self = fe_set_H64(Parameters)
            % Constructor.
            % Parameters:
            %           same as fe_set_3_manifold.
            % See discussion of constructors in <a href="matlab:helpwin 'fineale/classes/Contents'">classes/Contents</a>.
            if nargin <1
                Parameters.nfens=64;
            end
            Parameters.nfens=64;
            self = self@fe_set_3_manifold(Parameters);
        end
        
    end
    
    
    methods % concrete methods
        
        function N = bfun(self,param_coords)
            % Evaluate the basis function matrix for an 27-node brick.
            %
            % function val = bfun(self,param_coords)
            %
            %   Call as:
            %      N = bfun (g, pc)
            %   where
            %      g=FE set,
            %      pc=parametric coordinates, -1 < pc(j) < 1, length(pc)=3.
            %
            xi=param_coords(1);
            et=param_coords(2);
            ze=param_coords(3);
            
            A =1/4096; B =9/4096; C=81/4096; D =729/4096;
            e1 =(9*xi^2-1); e2 =(9*et^2-1); e3 =(9*ze^2-1);
            e4 =(xi^2-1); e5 =(et^2-1); e6 =(ze^2-1);
            e7 =(3*xi-1); e8 = (3*et-1); e9 =(3*ze-1);
            e10 =(3*xi+1); e11 = (3*et+1); e12 =(3*ze+1);
            e13=(xi+1); e14 =(xi-1); e15 =(et+1); e16 =(et-1); e17 =(ze+1); e18 =(ze-1);
            N=[-A*e18*e14*e16*e3*e1*e2;
                B*e18*e7*e16*e3*e4*e2;
                -B*e18*e10*e16*e3*e4*e2;
                A*e18*e13*e16*e3*e1*e2;
                B*e18*e14*e8*e3*e1*e5;
                -C*e18*e7*e8*e3*e4*e5;
                C*e18*e10*e8*e3*e4*e5;
                -B*e18*e13*e8*e3*e1*e5;
                -B*e18*e14*e11*e3*e1*e5;
                C*e18*e7*e11*e3*e4*e5;
                -C*e18*e10*e11*e3*e4*e5;
                B*e18*e13*e11*e3*e1*e5;
                A*e18*e14*e15*e3*e1*e2;
                -B*e18*e7*e15*e3*e4*e2;
                B*e18*e10*e15*e3*e4*e2;
                -A*e18*e13*e15*e3*e1*e2;
                B*e9*e14*e16*e6*e1*e2;
                -C*e9*e7*e16*e6*e4*e2;
                C*e9*e10*e16*e6*e4*e2;
                -B*e9*e13*e16*e6*e1*e2;
                -C*e9*e14*e8*e6*e1*e5;
                D*e9*e7*e8*e6*e4*e5;
                -D*e9*e10*e8*e6*e4*e5;
                C*e9*e13*e8*e6*e1*e5;
                C*e9*e14*e11*e6*e1*e5;
                -D*e9*e7*e11*e6*e4*e5;
                D*e9*e10*e11*e6*e4*e5;
                -C*e9*e13*e11*e6*e1*e5;
                -B*e9*e14*e15*e6*e1*e2;
                C*e9*e7*e15*e6*e4*e2;
                -C*e9*e10*e15*e6*e4*e2;
                B*e9*e13*e15*e6*e1*e2;
                -B*e12*e14*e16*e6*e1*e2;
                C*e12*e7*e16*e6*e4*e2;
                -C*e12*e10*e16*e6*e4*e2;
                B*e12*e13*e16*e6*e1*e2;
                C*e12*e14*e8*e6*e1*e5;
                -D*e12*e7*e8*e6*e4*e5;
                D*e12*e10*e8*e6*e4*e5;
                -C*e12*e13*e8*e6*e1*e5;
                -C*e12*e14*e11*e6*e1*e5;
                D*e12*e7*e11*e6*e4*e5;
                -D*e12*e10*e11*e6*e4*e5;
                C*e12*e13*e11*e6*e1*e5;
                B*e12*e14*e15*e6*e1*e2;
                -C*e12*e7*e15*e6*e4*e2;
                C*e12*e10*e15*e6*e4*e2;
                -B*e12*e13*e15*e6*e1*e2;
                A*e17*e14*e16*e3*e1*e2;
                -B*e17*e7*e16*e3*e4*e2;
                B*e17*e10*e16*e3*e4*e2;
                -A*e17*e13*e16*e3*e1*e2;
                -B*e17*e14*e8*e3*e1*e5;
                C*e17*e7*e8*e3*e4*e5;
                -C*e17*e10*e8*e3*e4*e5;
                B*e17*e13*e8*e3*e1*e5;
                B*e17*e14*e11*e3*e1*e5;
                -C*e17*e7*e11*e3*e4*e5;
                C*e17*e10*e11*e3*e4*e5;
                -B*e17*e13*e11*e3*e1*e5;
                -A*e17*e14*e15*e3*e1*e2;
                B*e17*e7*e15*e3*e4*e2;
                -B*e17*e10*e15*e3*e4*e2;
                A*e17*e13*e15*e3*e1*e2];
            
            %     % Symbolic derivation of the matrix of basis function values
            %     x=sym('x','real');
            %     xi=sym('xi','real');
            %     et=sym('et','real');
            %     ze=sym('ze','real');
            %     a=1/3;
            %     x0 =-1;
            %     l1 =eval(((x-1)*(x-a)*(x-(-a)))/((x0-1)*(x0-a)*(x0-(-a))));
            %     x0 =-a;
            %     l2 =eval(((x-1)*(x-a)*(x-(-1)))/((x0-1)*(x0-a)*(x0-(-1))));
            %     x0 =+a;
            %     l3 =eval(((x-1)*(x-(-a))*(x-(-1)))/((x0-1)*(x0-(-a))*(x0-(-1))));
            %     x0 =+1;
            %     l4 =eval(((x-a)*(x-(-a))*(x-(-1)))/((x0-a)*(x0-(-a))*(x0-(-1))));
            %     subs(l1,'x','et')
            %     xis = [subs(l1,'x','xi'),subs(l2,'x','xi'),subs(l3,'x','xi'),subs(l4,'x','xi')];
            %     ets = [subs(l1,'x','et'),subs(l2,'x','et'),subs(l3,'x','et'),subs(l4,'x','et')];
            %     zes = [subs(l1,'x','ze'),subs(l2,'x','ze'),subs(l3,'x','ze'),subs(l4,'x','ze')];
            %     xisets=(xis'*ets);
            %     N(:,:,1)=xisets*zes(1);
            %     N(:,:,2)=xisets*zes(2);
            %     N(:,:,3)=xisets*zes(3);
            %     N(:,:,4)=xisets*zes(4);
            %     for ze=[-1,-a,+a,+1]
            %         for et=[-1,-a,+a,+1]
            %             for xi=[-1,-a,+a,+1]
            %                 n=double(subs(subs(subs(N,'ze',num2str(ze)),'et',num2str(et)),'xi',num2str(xi)));
            %                 i=find(n>0.99);
            %                 disp([ char(simple(N(i))) ';'])
            %             end
            %         end
            %     end
            %
            %     x=sym('x','real');
            %     xi=sym('xi','real');
            %     et=sym('et','real');
            %     ze=sym('ze','real');
            %     a=1/3;
            %     x0 =-1;
            %     l1 =eval(((x-1)*(x-a)*(x-(-a)))/((x0-1)*(x0-a)*(x0-(-a))));
            %     x0 =-a;
            %     l2 =eval(((x-1)*(x-a)*(x-(-1)))/((x0-1)*(x0-a)*(x0-(-1))));
            %     x0 =+a;
            %     l3 =eval(((x-1)*(x-(-a))*(x-(-1)))/((x0-1)*(x0-(-a))*(x0-(-1))));
            %     x0 =+1;
            %     l4 =eval(((x-a)*(x-(-a))*(x-(-1)))/((x0-a)*(x0-(-a))*(x0-(-1))));
            %     subs(l1,'x','et')
            %     xis = [subs(l1,'x','xi'),subs(l2,'x','xi'),subs(l3,'x','xi'),subs(l4,'x','xi')];
            %     ets = [subs(l1,'x','et'),subs(l2,'x','et'),subs(l3,'x','et'),subs(l4,'x','et')];
            %     zes = [subs(l1,'x','ze'),subs(l2,'x','ze'),subs(l3,'x','ze'),subs(l4,'x','ze')];
            %     xisets=(xis'*ets);
            %     N(:,:,1)=xisets*zes(1);
            %     N(:,:,2)=xisets*zes(2);
            %     N(:,:,3)=xisets*zes(3);
            %     N(:,:,4)=xisets*zes(4);
            %     for ze=[-1,-a,+a,+1]
            %         for et=[-1,-a,+a,+1]
            %             for xi=[-1,-a,+a,+1]
            %                 n=double(subs(subs(subs(N,'ze',num2str(ze)),'et',num2str(et)),'xi',num2str(xi)));
            %                 i=find(n>0.99);
            %                 disp(['N(' num2str(i) ',:)=[' char(simple(diff(N(i),'xi'))) ',' char(simple(diff(N(i),'et'))) ',' char(simple(diff(N(i),'ze'))) '];'])
            %             end
            %         end
            %     end
            
        end
        
        
        function N = bfundpar (self, param_coords)
            % Evaluate the derivatives of the basis function matrix.
            %
            % function Nder = bfundpar (self, param_coords)
            %
            % Returns an array of NFENS rows, and DIM columns, where
            %    NFENS=number of nodes, and
            %    DIM=number of spatial dimensions.
            % Call as:
            %    Nder = bfundpar(g, pc)
            % where g=FE set
            %       pc=parametric coordinates, -1 < pc(j) < 1, length(pc)=3.
            %
            xi=param_coords(1);
            et=param_coords(2);
            ze=param_coords(3);
            A =1/4096; B =9/4096; C=81/4096; D =729/4096;
            xi2=xi^2; et2=et^2; ze2=ze^2;
            e1 =(9*xi2-1); e2 =(9*et2-1); e3 =(9*ze2-1);
            e4 =(xi2-1); e5 =(et2-1); e6 =(ze2-1);
            e7 =(3*xi-1); e8 = (3*et-1); e9 =(3*ze-1);
            e10 =(3*xi+1); e11 = (3*et+1); e12 =(3*ze+1);
            e13=(xi+1); e14 =(xi-1); e15 =(et+1); e16 =(et-1); e17 =(ze+1); e18 =(ze-1);
            e20=(-1+27*xi2-18*xi); e21=(-1+27*et2-18*et); e22=(27*ze2-18*ze-1);
            e23=(9*xi2-2*xi-3); e24=(9*xi2+2*xi-3); e25=(27*xi2-1+18*xi);
            e26=(-1-18*xi+27*xi2); e27=(-3+9*xi2-2*xi); e28=(-1+27*et2+18*et);
            e29=(-3+9*et2+2*et); e30=(27*et2-1-18*et); e31=(-1-18*et+27*et2);
            e32=(-1+27*ze2-18*ze); e33=(-3+9*xi2+2*xi); e34=(9*et2-2*et-3);
            e35=(9*et2-3-2*et); e36=(-3+9*et2-2*et); e37=(-1+27*xi2+18*xi);
            e38=(-18*xi+27*xi2-1); e39=(-3+2*et+9*et2); e40=(2*et-3+9*et2);
            e41=(9*et2+2*et-3); e42=(27*xi2-18*xi-1); e43=(27*et2-1+18*et);
            e44=(9*xi2-3+2*xi); e45=(27*xi2-1-18*xi); e46=(-1+18*et+27*et2);
            e47=(9*ze2-3-2*ze); e48=(9*ze2-2*ze-3); e49=(27*ze2+18*ze-1);
            e50=(-1+27*ze2+18*ze);
            N(1,:)=[-A*e18*e9*e12*e16*e8*e11*e20,-A*e18*e9*e12*e21*e14*e7*e10,-A*e22*e16*e8*e11*e14*e7*e10];
            N(2,:)=[B*e18*e9*e12*e16*e8*e11*e23,B*e18*e9*e12*e30*e14*e7*e13,B*e22*e16*e8*e11*e14*e7*e13];
            N(3,:)=[-B*e18*e9*e12*e16*e8*e11*e24,-B*e18*e9*e12*e31*e14*e10*e13,-B*e22*e16*e8*e11*e14*e10*e13];
            N(4,:)=[A*e18*e9*e12*e16*e8*e11*e25,A*e18*e9*e12*e31*e7*e10*e13,A*e32*e16*e8*e11*e7*e10*e13];
            N(5,:)=[B*e18*e9*e12*e16*e8*e15*e26,B*e18*e9*e12*e34*e14*e7*e10,B*e32*e16*e8*e15*e14*e7*e10];
            N(6,:)=[-C*e18*e9*e12*e16*e8*e15*e27,-C*e18*e9*e12*e35*e14*e7*e13,-C*e22*e16*e8*e15*e14*e7*e13];
            N(7,:)=[C*e18*e9*e12*e16*e8*e15*e33,C*e18*e9*e12*e36*e14*e10*e13,C*e22*e16*e8*e15*e14*e10*e13];
            N(8,:)=[-B*e18*e9*e12*e16*e8*e15*e37,-B*e18*e9*e12*e36*e7*e10*e13,-B*e22*e16*e8*e15*e7*e10*e13];
            N(9,:)=[-B*e18*e9*e12*e16*e11*e15*e38,-B*e18*e9*e12*e39*e14*e7*e10,-B*e22*e16*e11*e15*e14*e7*e10];
            N(10,:)=[C*e18*e9*e12*e16*e11*e15*e27,C*e18*e9*e12*e29*e14*e7*e13,C*e22*e16*e11*e15*e14*e7*e13];
            N(11,:)=[-C*e18*e9*e12*e16*e11*e15*e33,-C*e18*e9*e12*e40*e14*e10*e13,-C*e22*e16*e11*e15*e14*e10*e13];
            N(12,:)=[B*e18*e9*e12*e16*e11*e15*e37,B*e18*e9*e12*e41*e7*e10*e13,B*e22*e16*e11*e15*e7*e10*e13];
            N(13,:)=[A*e18*e9*e12*e8*e11*e15*e42,A*e18*e9*e12*e43*e14*e7*e10,A*(-1-18*ze+27*ze2)*e8*e11*e15*e14*e7*e10];
            N(14,:)=[-B*e18*e9*e12*e8*e11*e15*e23,-B*e18*e9*e12*e28*e14*e7*e13,-B*(-1-18*ze+27*ze2)*e8*e11*e15*e14*e7*e13];
            N(15,:)=[B*e18*e9*e12*e8*e11*e15*e33,B*e18*e9*e12*(27*et2+18*et-1)*e14*e10*e13,B*e22*e8*e11*e15*e14*e10*e13];
            N(16,:)=[-A*e18*e9*e12*e8*e11*e15*e37,-A*e18*e9*e12*e46*e7*e10*e13,-A*e22*e8*e11*e15*e7*e10*e13];
            N(17,:)=[B*e18*e9*e17*e16*e8*e11*e45,B*e18*e9*e17*e30*e14*e7*e10,B*e47*e16*e8*e11*e14*e7*e10];
            N(18,:)=[-C*e18*e9*e17*e16*e8*e11*e27,-C*e18*e9*e17*e30*e14*e7*e13,-C*e47*e16*e8*e11*e14*e7*e13];
            N(19,:)=[C*e18*e9*e17*e16*e8*e11*e33,C*e18*e9*e17*e30*e14*e10*e13,C*e47*e16*e8*e11*e14*e10*e13];
            N(20,:)=[-B*e18*e9*e17*e16*e8*e11*(27*xi2+18*xi-1),-B*e18*e9*e17*(-18*et+27*et2-1)*e7*e10*e13,-B*e48*e16*e8*e11*e7*e10*e13];
            N(21,:)=[-C*e18*e9*e17*e16*e8*e15*(-18*xi-1+27*xi2),-C*e18*e9*e17*e35*e14*e7*e10,-C*e48*e16*e8*e15*e14*e7*e10];
            N(22,:)=[D*e18*e9*e17*e16*e8*e15*(-3-2*xi+9*xi2),D*e18*e9*e17*e36*e14*e7*e13,D*e48*e16*e8*e15*e14*e7*e13];
            N(23,:)=[-D*e18*e9*e17*e16*e8*e15*(-3+2*xi+9*xi2),-D*e18*e9*e17*(-2*et+9*et2-3)*e14*e10*e13,-D*e48*e16*e8*e15*e14*e10*e13];
            N(24,:)=[C*e18*e9*e17*e16*e8*e15*(-1+18*xi+27*xi2),C*e18*e9*e17*(-2*et+9*et2-3)*e7*e10*e13,C*e48*e16*e8*e15*e7*e10*e13];
            N(25,:)=[C*e18*e9*e17*e16*e11*e15*e26,C*e18*e9*e17*e29*e14*e7*e10,C*e48*e16*e11*e15*e14*e7*e10];
            N(26,:)=[-D*e18*e9*e17*e16*e11*e15*(-2*xi-3+9*xi2),-D*e18*e9*e17*e40*e14*e7*e13,-D*e48*e16*e11*e15*e14*e7*e13];
            N(27,:)=[D*e18*e9*e17*e16*e11*e15*e33,D*e18*e9*e17*e40*e14*e10*e13,D*(-2*ze+9*ze2-3)*e16*e11*e15*e14*e10*e13];
            N(28,:)=[-C*e18*e9*e17*e16*e11*e15*(27*xi2+18*xi-1),-C*e18*e9*e17*e41*e7*e10*e13,-C*e48*e16*e11*e15*e7*e10*e13];
            N(29,:)=[-B*e18*e9*e17*e8*e11*e15*e45,-B*e18*e9*e17*(18*et-1+27*et2)*e14*e7*e10,-B*e48*e8*e11*e15*e14*e7*e10];
            N(30,:)=[C*e18*e9*e17*e8*e11*e15*(-2*xi+9*xi2-3),C*e18*e9*e17*e46*e14*e7*e13,C*e48*e8*e11*e15*e14*e7*e13];
            N(31,:)=[-C*e18*e9*e17*e8*e11*e15*e44,-C*e18*e9*e17*(18*et-1+27*et2)*e14*e10*e13,-C*e48*e8*e11*e15*e14*e10*e13];
            N(32,:)=[B*e18*e9*e17*e8*e11*e15*e25,B*e18*e9*e17*e28*e7*e10*e13,B*e48*e8*e11*e15*e7*e10*e13];
            N(33,:)=[-B*e18*e12*e17*e16*e8*e11*e42,-B*e18*e12*e17*(27*et2-18*et-1)*e14*e7*e10,-B*(9*ze2+2*ze-3)*e16*e8*e11*e14*e7*e10];
            N(34,:)=[C*e18*e12*e17*e16*e8*e11*(-3-2*xi+9*xi2),C*e18*e12*e17*(27*et2-18*et-1)*e14*e7*e13,C*(2*ze+9*ze2-3)*e16*e8*e11*e14*e7*e13];
            N(35,:)=[-C*e18*e12*e17*e16*e8*e11*(-3+2*xi+9*xi2),-C*e18*e12*e17*(-18*et+27*et2-1)*e14*e10*e13,-C*(2*ze+9*ze2-3)*e16*e8*e11*e14*e10*e13];
            N(36,:)=[B*e18*e12*e17*e16*e8*e11*(27*xi2+18*xi-1),B*e18*e12*e17*(-18*et+27*et2-1)*e7*e10*e13,B*(2*ze+9*ze2-3)*e16*e8*e11*e7*e10*e13];
            N(37,:)=[C*e18*e12*e17*e16*e8*e15*e38,C*e18*e12*e17*e34*e14*e7*e10,C*(9*ze2+2*ze-3)*e16*e8*e15*e14*e7*e10];
            N(38,:)=[-D*e18*e12*e17*e16*e8*e15*(-2*xi-3+9*xi2),-D*e18*e12*e17*(-3-2*et+9*et2)*e14*e7*e13,-D*(9*ze2+2*ze-3)*e16*e8*e15*e14*e7*e13];
            N(39,:)=[D*e18*e12*e17*e16*e8*e15*e33,D*e18*e12*e17*e35*e14*e10*e13,D*(9*ze2+2*ze-3)*e16*e8*e15*e14*e10*e13];
            N(40,:)=[-C*e18*e12*e17*e16*e8*e15*e37,-C*e18*e12*e17*(-2*et-3+9*et2)*e7*e10*e13,-C*(9*ze2+2*ze-3)*e16*e8*e15*e7*e10*e13];
            N(41,:)=[-C*e18*e12*e17*e16*e11*e15*e20,-C*e18*e12*e17*e40*e14*e7*e10,-C*(9*ze2-3+2*ze)*e16*e11*e15*e14*e7*e10];
            N(42,:)=[D*e18*e12*e17*e16*e11*e15*e27,D*e18*e12*e17*e29*e14*e7*e13,D*(9*ze2+2*ze-3)*e16*e11*e15*e14*e7*e13];
            N(43,:)=[-D*e18*e12*e17*e16*e11*e15*(2*xi+9*xi2-3),-D*e18*e12*e17*e29*e14*e10*e13,-D*(9*ze2+2*ze-3)*e16*e11*e15*e14*e10*e13];
            N(44,:)=[C*e18*e12*e17*e16*e11*e15*(-1+18*xi+27*xi2),C*e18*e12*e17*e29*e7*e10*e13,C*(9*ze2+2*ze-3)*e16*e11*e15*e7*e10*e13];
            N(45,:)=[B*e18*e12*e17*e8*e11*e15*e42,B*e18*e12*e17*e46*e14*e7*e10,B*(9*ze2+2*ze-3)*e8*e11*e15*e14*e7*e10];
            N(46,:)=[-C*e18*e12*e17*e8*e11*e15*(-3-2*xi+9*xi2),-C*e18*e12*e17*e28*e14*e7*e13,-C*(9*ze2+2*ze-3)*e8*e11*e15*e14*e7*e13];
            N(47,:)=[C*e18*e12*e17*e8*e11*e15*e44,C*e18*e12*e17*(27*et2+18*et-1)*e14*e10*e13,C*(9*ze2+2*ze-3)*e8*e11*e15*e14*e10*e13];
            N(48,:)=[-B*e18*e12*e17*e8*e11*e15*e25,-B*e18*e12*e17*e28*e7*e10*e13,-B*(9*ze2+2*ze-3)*e8*e11*e15*e7*e10*e13];
            N(49,:)=[A*e9*e12*e17*e16*e8*e11*e26,A*e9*e12*e17*(27*et2-18*et-1)*e14*e7*e10,A*e49*e16*e8*e11*e14*e7*e10];
            N(50,:)=[-B*e9*e12*e17*e16*e8*e11*e27,-B*e9*e12*e17*(27*et2-18*et-1)*e14*e7*e13,-B*(18*ze-1+27*ze2)*e16*e8*e11*e14*e7*e13];
            N(51,:)=[B*e9*e12*e17*e16*e8*e11*e44,B*e9*e12*e17*(-18*et+27*et2-1)*e14*e10*e13,B*(18*ze-1+27*ze2)*e16*e8*e11*e14*e10*e13];
            N(52,:)=[-A*e9*e12*e17*e16*e8*e11*(18*xi+27*xi2-1),-A*e9*e12*e17*e21*e7*e10*e13,-A*e49*e16*e8*e11*e7*e10*e13];
            N(53,:)=[-B*e9*e12*e17*e16*e8*e15*e38,-B*e9*e12*e17*e36*e14*e7*e10,-B*e49*e16*e8*e15*e14*e7*e10];
            N(54,:)=[C*e9*e12*e17*e16*e8*e15*(-2*xi+9*xi2-3),C*e9*e12*e17*e36*e14*e7*e13,C*(-1+18*ze+27*ze2)*e16*e8*e15*e14*e7*e13];
            N(55,:)=[-C*e9*e12*e17*e16*e8*e15*(2*xi-3+9*xi2),-C*e9*e12*e17*e34*e14*e10*e13,-C*(-1+18*ze+27*ze2)*e16*e8*e15*e14*e10*e13];
            N(56,:)=[B*e9*e12*e17*e16*e8*e15*(-1+18*xi+27*xi2),B*e9*e12*e17*e35*e7*e10*e13,B*e49*e16*e8*e15*e7*e10*e13];
            N(57,:)=[B*e9*e12*e17*e16*e11*e15*e26,B*e9*e12*e17*e40*e14*e7*e10,B*e49*e16*e11*e15*e14*e7*e10];
            N(58,:)=[-C*e9*e12*e17*e16*e11*e15*(-2*xi-3+9*xi2),-C*e9*e12*e17*e40*e14*e7*e13,-C*e50*e16*e11*e15*e14*e7*e13];
            N(59,:)=[C*e9*e12*e17*e16*e11*e15*e24,C*e9*e12*e17*e40*e14*e10*e13,C*e50*e16*e11*e15*e14*e10*e13];
            N(60,:)=[-B*e9*e12*e17*e16*e11*e15*(18*xi+27*xi2-1),-B*e9*e12*e17*e41*e7*e10*e13,-B*e49*e16*e11*e15*e7*e10*e13];
            N(61,:)=[-A*e9*e12*e17*e8*e11*e15*(-18*xi-1+27*xi2),-A*e9*e12*e17*(27*et2+18*et-1)*e14*e7*e10,-A*(18*ze-1+27*ze2)*e8*e11*e15*e14*e7*e10];
            N(62,:)=[B*e9*e12*e17*e8*e11*e15*(-2*xi-3+9*xi2),B*e9*e12*e17*(27*et2+18*et-1)*e14*e7*e13,B*e49*e8*e11*e15*e14*e7*e13];
            N(63,:)=[-B*e9*e12*e17*e8*e11*e15*e44,-B*e9*e12*e17*(18*et+27*et2-1)*e14*e10*e13,-B*e49*e8*e11*e15*e14*e10*e13];
            N(64,:)=[A*e9*e12*e17*e8*e11*e15*e37,A*e9*e12*e17*e28*e7*e10*e13,A*e49*e8*e11*e15*e7*e10*e13];
            
        end
        
        
        function draw (self, gv, context)
            % Produce a graphic representation.
            %
            % function draw (self, gv, context)
            %
            %
            % Input arguments
            % self = self
            % gv = graphic viewer
            % context = struct
            % with mandatory fields
            %    x=reference geometry nodal_field
            %    u=displacement nodal_field
            % with optional fields
            %    facecolor = color of the faces, solid
            %    colorfield =nodal_field with vertex colors
            %       only one of the facecolor and colorfield may be supplied
            %    shrink = shrink factor
            %
            % these are all defined in clockwise order, which is what Matlab expects
            context.faces=[1,5,6,2;2,6,7,3;3,7,8,4;5,9,10,6;6,10,11,7;7,11,12,8;9,13,14,10;10,14,15,11;11,15,16,12;
                49,50,54,53;50,51,55,54;51,52,56,55;53,54,58,57;54,55,59,58;55,56,60,59;57,58,62,61;58,59,63,62;59,60,64,63;
                4,8,24,20;8,12,28,24;12,16,32,28;20,24,40,36;24,28,44,40;28,32,48,44;36,40,56,52;40,44,60,56;44,48,64,60;
                16,15,31,32;15,14,30,31;14,13,29,30;32,31,47,48;31,30,46,47;30,29,45,46;48,47,63,64;47,46,62,63;46,45,61,62;
                1,2,18,17;2,3,19,18;3,4,20,19;17,18,34,33;18,19,35,34;19,20,36,35;33,34,50,49;34,35,51,50;35,36,52,51;
                13,9,25,29;9,5,21,25;5,1,17,21;29,25,41,45;25,21,37,41;21,17,33,37;45,41,57,61;41,37,53,57;37,33,49,53];
            context.edges= [1,2,3,4,8,12,16,15,14,13,9,5,1;
                49,50,51,52,56,60,64,63,62,61,57,53,49;
                4,8,12,16,32,48,64,60,56,52,36,20,4;
                1,17,33,49,53,57,61,45,29,13,9,5,1];
            node_pc(1,:)=[-1,-1,-1];
            node_pc(2,:)=[-0.333333333333333,-1,-1];
            node_pc(3,:)=[0.333333333333333,-1,-1];
            node_pc(4,:)=[1,-1,-1];
            node_pc(5,:)=[-1,-0.333333333333333,-1];
            node_pc(6,:)=[-0.333333333333333,-0.333333333333333,-1];
            node_pc(7,:)=[0.333333333333333,-0.333333333333333,-1];
            node_pc(8,:)=[1,-0.333333333333333,-1];
            node_pc(9,:)=[-1,0.333333333333333,-1];
            node_pc(10,:)=[-0.333333333333333,0.333333333333333,-1];
            node_pc(11,:)=[0.333333333333333,0.333333333333333,-1];
            node_pc(12,:)=[1,0.333333333333333,-1];
            node_pc(13,:)=[-1,1,-1];
            node_pc(14,:)=[-0.333333333333333,1,-1];
            node_pc(15,:)=[0.333333333333333,1,-1];
            node_pc(16,:)=[1,1,-1];
            node_pc(17,:)=[-1,-1,-0.333333333333333];
            node_pc(18,:)=[-0.333333333333333,-1,-0.333333333333333];
            node_pc(19,:)=[0.333333333333333,-1,-0.333333333333333];
            node_pc(20,:)=[1,-1,-0.333333333333333];
            node_pc(21,:)=[-1,-0.333333333333333,-0.333333333333333];
            node_pc(22,:)=[-0.333333333333333,-0.333333333333333,-0.333333333333333];
            node_pc(23,:)=[0.333333333333333,-0.333333333333333,-0.333333333333333];
            node_pc(24,:)=[1,-0.333333333333333,-0.333333333333333];
            node_pc(25,:)=[-1,0.333333333333333,-0.333333333333333];
            node_pc(26,:)=[-0.333333333333333,0.333333333333333,-0.333333333333333];
            node_pc(27,:)=[0.333333333333333,0.333333333333333,-0.333333333333333];
            node_pc(28,:)=[1,0.333333333333333,-0.333333333333333];
            node_pc(29,:)=[-1,1,-0.333333333333333];
            node_pc(30,:)=[-0.333333333333333,1,-0.333333333333333];
            node_pc(31,:)=[0.333333333333333,1,-0.333333333333333];
            node_pc(32,:)=[1,1,-0.333333333333333];
            node_pc(33,:)=[-1,-1,0.333333333333333];
            node_pc(34,:)=[-0.333333333333333,-1,0.333333333333333];
            node_pc(35,:)=[0.333333333333333,-1,0.333333333333333];
            node_pc(36,:)=[1,-1,0.333333333333333];
            node_pc(37,:)=[-1,-0.333333333333333,0.333333333333333];
            node_pc(38,:)=[-0.333333333333333,-0.333333333333333,0.333333333333333];
            node_pc(39,:)=[0.333333333333333,-0.333333333333333,0.333333333333333];
            node_pc(40,:)=[1,-0.333333333333333,0.333333333333333];
            node_pc(41,:)=[-1,0.333333333333333,0.333333333333333];
            node_pc(42,:)=[-0.333333333333333,0.333333333333333,0.333333333333333];
            node_pc(43,:)=[0.333333333333333,0.333333333333333,0.333333333333333];
            node_pc(44,:)=[1,0.333333333333333,0.333333333333333];
            node_pc(45,:)=[-1,1,0.333333333333333];
            node_pc(46,:)=[-0.333333333333333,1,0.333333333333333];
            node_pc(47,:)=[0.333333333333333,1,0.333333333333333];
            node_pc(48,:)=[1,1,0.333333333333333];
            node_pc(49,:)=[-1,-1,1];
            node_pc(50,:)=[-0.333333333333333,-1,1];
            node_pc(51,:)=[0.333333333333333,-1,1];
            node_pc(52,:)=[1,-1,1];
            node_pc(53,:)=[-1,-0.333333333333333,1];
            node_pc(54,:)=[-0.333333333333333,-0.333333333333333,1];
            node_pc(55,:)=[0.333333333333333,-0.333333333333333,1];
            node_pc(56,:)=[1,-0.333333333333333,1];
            node_pc(57,:)=[-1,0.333333333333333,1];
            node_pc(58,:)=[-0.333333333333333,0.333333333333333,1];
            node_pc(59,:)=[0.333333333333333,0.333333333333333,1];
            node_pc(60,:)=[1,0.333333333333333,1];
            node_pc(61,:)=[-1,1,1];
            node_pc(62,:)=[-0.333333333333333,1,1];
            node_pc(63,:)=[0.333333333333333,1,1];
            node_pc(64,:)=[1,1,1];
            if isfield(context,'shrink')
                c=ones(size(node_pc,1),1)*(sum(node_pc)/size(node_pc,1));
                pc=c+context.shrink*(node_pc-c);
                for j=1:size(pc,1)
                    context.shrunk_pc_N{j}=bfun(self,pc(j,:));
                end
            end
            draw@fe_set_3_manifold (self, gv, context);
        end
        
        function   draw_isosurface(self, gv, context)
            % Draw isosurface within a three-dimensional manifold finite element.
            %
            % function   draw_isosurface(self, gv, context)
            %
            % Input arguments
            % self = descendent of fe_set_3_manifold
            % gv = graphic viewer
            % context = struct
            % with mandatory fields
            %    x=reference geometry nodal_field
            %    u=displacement field
            % with optional fields
            %    scalarfield =field with scalar vertex data,
            %    isovalue = value of the isosurface
            %    color = what color should the isosurface be?  Default: [1,1,1]/2;
            %    edgealpha= transparency value
            %    linewidth=line width (default 2)
            conns = self.conn; % connectivity
            [f,v]= isosurf_faces(conns, context.geom, context.u, context.scalarfield, context.isovalue, 5);
            
            if (isfield(context,'color'))
                color =context.color;
            else
                color =[1,1,1]/2;
            end
            
            facealpha = 1;
            if (isfield(context,'facealpha'))
                facealpha = context.facealpha;
            end
            patch('Faces',f,'Vertices',v,'FaceColor',color,'EdgeColor','none','FaceAlpha', facealpha);
        end
        
        
        
        function val =boundary_conn(self)
            % Get boundary connectivity.
            conn =self.conn;
            val=[conn(:,[1,5,9,13,2,6,10,14,3,7,11,15,4,8,12,16]);
                conn(:,[4,8,12,16,20,24,28,32,36,40,44,48,52,56,60,64]);
                conn(:,[16,15,14,13,32,31,30,29,48,47,46,45,64,63,62,61]);
                conn(:,[13,9,5,1,29,25,21,17,45,41,37,33,61,57,53,49]);
                conn(:,[64,63,62,61,60,59,58,57,56,55,54,53,52,51,50,49]);
                conn(:,[1,2,3,4,17,18,19,20,33,34,35,36,49,50,51,52])];
        end
        
        function val =boundary_fe(self)
            % Get the constructor of the class of the  boundary finite element.
            val = @fe_set_Q16;
        end
        
    end
    
end


