% Compute load vector for nonzero essential boundary conditions.
%
% function F = nz_ebc_loads(self, assembler, geom, u)
% 
%
% Return the assembled system vector F.
%    Arguments
%     assembler =  descendent of sysvec_assembler
%     geom=geometry field
%     u=displacement field
%
function F = nz_ebc_loads(self, assembler, geom, u)
    fes = self.fes;% grab the finite elements to work on
    %     Evaluate the nodal basis function gradients
    bfun_gradients = nodal_bfun_gradients (self, geom);
    % Material orientation
    Rm_constant = is_material_orientation_constant(self);% if not constant, need to compute  at each point
    if (~Rm_constant)
        Rmh = self.Rm;% handle to a function  to evaluate Rm
    else
        Rm = self.Rm;
    end    
    % Material
    mat = self.material;
    D_constant = are_tangent_moduli_constant (mat);
    if (D_constant)
        D = tangent_moduli(mat,[]);
        [Princd,W]=eig(D);
        [ignore,ix]=sort(diag(W));
        Wmean=(W(ix(1),ix(1))+W(ix(2),ix(2)))/2;
    end
    % Retrieve data for efficiency
    conns = fes.conn; % connectivity
    labels = fes.label; % connectivity
    xs =geom.values;
    % Prepare assembler
    start_assembly(assembler, u.nfreedofs);
    for nix=1:length(bfun_gradients)
        pu = reshape(u, gather_fixed_values(u, bfun_gradients{nix}.patchconn));
        if norm (pu) ~= 0
            if bfun_gradients{nix}.Vpatch>0 % This node has a element patch in this block
                np=length(bfun_gradients{nix}.patchconn);
                c = xs(nix, :); % coordinates of central node
                D = tangent_moduli (mat, struct ('xyz',c));% material tangent
                if (~Rm_constant)% do I need to evaluate the local material orientation?
                    if (~isempty(labels )),  Rm =Rmh(c,J,labels(i));
                    else,                    Rm =Rmh(c,J,[]);                end
                end
                Bnodal= self.hBlmat (self,bfun_gradients{nix}.Nspd,c,Rm);
                if (~D_constant)
                    D = tangent_moduli (mat, struct ('xyz',c));% material tangent
                    [Princd,W]=eig(D);
                    [ignore,ix]=sort(diag(W));
                    Wmean=(W(ix(1),ix(1))+W(ix(2),ix(2)))/2;
                end
                Ke = Bnodal'*D*Bnodal * bfun_gradients{nix}.Vpatch;
                if (self.stabfact>0)
                    lx = xs(bfun_gradients{nix}.patchconn,:); % coordinates of nodes
                    Phi =self.hPhi(self,u.dim,np,lx,c);
                    lastwarn('','')
                    A1 =Phi* ((Phi'*Phi)\Phi');
                    [msgstr, msgid] = lastwarn;
                    if (isempty(msgid))
                        Ke=Ke + self.stabfact*Wmean*(eye(size(A1))-A1);
                    else
                        disp(msgid)
                    end
                end
                eqnums =reshape(u,gather_dofnums(u,bfun_gradients{nix}.patchconn));
                Fe =  -Ke*pu;
                assemble(assembler, Fe, eqnums);
            end
        end    
    end
    F = make_vector (assembler);
end
