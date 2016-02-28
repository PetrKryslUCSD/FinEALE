% Compute the stiffness matrix by computing and assembling the
% matrices of the individual FEs.
%
% function K = stiffness (self, assembler, geom, u)
%
% Return a stiffness matrix.
%     assembler = descendent of the sysmat_assembler class
%     geom=geometry field
%     u=displacement field
function K = stiffness_stab_filter (self, assembler, geom, u)
% Compute the nodal basis function gradients.
% Return the cell array of structures with attributes
%		 bfun_gradients{nix}.Nspd= basis function gradient matrix
%        bfun_gradients{nix}.Vpatch= nodal patch volume
%        bfun_gradients{nix}.patchconn= nodal patch connectivity
%
fes = self.fes;% grab the finite elements to work on
%     Evaluate the nodal basis function gradients
bfun_gradients = nodal_bfun_gradients (self, geom);
% Material orientation
Rm_identity = is_material_orientation_identity(self);% if identity, work not needed
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
Kedim =3*u.dim*fes.nfens;% rough estimate
start_assembly(assembler, Kedim, Kedim, size(conns,1), u.nfreedofs, u.nfreedofs);
for nix=1:length(bfun_gradients)
    if bfun_gradients{nix}.Vpatch~=0 % This node has a element patch in this block
        np=length(bfun_gradients{nix}.patchconn);
        c = xs(nix, :); % coordinates of central node
        if (~D_constant)
            D = tangent_moduli(mat,struct('xyz',c));
            [Princd,W]=eig(D);
            [ignore,ix]=sort(diag(W));
            Wmean=(W(ix(1),ix(1))+W(ix(2),ix(2)))/2;
        end
        lx = xs(bfun_gradients{nix}.patchconn,:); % coordinates of nodes
        Phi =self.hPhi(self,u.dim,np,lx,c);
        lastwarn('','')
        A1 =Phi* ((Phi'*Phi)\Phi');
        [msgstr, msgid] = lastwarn;
        if (isempty(msgid))
            Kstab=self.stabfact*Wmean*(eye(size(A1))-A1);
            eqnums =reshape(u,gather_dofnums(u,bfun_gradients{nix}.patchconn));
            assemble_symmetric(assembler, Kstab, eqnums);
        else
            disp(msgid)
        end
        
    end
end
K = make_matrix (assembler);
end
