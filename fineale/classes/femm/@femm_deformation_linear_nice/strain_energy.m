function SE = strain_energy (self, geom, u)
% Compute the strain energy by computing the strain energy of the individual FEs.
%
% function SE = strain_energy (self, geom, u)
%
% Return a double value for the sum of the element energies.
%     geom=geometry field
%     u=displacement field
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
end
% Retrieve data for efficiency
conns = fes.conn; % connectivity
labels = fes.label; % connectivity
xs =geom.values;
% Prepare assembler
Kedim =3*u.dim*fes.nfens;% rough estimate
% Initialize the accumulator.
SE = 0;
for nix=1:length(bfun_gradients)
    if bfun_gradients{nix}.Vpatch~=0 % This node has a element patch in this block
        np=length(bfun_gradients{nix}.patchconn);
        c = xs(nix, :); % coordinates of central node
        D = tangent_moduli (mat, struct ('xyz',c));% material tangent
        if (~Rm_constant)% do I need to evaluate the local material orientation?
            if (~isempty(labels )),  Rm =Rmh(c,J,labels(i));
            else,                    Rm =Rmh(c,J,[]);                end
        end
        Bnodal= self.hBlmat (self,bfun_gradients{nix}.Nspd,c,Rm);
        if (~D_constant)
            D = tangent_moduli(mat,struct('xyz',c));
        end
        Ke = Bnodal'*D*Bnodal * bfun_gradients{nix}.Vpatch;
        U =reshape(u,gather_values(u,bfun_gradients{nix}.patchconn));
        SE = SE + U'*Ke*U/2;
    end
end
end


