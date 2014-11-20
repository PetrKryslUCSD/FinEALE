function H = surface_normal_spring_stiffness(self, assembler, geom, u)
% Compute the stiffness matrix of surface normal spring.
%
% function H = surface_normal_spring_stiffness(self, assembler, geom, u)
%
% Arguments
%           self  = heat diffusion model  
%           assembler = descendent of the sysmat_assembler class
%           geom=geometry field
%           u=displacement field
%
% Returns H as a matrix.
%
% Rationale: consider continuously distributed springs between the surface of 
% the solid body and the 'ground', in the direction normal to the surface. 
% If the spring coefficient becomes large, we have an approximate
% method of enforcing the normal displacement to the surface.
 error('Base class does not implement this behaviour!');
end