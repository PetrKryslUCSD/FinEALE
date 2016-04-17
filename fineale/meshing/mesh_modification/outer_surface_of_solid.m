function osfes=outer_surface_of_solid(fens,bdry_fes,parametric_centroid)
% Select the geometric cells that form the outer boundary surface.
%
% function osfes=outer_surface_of_solid(fens,bdry_fes,parametric_centroid)
%
%     fens=finite element nodes,
%     bdry_fes=boundary geometric cells,
%     parametric_centroid=where is the centroid for the geometric cells
%          in parametric coordinates (for instance for three node triangle,
%          parametric_centroid =[1/3,1/3])

xyz=fens.xyz;
% This is the centroid of the cloud of nodes
centroid =mean(xyz);
% for axially symmetric geometry placed the centroid on the axis
if (bdry_fes.axisymm)
    centroid(1)=0;
end

% How do we compute the normal to the  boundary face?
if (size(xyz,2)==3)
    compute_normal=@(J)skewmat(J(:,1))*J(:,2);
elseif (size(xyz,2)==2)
    compute_normal=@(J)[J(2);-J(1)];
end

angtol=0.001;

% If we have 2-D axial symmetry: eliminate  boundary cells that are on the
% axis of symmetry
onaxisl  = [];
if (bdry_fes.axisymm)
    onaxisl = fe_select(fens,bdry_fes,struct ('box',[0,0,-inf,inf],'inflate',norm(max(xyz)-min(xyz))/10000));
    % Prune the  list of cells to include by flooding below
    bdry_fes=subset(bdry_fes,setdiff(1:count(bdry_fes),onaxisl));;
end

% Now we will try to find  one  boundary cell  so that  all the nodes lie
% in the half space  behind it (that is in the other half space from the
% one into which  its  normal is pointing), and the half space also
% includes the centroid.
conns = bdry_fes.conn; % connectivity
start =[];
for j=1:count(bdry_fes)
    x=xyz(conns(j,:),:);
    N= bfun(bdry_fes,parametric_centroid);
    c=(x'*N)';
    gradNpar = bfundpar(bdry_fes,parametric_centroid);
    J = x' * gradNpar;
    n=compute_normal(J);
    n=n/norm(n);
    if (all_behind(xyz,centroid,c,n,conns(j,:),angtol))
        start=j;
        break;
    end
end 

% Could not find a single  surface cell  so that all nodes are behind it:
% failure.
if (isempty(start))
    osfes=[]; return
end

startfen =conns(start,1);

% Select all cells connected together
osfesl = fe_select(fens,bdry_fes,...
    struct ('flood', true, 'startfen', startfen));
osfes=subset(bdry_fes,osfesl);
return;

%
%     gv=drawmesh({fens,osfes},'fes','facecolor','r','facealpha',0.3)
%     gv=drawmesh({fens,osfes},'gv',gv,'fes','facecolor','y','facealpha',0.3)
%     set_graphics_defaults
%     labels('x','y','z')

function AllBehind=  all_behind(xyz,centroid,c,n,conn,tol)
    AllBehind = true;
    %     Check the centroid first
    v=centroid-c;
    vn=norm(v);
    if (vn>0)
        d=dot(n,v)/vn;
        if (d>tol)
            AllBehind = false; return
        end
    end
    % Now check all the vertices
    % except for the nodes which are  connected by the tested cell
    list=setdiff(1:size(xyz,1),conn);
    for kl=1:length(list)
        k=list(kl);
        v=xyz(k,:)-c;
        vn=norm(v);
        if (vn>0)
            d=dot(n,v)/vn;
            if (d>tol)
                AllBehind = false; return
            end
        end
    end 
end
 
end

