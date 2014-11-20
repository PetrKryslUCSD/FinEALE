% Create a field from quantities at integration points using 
% the Superconvergent Patch Recovery (SPR).
%
% function fld = field_from_integration_points_spr(self, geom, T, output, component, options)
%
% Input arguments
%     geom     - reference geometry field
%     u        - displacement field
%     dT       - temperature difference field
%     output   - this is what you would assign to the 'output' attribute of
%                the context argument of the material update() method.
%     component- component of the 'output' array: see the material update()
%                method.
%     options  - struct with fields recognized by update() method of the 
%                material; optional argument
% Output argument
%     fld - the new field that can be used to map values the colors
%
% Create a field from quantities at integration points using 
% the Superconvergent Patch Recovery (SPR) as formulated by Zienkiewicz and 
% Zhu. Linear polynomial Least-squares fit is used for all elements.
%
function fld = field_from_integration_points_spr(self, geom, T, output, component, options)
    gcells = get(self.femmlock,'gcells');
    % Make the inverse map from finite element nodes to gcells
    f2gmap=fenode_to_gcell_map (struct ('gcells',gcells));
    gmap=get (f2gmap,'gcell_map');
    % Container of intermediate results
    sum_inv_dist =0*(1:length(gmap))';
    sum_quant_inv_dist =zeros(length(gmap),length(component));
    fld =field(struct ('name',['fld'], 'dim', length(component), 'data',zeros(get(T,'nfens'),length(component))));
    nvals = gather (fld, (1: get(fld,'nfens')),'values','noreshape');
    % This is an inverse-distance interpolation inspector.
    x= []; conn = []; spr_data= {}; spr_data_idx=0;
    function idat=idi(idat, out, xyz, u, pc)
        spr_data_idx=spr_data_idx+1;
        spr_data{spr_data_idx}.xyz=xyz;
        spr_data{spr_data_idx}.out=reshape(out(component),1,length(component));;
    end
     conns = get(gcells, 'conn'); % connectivity
   % Loop over cells to interpolate to nodes
    idat=0;
    for i=1:length(gmap)
        pgl=gmap{i};
        spr_data_idx=0;
        xc=gather(geom,i,'values','noreshape');
        for k=1:length(pgl)
            conn = conns(pgl(k),:);
            x=gather(geom,conn,'values','noreshape');
            idat = inspect_integration_points(self, geom, T, ...
                pgl(k), struct ('output',output), @idi, idat);
        end
        nvals(i,:)=spr(spr_data,spr_data_idx,xc);
    end
    % Make the field
    fld = scatter(fld, (1: get(fld,'nfens')),nvals);
    return;
end

function out=spr(spr_data,n,xc)
    if (n==0)
        error('No data');
    elseif (n==1)
        out=spr_data{1}.out;
    else
        dim=length(spr_data{1}.xyz);
        if (dim==1)
            if (n>=2)
                na=2;
                A=zeros(na,na);
                nc=length(spr_data{1}.out);
                b=zeros(na,nc);
                for k=1:n
                    xk=spr_data{k}.xyz-xc;
                    pk= [1,xk]; 
                    A=A+pk'*pk;
                    b=b+pk'*spr_data{k}.out;
                end 
                a=A\b;
                xk=xc-xc;
                p= [1,xk];
                out=p*a;
            else
                out=spr_data{1}.out/n;
                for k=2:n
                    out=out+spr_data{k}.out/n;
                end 
            end
        elseif (dim==2)
% if (n>=6)
%                 na=6;
%                 A=zeros(na,na);
%                 nc=length(spr_data{1}.out);
%                 b=zeros(na,nc);
%                 for k=1:n
%                     xk=spr_data{k}.xyz-xc;
%                     pk= [1,xk,xk.^2,prod(xk)]; 
%                     A=A+pk'*pk;
%                     b=b+pk'*spr_data{k}.out;
%                 end 
%                 a=A\b;
%                 xk=xc-xc;
%                 p= [1,xk,xk.^2,prod(xk)];
%                 out=p*a;
%             else
            if (n>=3)
                na=3;
                A=zeros(na,na);
                nc=length(spr_data{1}.out);
                b=zeros(na,nc);
                for k=1:n
                    xk=spr_data{k}.xyz-xc;
                    pk= [1,xk]; 
                    A=A+pk'*pk;
                    b=b+pk'*spr_data{k}.out;
                end 
                a=A\b;
                xk=xc-xc;
                p= [1,xk];
                out=p*a;
            else
                out=spr_data{1}.out/n;
                for k=2:n
                    out=out+spr_data{k}.out/n;
                end 
            end
        elseif (dim==3)
            if (n>=4)
                na=4;
                A=zeros(na,na);
                nc=length(spr_data{1}.out);
                b=zeros(na,nc);
                for k=1:n
                    xk=spr_data{k}.xyz-xc;
                    pk= [1,xk]; 
                    A=A+pk'*pk;
                    b=b+pk'*spr_data{k}.out;
                end 
                a=A\b;
                xk=xc-xc;
                p= [1,xk];
                out=p*a;
            else
                out=spr_data{1}.out/n;
                for k=2:n
                    out=out+spr_data{k}.out/n;
                end 
            end
        else
            error('Unknown space and dimension');
        end
    end
end