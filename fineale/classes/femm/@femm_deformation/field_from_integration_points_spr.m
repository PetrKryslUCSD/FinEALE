function fld = field_from_integration_points_spr(self, ...
        geom, u, dT, output, component,  options)
% Create a field from quantities at integration points using 
% the Superconvergent Patch Recovery (SPR).
%
% function fld = field_from_integration_points_spr(self, ...
%         geom, u, dT, output, component, options)
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
    fes = self.fes;
    % Make the inverse map from finite element nodes to gcells
    fen2fe_map =node_to_element_map(self);
    gmap=fen2fe_map.map;
    % Container of intermediate results
    sum_inv_dist =0*(1:length(gmap))';
    sum_quant_inv_dist =zeros(length(gmap),length(component));
    fld =nodal_field(struct ('name',['fld'], 'dim', length(component), 'data',zeros(u.nfens,length(component))));
    nvals = fld.values;
    % connectivity
    conns = fes.conn; % connectivity
    xs =geom.values;;
    % This is an super convergent patch recovery interpolation inspector.
    x= []; conn = []; spr_data= {}; spr_data_idx=0;
    function idat=idi(idat, out, xyz, u, pc)
        spr_data_idx=spr_data_idx+1;
        spr_data{spr_data_idx}.xyz=xyz;
        spr_data{spr_data_idx}.out=reshape(out(component),1,length(component));;
    end
    if exist( 'options','var' ) && (~isempty(options))
        context = options;
    end
    context.output = output;
    % Loop over cells to interpolate to nodes
    idat=0;
    for i=1:length(gmap)
        pgl=gmap{i};
        spr_data_idx=0;
        xc=xs(i,:);
        for k=1:length(pgl)
            conn = conns(pgl(k),:);
            x=xs(conn,:);
            idat = inspect_integration_points(self, geom, u, dT, ...
                pgl(k), context, @idi, idat);
        end
        nvals(i,:)=spr(spr_data,spr_data_idx,xc);
    end
    % Make the field
    fld = scatter(fld, (1: fld.nfens),nvals);
    return;
end

function out=spr(spr_data,n,xc)
    if (n==0)
        out=0;
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