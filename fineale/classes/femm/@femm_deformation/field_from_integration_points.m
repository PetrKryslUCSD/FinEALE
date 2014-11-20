function fld = field_from_integration_points(self, ...
        geom, u, dT, output, component, options)
% Create a field from quantities at integration points.
%
% function fld = field_from_integration_points(self, ...
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
    fes = self.fes;
    % Make the inverse map from finite element nodes to gcells
    fen2fe_map =node_to_element_map(self);
    gmap=fen2fe_map.map;
    % Container of intermediate results
    sum_inv_dist =0*(1:length(gmap))';
    sum_quant_inv_dist =zeros(length(gmap),length(component));
    fld =nodal_field(struct ('name',['fld'], 'dim', length(component), 'data',zeros(u.nfens,length(component))));
    nvals = fld.values;
    % This is an inverse-distance interpolation inspector.
    x= []; conn = [];
    function idat=idi(idat, out, xyz, u, pc)
        d=x-ones(length(conn),1)*xyz;
        ld=reshape(d',prod(size(d)),1);
        d =sum(reshape(ld.*ld,size(d,2), size(d,1)));
        zi = find(d==0);
        if (~isempty(zi))
            nzi = find(d>0);
            d(zi) =min(d(nzi))/1e9;
        end
        invd =1./d;
        quant=reshape(out(component),1,length(component));
        sum_quant_inv_dist(conn,:)=sum_quant_inv_dist(conn,:) + invd'*quant;
        sum_inv_dist(conn)=sum_inv_dist(conn)+invd';
    end
    conns = fes.conn; % connectivity
    if (exist( 'options','var' ) && (~isempty(options)))
        context = options;
    end
    context.output = output;
    % Loop over cells to interpolate to nodes
    idat=0;
    for i=1:count(fes)
        conn = conns(i,:);
        x=gather_values(geom,conn);
        idat = inspect_integration_points(self, geom, u, dT, ...
            i, context, @idi, idat);
    end
    % compute the data array
    nzi = find (sum_inv_dist);
    for j=1:length(component)
        nvals(nzi,j)=sum_quant_inv_dist(nzi,j)./sum_inv_dist(nzi);
    end
    % Make the field
    fld = scatter(fld, (1: fld.nfens),nvals);
    return;
end
