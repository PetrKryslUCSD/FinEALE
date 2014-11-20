function model_data =heat_diffusion_adaptive_2D_steady_state(region_definition,make_model_data,adaptivity_options)
% Adaptive steady-state heat conduction solver.
%
% function model_data =heat_diffusion_adaptive_2D_steady_state(model_data)
%
%
% Arguments: 
%
% model_data = struct  with fields as discussed for  the algorithm
% heat_diffusion_steady_state. Additionally, the following fields 
% need to be supplied: 
%
% model_data.region_definition= cell array to describe the region geometry
% % region_definition={'REGEXT 0 0 0.6 1.0',...
% %     ['curve 1 line 0 0 0.6 0'],...
% %     ['curve 2 line 0.6 0 0.6 0.2'],...
% %     ['curve 3 line 0.6 0.2 0.6 1.0'],...
% %     ['curve 4 line 0.6 1.0 0 1.0'],...
% %     ['curve 5 line 0 1.0 0 0'],...
% %     ['subregion 1  property 1 boundary 1 2 3 4 5']};
% Caveat: Only as single-region geometry is supported at the moment.
%
% adaptivity_options.initial_mesh_size= initial mesh size,
% adaptivity_options.mesh_options= mesh options (see targe2_mesher())
 %
% Output:
% model_data= the struct on input is augmented as discussed for the
% algorithm heat_diffusion_steady_state   
% 

% Check the control parameters

nadapt = 3; % How many adaptive steps should we make?
if (isfield(adaptivity_options,'nadapt'))
    nadapt  =adaptivity_options.nadapt;;
end
observer =@(step,model_data) disp(['Adaptive step ' num2str(step)]);
if ( isfield(adaptivity_options,'observer'))
    observer  =adaptivity_options.observer;;
end

% Generate the initial mesh
[fens,fes,groups,edge_fes,edge_groups]...
    =targe2_mesher(cat(2,region_definition,{['m-ctl-point constant ' num2str(adaptivity_options.initial_mesh_size)]}), 1.0,...
    adaptivity_options.mesh_options);

step=0;
while true
    step=step+1;
    
    model_data= make_model_data (fens,fes,groups,edge_fes,edge_groups);
    
    model_data =heat_diffusion_steady_state(model_data);
    
    
    if ~isempty(observer)% report the progress
            observer (step,model_data);
    end
        
    % Have we gone through all the adaptivity attempts?   If so, we are done.
    if  step == nadapt,
        break;
    end
   
    nodal_fluxes = field_from_integration_points(model_data.region{1}.femm, model_data.geom, model_data.temp,  'flux', 1:2);
    elerrs = flux_L2_error (model_data.region{1}.femm, model_data.geom, model_data.temp, nodal_fluxes);
    total_err=sqrt(sum(elerrs.^2));
    targeterr=sqrt(total_err^2/adaptivity_options.targetnel) ;
    
    [hcurs, hests] =T3_mesh_sizes(fes.conn,model_data.geom.values,targeterr,elerrs,adaptivity_options.convergence_rate);
    
    [fens,fes,groups,edge_fes,edge_groups] ...
        = targe2_mesher_adapt(region_definition,fes.conn,model_data.geom.values,hests,1.0,adaptivity_options.mesh_options);
end

end
