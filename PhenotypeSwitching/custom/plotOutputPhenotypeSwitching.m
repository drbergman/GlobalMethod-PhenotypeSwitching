function figs = plotOutputPhenotypeSwitching(out,input,nsamps,summary_flags)

% for i = 1:length(input)
%     input(i) = input(i).parameter_setup_fn(input(i));
% end
nsubcohorts = numel(out)/nsamps;
TME_size = input(1).TME_size;
Tend = input(1).Tend;
save_dt = input(1).save_dt;
N = prod(TME_size);
nt = round(Tend/(save_dt)+1);
t = linspace(0,Tend,nt);
[NumCells,InternSubstrates,AllConcentrations] = normalizeOutputPhenotypeSwitching(out,input,nsubcohorts,nsamps);

%% patch plot
if summary_flags.show_patch_plots || summary_flags.make_svgs
    for si = 1:min(1,nsamps)
        for i = 1:nsubcohorts
%             if input(i).substrates(1).method=="global"
%                 Regions = setupRegions(input(i).substrates(1),input(i).TME_size,input(i).bv_location);
%             else
%                 Regions = reshape(1:N,TME_size);
%             end
            if input(i).substrates(1).method=="local" && isempty(input(i).substrates(1).solver.regions)
                input(i).substrates(1).solver.regions = reshape(1:prod(input(i).TME_size),input(i).TME_size);
            end
            actually_make_svgs = summary_flags.make_svgs && si==1;
            if input(i).ntypes == 3
                cell_labels = ["Proliferating","Quiescent","Necrotic"];
            else
                cell_labels = ["Proliferating","Quiescent"];
            end
            if length(TME_size)==3
                scatterPlot(input(i).TME_size,[1,1,1],out(i,si),1,input(i).TME_size,actually_make_svgs,1,cell_labels,100,input(i).substrates(1).method,"octant")
            else
                patchPlot4(TME_size,[1,1],out(i,si),1,TME_size,actually_make_svgs,1,cell_labels,100,input(i).substrates(1).method,1,1:length(input(i).substrates),input(i).substrates)
            end
        end
    end

end

%% plotting
disp_names = strings(nsubcohorts,1);
for i = 1:nsubcohorts
    if input(i).substrates(1).method=="local"
        disp_names(i) = 'local';
    elseif input(i).substrates(1).pars.regions_type == "concentric_circles"
        disp_names(i) = sprintf('m = %d',input(i).substrates(1).pars.num_concentric_circle_regions);
    elseif input(i).substrates(1).pars.regions_type == "rectangular_grid"
        if length(input(i).substrates(1).pars.rectangle_compartment_size) == 2
            disp_names(i) = sprintf('m = (%d,%d)',input(i).substrates(1).pars.rectangle_compartment_size(1),input(i).substrates(1).pars.rectangle_compartment_size(2));
        else
            disp_names(i) = sprintf('m = (%d,%d,%d)',input(i).substrates(1).pars.rectangle_compartment_size(1),input(i).substrates(1).pars.rectangle_compartment_size(2),input(i).substrates(1).pars.rectangle_compartment_size(3));
        end
    else
        disp_names(i) = sprintf('#%d',i);
    end
end

figs(1).handle=figure; 
figs(1).name = "summary";
nr = 2;
nc = 2*input(1).ntypes;
ax_vals = gobjects(1,nc);
ax_rels = gobjects(1,nc);
for type_ind = 1:nc
    ax_vals(type_ind) = subplot(nr,nc,type_ind); hold on;
    ax_rels(type_ind) = subplot(nr,nc,r2c(nr,nc,[2,type_ind])); hold on;
end
set(ax_vals(1:input(1).ntypes),'YLim',[0 N])
cmap = cubehelix(nsubcohorts,1.51,0.93,2.11,0.78,[0.12,0.54],[0.22,0.68]);
vbar_abm = mean(NumCells(:,1,:,:),3);
pa_vals = gobjects(nsubcohorts,nc);
pa_rels = gobjects(nsubcohorts,nc);
ls_vals = gobjects(nsubcohorts,nc);
ls_rels = gobjects(nsubcohorts,nc);
for i = 1:nsubcohorts

    for type_ind = 1:input(1).ntypes

        v = NumCells(:,i,:,type_ind);
        vbar = mean(v,3);
        dv = v - vbar_abm(:,:,:,type_ind);
        rel_dv = 100*dv./vbar_abm(:,:,:,type_ind);
        sd = std(dv,[],3);
        rel_sd = std(rel_dv,[],3);
        rel_dv_bar = mean(rel_dv,3);
        pa_vals(i,type_ind) = patch(ax_vals(type_ind),[t,flip(t)],reshape([vbar,flip(vbar)]+[-sd,flip(sd)],1,[]),cmap(i,:),...
            'FaceAlpha',0.2,'DisplayName',disp_names(i),'EdgeColor','none');
        pa_rels(i,type_ind) = patch(ax_rels(type_ind),[t,flip(t)],reshape([rel_dv_bar,flip(rel_dv_bar)]+[-rel_sd,flip(rel_sd)],1,[]),cmap(i,:),...
            'FaceAlpha',0.2,'DisplayName',disp_names(i),'EdgeColor','none');

        ls_vals(i,type_ind) = plot(ax_vals(type_ind),t,vbar,'Color',cmap(i,:),'LineWidth',2,'DisplayName',disp_names(i));
        ls_rels(i,type_ind) = plot(ax_rels(type_ind),t,rel_dv_bar,'Color',cmap(i,:),'LineWidth',2,'DisplayName',disp_names(i));

    end
    
    for type_ind = 1:input(1).ntypes
        isbar_abm = mean(InternSubstrates(:,1,:,type_ind),3);
        is = InternSubstrates(:,i,:,type_ind);
        isbar = mean(is,3);
        dis = is - isbar_abm;
        rel_dis = 100*dis./isbar_abm;
        sd = std(dis,[],3);
        rel_sd = std(rel_dis,[],3);
        rel_dis_bar = mean(rel_dis,3);
        pa_vals(i,input(1).ntypes+type_ind) = patch(ax_vals(input(1).ntypes+type_ind),[t,flip(t)],reshape([isbar,flip(isbar)]+[-sd,flip(sd)],1,[]),cmap(i,:),...
            'FaceAlpha',0.2,'DisplayName',disp_names(i),'EdgeColor','none');
        pa_rels(i,input(1).ntypes+type_ind) = patch(ax_rels(input(1).ntypes+type_ind),[t,flip(t)],reshape([rel_dis_bar,flip(rel_dis_bar)]+[-rel_sd,flip(rel_sd)],1,[]),cmap(i,:),...
            'FaceAlpha',0.2,'DisplayName',disp_names(i),'EdgeColor','none');

        ls_vals(i,input(1).ntypes+type_ind) = plot(ax_vals(input(1).ntypes+type_ind),t,isbar,'Color',cmap(i,:),'LineWidth',2,'DisplayName',disp_names(i));
        ls_rels(i,input(1).ntypes+type_ind) = plot(ax_rels(input(1).ntypes+type_ind),t,rel_dis_bar,'Color',cmap(i,:),'LineWidth',2,'DisplayName',disp_names(i));
    end

end
cell_type = input(1).type_names;
for ci = 1:nc
    ax_vals(ci).Children = [flip(ls_vals(:,ci));flip(pa_vals(:,ci))];
    ax_rels(ci).Children = [flip(ls_rels(:,ci));flip(pa_rels(:,ci))];
end
for type_ind = 1:input(1).ntypes
    title(ax_vals(type_ind),cell_type(type_ind));
%     title(ax_rels(type_ind),[cell_type(type_ind),sprintf(' relative to %s',disp_names(1))])
    if type_ind==1
        ylabel(ax_rels(type_ind),["% Diff of above",sprintf(' from %s',disp_names(1))])
    end
    title(ax_vals(input(1).ntypes+type_ind),['Internalized Substrates of',cell_type(type_ind)]);
%     title(ax_rels(input(1).ntypes+type_ind),['Internalized Substrates of',cell_type(type_ind),sprintf(' relative to %s',disp_names(1))])
%     ylabel(ax_rels(input(1).ntypes+type_ind),["% Diff of above",sprintf(' from %s',disp_names(1))])
end
set([ax_vals,ax_rels],'FontSize',20)


legend(ax_rels(2),pa_rels(:,2),'location','best')
set(ax_vals,'XLim',[0 Tend])
set(ax_rels,'XLim',[0 Tend])

% prepare print parameters
figs(1).fontsize = 6;
figs(1).print_position = [0 0 6 2];

%% plot substrate concentration data
figs(2).handle = figure;
figs(2).name = "substrates";
nr = 2; % actual values and relative values
nc = input(1).ntypes+3; % [in Proliferating,in quiescent,at bv,free,total]
substrate_vals_ax = gobjects(1,nc);
substrate_rels_ax = gobjects(1,nc);
for i = 1:nc
    substrate_vals_ax(i) = subplot(nr,nc,i); hold on;
    substrate_rels_ax(i) = subplot(nr,nc,r2c(nr,nc,[2,i])); hold on;
end

cmap = cubehelix(nsubcohorts,1.51,0.93,2.11,0.78,[0.12,0.54],[0.22,0.68]);
cbar_abm = mean(AllConcentrations(:,1,:,:),3);
substrate_pa_vals = gobjects(nsubcohorts,nc);
substrate_pa_rels = gobjects(nsubcohorts,nc);
substrate_ls_vals = gobjects(nsubcohorts,nc);
substrate_ls_rels = gobjects(nsubcohorts,nc);

for type_ind = 1:input(1).ntypes
    copyobj( ax_vals(input(1).ntypes+type_ind).Children, substrate_vals_ax(type_ind)   )
    copyobj( ax_rels(input(1).ntypes+type_ind).Children, substrate_rels_ax(type_ind)   )
end
for i = 1:nsubcohorts

    for j = 1:size(AllConcentrations,4)
        c = AllConcentrations(:,i,:,j);
        cbar = mean(c,3);
        dc = c - cbar_abm(:,:,:,j);
        rel_dc = 100*dc./cbar_abm(:,:,:,j);
        sd = std(dc,[],3);
        rel_sd = std(rel_dc,[],3);

        rel_dc_bar = mean(rel_dc,3);
        substrate_pa_vals(i,input(1).ntypes+j) = patch(substrate_vals_ax(input(1).ntypes+j),[t,flip(t)],reshape([cbar,flip(cbar)]+[-sd,flip(sd)],1,[]),cmap(i,:),...
            'FaceAlpha',0.2,'DisplayName',disp_names(i),'EdgeColor','none');
        substrate_pa_rels(i,input(1).ntypes+j) = patch(substrate_rels_ax(input(1).ntypes+j),[t,flip(t)],reshape([rel_dc_bar,flip(rel_dc_bar)]+[-rel_sd,flip(rel_sd)],1,[]),cmap(i,:),...
            'FaceAlpha',0.2,'DisplayName',disp_names(i),'EdgeColor','none');

        substrate_ls_vals(i,input(1).ntypes+j) = plot(substrate_vals_ax(input(1).ntypes+j),t,cbar,'Color',cmap(i,:),'LineWidth',2,'DisplayName',disp_names(i));
        substrate_ls_rels(i,input(1).ntypes+j) = plot(substrate_rels_ax(input(1).ntypes+j),t,rel_dc_bar,'Color',cmap(i,:),'LineWidth',2,'DisplayName',disp_names(i));
    end
end
for ci = (input(1).ntypes+1):nc % I think 1:2 should be already handled by the copying
    substrate_vals_ax(ci).Children = [flip(substrate_ls_vals(:,ci));flip(substrate_pa_vals(:,ci))];
    substrate_rels_ax(ci).Children = [flip(substrate_ls_rels(:,ci));flip(substrate_pa_rels(:,ci))];
end
% conc_type = {{'Avg Proliferating','Internalized Substrate'},{'Avg Quiescent','Internalized Substrate'},{'Avg [Free','Substrate] at BV'},{'Avg [Free','Substrate]'},{'Avg Total','Substrate'}};
conc_type = cell(nc,1);
for type_ind = 1:input(1).ntypes
    conc_type{type_ind} = [sprintf("Avg %s",input(1).type_names(type_ind)),"Intern Substrate"];
end
conc_type(end-2:end) = {{'Avg [Free','Substrate] at BV'},{'Avg [Free','Substrate]'},{'Avg Total','Substrate'}};
for ci = 1:nc
    title(substrate_vals_ax(ci),conc_type{ci});
%     title(substrate_rels_ax(ci),[conc_type{ci},sprintf(' relative to %s',disp_names(1))])
    if ci==1
        ylabel(substrate_rels_ax(ci),["% Diff of above",sprintf(' from %s',disp_names(1))])
    end
end
set([substrate_vals_ax,substrate_rels_ax],'FontSize',20)


legend(substrate_rels_ax(input(1).ntypes+1),substrate_pa_rels(:,input(1).ntypes+1),'location','best')
set(substrate_vals_ax,'XLim',[0 Tend])
set(substrate_rels_ax,'XLim',[0 Tend])

% prepare print parameters
figs(2).fontsize = 6;
figs(2).print_position = [0 0 8 2];

%% plot slice of concentration over time
regions_types_all = arrayfun(@(i) [input(i).substrates(1).pars.regions_type],1:nsubcohorts);
rg_and_floor_ind = find(regions_types_all=="rectangular_grid" & [input.bv_location]=="floor");
if ~isempty(rg_and_floor_ind)
    n_rg = length(rg_and_floor_ind);
    C_line = zeros(nt,TME_size(end),n_rg,nsamps);
    out_rg = out(rg_and_floor_ind,:);

    figs(end+1).handle = figure;
    figs(end).name = "concentration_mesh";
    mesh_ax = gobjects(n_rg,1);
    nr = floor(sqrt(n_rg));
    nc = ceil(n_rg/nr);
    for i = 1:n_rg
        for si = 1:nsamps
            % since this blood vessel setup assumes that the n-1
            % dimensional floor (last dimension is at smallest value) is a
            % blood vessel, average across the 1:(n-1) dimensions. since
            % time is actually in the first dimension, this means use mean
            % on the 2:n dimensions
            if input(rg_and_floor_ind(i)).substrates(1).method=="global"
%                 Regions = setupRegions(input(rg_and_floor_ind(i)).substrates(1),input(i).TME_size,input(i).bv_location);
                C_line(:,:,i,si) = interp1(out_rg(i,si).t,squeeze(mean(reshape(out_rg(i,si).C{1}(:,input(rg_and_floor_ind(i)).substrates(1).solver.regions),[length(out_rg(i,si).t),size(input(rg_and_floor_ind(i)).substrates(1).solver.regions)]),2:length(TME_size))),t,'linear');
            else
                C_line(:,:,i,si) = interp1(out_rg(i,si).t,squeeze(mean(out_rg(i,si).C{1},2:length(TME_size))),t,'linear');
            end
        end

        z_data = mean(C_line(:,:,i,:),4);
        mesh_ax(i) = subplot(nr,nc,i);
        mesh(mesh_ax(i),t,1:TME_size(end),z_data')
        xlabel('t')
        if length(TME_size)==2
            ylabel('y')
        elseif length(TME_size)==3
            ylabel('z')
        end
        zlabel('C')
        title(disp_names(rg_and_floor_ind(i)))
    end
    set(mesh_ax,'FontSize',20,'XLim',t([1 end]),'YLim',[1,TME_size(end)])
    normalizeZLims(mesh_ax)
    set(mesh_ax,'View',[146,38])

    % prepare print parameters
    figs(3).fontsize = 6;
    figs(3).print_position = [0 0 4 3];
end



%% plot event data
event_names = {'prolifs','contact_inhibitions','spontaneous_apoptosis'};
nr = length(event_names);
nc = input(1).ntypes;
fig_temp = gobjects(2,1);
ax_event = gobjects(nr,nc,2);
pa_event_vals = gobjects(nsubcohorts,input(1).ntypes);
pa_event_rels = gobjects(nsubcohorts,input(1).ntypes);
for fi = 1:2
    fig_temp(fi) = figure;
    for ri = 1:nr
        for ci = 1:nc
            ax_event(ri,ci,fi) = subplot(nr,nc,r2c(nr,nc,[ri,ci]));
            hold on
        end
    end
end

V_event = zeros(nt,nsubcohorts,nsamps,2,length(event_names));
for ei = 1:length(event_names)
    for type_ind = 1:input(1).ntypes
        for i = 1:nsubcohorts
            for si = 1:nsamps
                if size(out(i,si).(event_names{ei}),2)>=type_ind
                    V_event(:,i,si,type_ind,ei) = interp1(out(i,si).t,out(i,si).(event_names{ei})(:,type_ind),t,'linear',out(i,si).(event_names{ei})(end,type_ind));
                end
            end

            if i==1
                vbar_abm = mean(V_event(:,1,:,type_ind,ei),3);
            end

            v = V_event(:,i,:,type_ind,ei);
            vbar = mean(v,3);
            dv = v - vbar_abm;
            rel_dv = 100*dv./vbar_abm;
            sd = std(dv,[],3);
            rel_sd = std(rel_dv,[],3);
            rel_dv_bar = mean(rel_dv,3);
            pa_event_vals(i,type_ind) = patch(ax_event(ei,type_ind,1),[t,flip(t)],reshape([vbar,flip(vbar)]+[-sd,flip(sd)],1,[]),cmap(i,:),...
                'FaceAlpha',0.2,'DisplayName',disp_names(i),'EdgeColor','none');
            pa_event_rels(i,type_ind) = patch(ax_event(ei,type_ind,2),[t,flip(t)],reshape([rel_dv_bar,flip(rel_dv_bar)]+[-rel_sd,flip(rel_sd)],1,[]),cmap(i,:),...
                'FaceAlpha',0.2,'DisplayName',disp_names(i),'EdgeColor','none');

            plot(ax_event(ei,type_ind,1),t,mean(V_event(:,i,:,type_ind,ei),3),'Color',cmap(i,:),'LineWidth',2);
            plot(ax_event(ei,type_ind,2),t(1:end-1),mean(diff(V_event(:,i,:,type_ind,ei))./NumCells(1:end-1,i,:,type_ind),3)/(t(2)-t(1)),'Color',cmap(i,:),'LineWidth',2);
            title(ax_event(ei,type_ind,1),regexprep(event_names{ei},'_',' '))
            title(ax_event(ei,type_ind,2),['relative ',regexprep(event_names{ei},'_',' ')])
%             ylabel(ax_event(ei,type_ind,2),["% Diff of above",sprintf(' from %s',disp_names(1))])
        end
    end

end

set(ax_event,'XLim',[0 Tend])
figs(end+1).handle = fig_temp(1);
figs(end).name = "events_cumulative";
figs(end+1).handle = fig_temp(2);
figs(end).name = "events_perstep";

% prepare print parameters
for fi = length(figs)-1:length(figs)
    figs(fi).fontsize = 6;
    figs(fi).print_position = [0 0 4 3];
end