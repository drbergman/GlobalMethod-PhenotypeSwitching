function FigurePanels(out,input,codensity_computation_on)

%%
TME_size = input(1).TME_size;
Tend = input(1).Tend;
save_dt = input(1).save_dt;
nsamps = numel(out)/numel(input);
N = prod(TME_size);
nt = round(Tend/(save_dt)+1);
t = linspace(0,Tend,nt)/(60*24);
if isfield(input(1).substrates,"oxygen")
    for i = 1:numel(input)
        input(i).substrates.method = input(i).substrates.oxygen.method;
    end
end
if ~iscell(out(1).C)
    for i = 1:numel(out)
        out(i).C = {out(i).C};
    end
end
[NumCells,InternSubstrates,AllConcentrations] = normalizeOutputPhenotypeSwitching(out,input,size(out,1),nsamps);
for mi = 1:numel(input)
    disp_names(mi) = input(mi).substrates(1).method;
end

%%
% method_color = cubehelix(2,2,-1.5,1,1,[0,.6],[.4,.6]);
method_color = [ 0.0707         0    0.2163
    0.3854    0.7045    0.6247];
font_size = 6;
res_factor = 6;
single_panel_fig_pos = [0 0 1 0.75];
one_x_two_panel_fig_pos = [0 0 2 0.75];
two_x_two_panel_fig_pos = [0 0 2 1.5];
f = gobjects(0,1);

%% plot healthy populations relative to local method
f(end+1) = figure("Name","healthy_rel");
ax = gca;
hold on;


method_means = mean(NumCells(:,:,:,1),3);
d_numcells = NumCells(:,:,:,1) - method_means(:,1,:);
rel_d = 100*d_numcells./method_means(:,1,1);
rel_sd = std(rel_d,[],3);
rel_d_bar = mean(rel_d,3);

pa_rels = gobjects(size(out,1),1);
ls_rels = gobjects(size(out,1),1);
for mi = 1:size(out,1)
    pa_rels(mi) = patch(ax,[t,flip(t)],reshape([rel_d_bar(:,mi,1),flip(rel_d_bar(:,mi,1))]+[-rel_sd(:,mi,1),flip(rel_sd(:,mi,1))],1,[]),method_color(mi,:),...
        'FaceAlpha',0.2,'DisplayName',disp_names(mi),'EdgeColor','none');

    ls_rels(mi) = plot(ax,t,rel_d_bar(:,mi,1),'Color',method_color(mi,:),'LineWidth',2,'DisplayName',disp_names(mi));
end

legend(pa_rels,"Location","best")

xlabel('Time (d)')
ylabel(["Percent","Difference"])
title(["Proliferating","Tumor Cells"])

ax.FontSize = font_size*res_factor;
f(end).Units = "inches";
f(end).Position = single_panel_fig_pos*res_factor;

%% plot quiescent populations relative to local method
f(4) = figure("Name","quiesc_rel");
ax = gca;
hold on;

method_means = mean(NumCells(:,:,:,2),3);
d_numcells = NumCells(:,:,:,2) - method_means(:,1,:);
rel_d = 100*d_numcells./method_means(:,1,1);
rel_d(isnan(rel_d)) = 0; % when there are no quiescent cells, this relative difference is nan
rel_sd = std(rel_d,[],3);
rel_d_bar = mean(rel_d,3);

pa_rels = gobjects(size(out,1),1);
ls_rels = gobjects(size(out,1),1);
for mi = 1:size(out,1)
    pa_rels(mi) = patch(ax,[t,flip(t)],reshape([rel_d_bar(:,mi,1),flip(rel_d_bar(:,mi,1))]+[-rel_sd(:,mi,1),flip(rel_sd(:,mi,1))],1,[]),method_color(mi,:),...
        'FaceAlpha',0.2,'DisplayName',disp_names(mi),'EdgeColor','none');

    ls_rels(mi) = plot(ax,t,rel_d_bar(:,mi,1),'Color',method_color(mi,:),'LineWidth',2,'DisplayName',disp_names(mi));
end

legend(pa_rels,"Location","best")

xlabel('Time (d)')
ylabel(["Percent","Difference"])
title(["Quiescent","Tumor Cells"])

ax.FontSize = font_size*res_factor;
f(4).Units = "inches";
f(4).Position = single_panel_fig_pos*res_factor;

%% plot wall times
f(5) = figure("Name","wall_times");
ax = gca;
hold on;

wall_times = arrayify(out,'wall_time',1);
for mi = 1:size(out,1)
    histogram(wall_times(mi,:)/60,'FaceColor',method_color(mi,:),'DisplayName',disp_names(mi))
end

legend("Location","north")

xlabel('Wall Time (min)')
ylabel("Frequency")
title("Wall Time by Method")

ax.FontSize = font_size*res_factor;
f(5).Units = "inches";
f(5).Position = single_panel_fig_pos*res_factor;

    %% mean y values
f(7) = figure("Name","mean_y_vals_patch");
ax = gca;
hold on;
ybar = zeros(nt,nsamps,2);
pa = gobjects(size(out,1),1);
for mi = 1:size(out,1)
    for si = 1:nsamps
        y_temp = zeros(length(out(mi,si).t),1);
        for ti = 1:length(out(mi,si).t)
            [~,yy] = ind2sub(input(mi).TME_size,find(out(mi,si).cell_types_per_compartment(ti,:,:)==1));
            
            y_temp(ti) = mean(yy);
        end
        ybar(:,si,mi) = interp1(out(mi,si).t/(60*24),y_temp,t);
    end
    ybar(:,:,mi) = ybar(:,:,mi)*input(mi).lattice_h;

    sd = std(ybar(:,:,mi),[],2);
    method_mean = mean(ybar(:,:,mi),2);
    pa(mi) = patch([t,flip(t)],[method_mean-sd;flip(method_mean+sd)]',method_color(mi,:),'FaceAlpha',0.2,'EdgeColor','none','DisplayName',disp_names(mi));
    plot(t,method_mean,'Color',method_color(mi,:),'LineWidth',2)
    

end

legend(pa,'Location','northeast')



xlabel('Time (d)')
ylabel("Distance (\mum)")
title(["Mean Distance","to Vasculature"])

ax.FontSize = font_size*res_factor;
f(7).Units = "inches";
f(7).Position = single_panel_fig_pos*res_factor;


%% mean y values - quiescent
f(8) = figure("Name","mean_y_vals_patch_quiesc");
ax = gca;
hold on;
ybar = zeros(nt,nsamps,2);
pa = gobjects(size(out,1),1);
for mi = 1:size(out,1)
    for si = 1:nsamps
        y_temp = zeros(length(out(mi,si).t),1);
        for ti = 1:length(out(mi,si).t)
            [~,yy] = ind2sub(size(out(mi,si).cell_types_per_compartment,2:3),find(out(mi,si).cell_types_per_compartment(ti,:,:)==2));
            
            y_temp(ti) = mean(yy);
        end
        ybar(:,si,mi) = interp1(out(mi,si).t/(60*24),y_temp,t);
    end
    ybar(:,:,mi) = ybar(:,:,mi)*input(mi).lattice_h;

    sd = std(ybar(:,:,mi),[],2);
    method_mean = mean(ybar(:,:,mi),2);
    pa(mi) = patch([t,flip(t)],[method_mean-sd;flip(method_mean+sd)]',method_color(mi,:),'FaceAlpha',0.2,'EdgeColor','none','DisplayName',disp_names(mi));
    plot(t,method_mean,'Color',method_color(mi,:),'LineWidth',2)
    

end

legend(pa,'Location','northeast')



xlabel('Time (d)')
ylabel("Distance (\mum)")
title(["Mean Distance","to Vasculature"])

ax.FontSize = font_size*res_factor;
f(8).Units = "inches";
f(8).Position = single_panel_fig_pos*res_factor;


%% internalized substrate healthy

f(9) = figure("Name","healthy_intern_substrate");
ax = gca;
hold on;

for mi = 1:size(out,1)
    for si = 1:nsamps
        plot(t,InternSubstrates(:,mi,si,1),'color',method_color(mi,:),'LineWidth',1)
    end
end

xlabel('Time (d)')
ylabel('Concentration')
title(["Proliferating","Internalized Oxygen"])

ax.FontSize = font_size*res_factor;
f(9).Units = "inches";
f(9).Position = single_panel_fig_pos*res_factor;

%% plot all quiescent tumor cell populations
f(10) = figure("Name","quiesc_intern_substrate");
ax = gca;
hold on;

for mi = 1:size(out,1)
    for si = 1:nsamps
        plot(t,InternSubstrates(:,mi,si,2),'color',method_color(mi,:),'LineWidth',1)
    end
end

xlabel('Time (d)')
ylabel('Concentration')
title(["Quiescent","Internalized Oxygen"])

ax.FontSize = font_size*res_factor;
f(10).Units = "inches";
f(10).Position = single_panel_fig_pos*res_factor;

%% plot healthy populations relative to local method
f(11) = figure("Name","healthy_intern_substrate_rel");
ax = gca;
hold on;


method_means = mean(InternSubstrates(:,:,:,1),3);
d_intern_sub = InternSubstrates(:,:,:,1) - method_means(:,1,:);
rel_d = 100*d_intern_sub./method_means(:,1,1);
rel_sd = std(rel_d,[],3);
rel_d_bar = mean(rel_d,3);

pa_rels = gobjects(size(out,1),1);
ls_rels = gobjects(size(out,1),1);
for mi = 1:size(out,1)
    pa_rels(mi) = patch(ax,[t,flip(t)],reshape([rel_d_bar(:,mi,1),flip(rel_d_bar(:,mi,1))]+[-rel_sd(:,mi,1),flip(rel_sd(:,mi,1))],1,[]),method_color(mi,:),...
        'FaceAlpha',0.2,'DisplayName',disp_names(mi),'EdgeColor','none');

    ls_rels(mi) = plot(ax,t,rel_d_bar(:,mi,1),'Color',method_color(mi,:),'LineWidth',2,'DisplayName',disp_names(mi));
end

legend(pa_rels,"Location","best")

xlabel('Time (d)')
ylabel(["Percent","Difference"])
title(["Proliferating","Internalized Oxygen"])

ax.FontSize = font_size*res_factor;
f(11).Units = "inches";
f(11).Position = single_panel_fig_pos*res_factor;

%% plot quiescent populations relative to local method
f(12) = figure("Name","quiesc_intern_substrate_rel");
ax = gca;
hold on;

method_means = mean(InternSubstrates(:,:,:,2),3);
d_intern_sub = InternSubstrates(:,:,:,2) - method_means(:,1,:);
rel_d = 100*d_intern_sub./method_means(:,1,1);
rel_d(isnan(rel_d)) = 0; % when there are no quiescent cells, this relative difference is nan
rel_sd = std(rel_d,[],3);
rel_d_bar = mean(rel_d,3);

pa_rels = gobjects(size(out,1),1);
ls_rels = gobjects(size(out,1),1);
for mi = 1:size(out,1)
    pa_rels(mi) = patch(ax,[t,flip(t)],reshape([rel_d_bar(:,mi,1),flip(rel_d_bar(:,mi,1))]+[-rel_sd(:,mi,1),flip(rel_sd(:,mi,1))],1,[]),method_color(mi,:),...
        'FaceAlpha',0.2,'DisplayName',disp_names(mi),'EdgeColor','none');

    ls_rels(mi) = plot(ax,t,rel_d_bar(:,mi,1),'Color',method_color(mi,:),'LineWidth',2,'DisplayName',disp_names(mi));
end

legend(pa_rels,"Location","best")

xlabel('Time (d)')
ylabel(["Percent","Difference"])
title(["Quiescent","Internalized Oxygen"])

ax.FontSize = font_size*res_factor;
f(12).Units = "inches";
f(12).Position = single_panel_fig_pos*res_factor;

%% average concentration in microenvironment
f(13) = figure("Name","avg_conc");
ax = gca;
hold on;

method_means = mean(AllConcentrations(:,:,:,2),3);
method_sd = std(AllConcentrations(:,:,:,2),[],3);

for mi = 1:size(out,1)
    patch(ax,[t,flip(t)],[method_means(:,mi) + method_sd(:,mi);flip(method_means(:,mi) - method_sd(:,mi))],method_color(mi,:),...
        'FaceAlpha',0.2,'DisplayName',disp_names(mi),'EdgeColor','none');
    plot(t,method_means(:,mi),'color',method_color(mi,:),'LineWidth',1)
end

xlabel('Time (d)')
ylabel('Concentration')
title(["Average","Free Oxygen"])

ax.FontSize = font_size*res_factor;
f(13).Units = "inches";
f(13).Position = single_panel_fig_pos*res_factor;

%% average concentration in microenvironment relative to local average
f(14) = figure("Name","avg_conc_rel");
ax = gca;
hold on;

method_means = mean(AllConcentrations(:,:,:,2),3);
d_avg_conc = AllConcentrations(:,:,:,2) - method_means(:,1);
rel_d = 100*d_avg_conc./method_means(:,1,1);
rel_d(isnan(rel_d)) = 0; % when there are no quiescent cells, this relative difference is nan
rel_sd = std(rel_d,[],3);
rel_d_bar = mean(rel_d,3);

pa_rels = gobjects(size(out,1),1);
ls_rels = gobjects(size(out,1),1);
for mi = 1:size(out,1)
    pa_rels(mi) = patch(ax,[t,flip(t)],reshape([rel_d_bar(:,mi,1),flip(rel_d_bar(:,mi,1))]+[-rel_sd(:,mi,1),flip(rel_sd(:,mi,1))],1,[]),method_color(mi,:),...
        'FaceAlpha',0.2,'DisplayName',disp_names(mi),'EdgeColor','none');

    ls_rels(mi) = plot(ax,t,rel_d_bar(:,mi,1),'Color',method_color(mi,:),'LineWidth',2,'DisplayName',disp_names(mi));
end

legend(pa_rels,"Location","best")

xlabel('Time (d)')
ylabel(["Percent","Difference"])
title(["Average","Free Oxygen"])

ax.FontSize = font_size*res_factor;
f(14).Units = "inches";
f(14).Position = single_panel_fig_pos*res_factor;

%% plot event data
last_fig_num = 14;
event_names = {'prolifs','contact_inhibitions','spontaneous_apoptosis'};
title_names = ["Proliferations","Contact Inhibitions","Spontaneous Apoptosis"];
nfigs = length(event_names);
fig_temp = gobjects(nfigs,2);
ax_event = gobjects(nfigs,2);
pa_event_vals = gobjects(numel(input),1);
pa_event_rels = gobjects(numel(input),1);
for ei = 1:nfigs
    for vi = 1:2 % value index
        fig_temp(ei,vi) = figure;
        fig_temp(ei,vi).Name = event_names{ei};
        if vi == 2
            fig_temp(ei,vi).Name = [fig_temp(ei,vi).Name,'_rel'];
        end
        ax_event(ei,vi) = gca;
        hold on
    end
end

V_event = zeros(nt,numel(input),nsamps,length(event_names));
for ei = 1:length(event_names)
    for mi = 1:numel(input)
        for si = 1:nsamps
            if size(out(mi,si).(event_names{ei}),2)>=1
                V_event(:,mi,si,ei) = interp1(out(mi,si).t,out(mi,si).(event_names{ei})(:,1),t*24*60,'linear',out(mi,si).(event_names{ei})(end,1));
            end
        end
        method_means(:,mi,ei) = mean(V_event(:,mi,:,ei),3);
    end
end

for ei = 1:length(event_names)
    for mi = 1:numel(input)
        dv = V_event(:,mi,:,ei) - method_means(:,mi,ei);
        rel_dv = 100*dv./method_means(:,1,ei);
        rel_dv(isnan(rel_dv)) = 0; % when there are no quiescent cells, this relative difference is nan
        sd = std(V_event(:,mi,:,ei),[],3);
        rel_sd = std(rel_dv,[],3);
        rel_dv_bar = mean(rel_dv,3);
        pa_event_vals(mi) = patch(ax_event(ei,1),[t,flip(t)],reshape([method_means(:,mi,ei),flip(method_means(:,mi,ei))]+[-sd,flip(sd)],1,[]),method_color(mi,:),...
            'FaceAlpha',0.2,'DisplayName',disp_names(mi),'EdgeColor','none');
        pa_event_rels(mi) = patch(ax_event(ei,2),[t,flip(t)],reshape([rel_dv_bar,flip(rel_dv_bar)]+[-rel_sd,flip(rel_sd)],1,[]),method_color(mi,:),...
            'FaceAlpha',0.2,'DisplayName',disp_names(mi),'EdgeColor','none');

        plot(ax_event(ei,1),t,method_means(:,mi,ei),'Color',method_color(mi,:),'LineWidth',2);
        plot(ax_event(ei,2),t,rel_dv_bar,'Color',method_color(mi,:),'LineWidth',2);
        title(ax_event(ei,1),title_names(ei))
        title(ax_event(ei,2),["Percent Difference in",title_names(ei)])
    end
    legend(ax_event(ei,1),pa_event_vals,'location','best')
    legend(ax_event(ei,2),pa_event_rels,'location','best')
end

set(ax_event,'FontSize',font_size*res_factor);
for fi = 1:numel(fig_temp)
    f(last_fig_num+fi) = fig_temp(fi);
    f(last_fig_num+fi).Units = "inches";
    f(last_fig_num+fi).Position = single_panel_fig_pos*res_factor;
end

%% heterogeneity plot
if length(input(1).TME_size)==2 % only do this for the 2D runs
        f(end+1) = figure("Name","heterogeneity_plot");
        ax(2) = subplot(1,2,2);
        ax(1) = subplot(1,2,1);
        N = zeros(nt,TME_size(2),2,2,nsamps);
        for mi = 1:size(out,1)
            for si = 1:nsamps
                for type_ind = 1:2
                    N(:,:,type_ind,mi,si) = sum(out(mi,si).cell_types_per_compartment==type_ind,2);
                end
            end
        end
        N = permute(N,[1,2,5,3,4]);
        N = reshape(N,[],2,2);

        for mi = 1:size(out,1)
            ind = N(:,1,mi)==0 | N(:,2,mi)==0;
            N(ind,:,mi) = NaN;
        end
        for mi = 1:size(out,1)
            ind = ~isnan(N(:,1,mi));

            histogram2(ax(mi),N(ind,1,mi),N(ind,2,mi),0:100,0:100,'DisplayStyle','tile','EdgeColor','none');

        end

        xlabel(ax,'Proliferating')
        ylabel(ax,"Quiescent")

        set(ax,'FontSize',font_size*res_factor)
        f(end).Units = "inches";
        f(end).Position = one_x_two_panel_fig_pos*res_factor;
end

%% codensities


if codensity_computation_on
    
ntypes = input(1).ntypes;
tt = unique(round(linspace(1,nt,100)));
ncod = length(tt);
COD = zeros(ncod,ntypes,ntypes,2,nsamps);

for mi = 1:size(out,1)
    for si = 1:nsamps
        NCPC = out(mi,si).cell_types_per_compartment;
        NCPC = sliceof(NCPC,1,tt);
        for ti = 1:ncod
            COD(ti,:,:,mi,si) = codensityFunction(squeeze(NCPC(ti,:,:,:)),10,ntypes);
        end
    end
end
%%
f(end+1) = figure("Name","codensities");
for this_i = 1:ntypes
    for other_i = 1:ntypes
        ax(other_i,this_i) = subplot(ntypes,ntypes,r2c(ntypes,ntypes,[other_i,this_i]));
        hold on;
    end
end
for mi = 1:size(out,1)
    for oi = 1:2
        for ti = 1:2
            [x,y_mean,pc] = my_patchPlot(t(tt),squeeze(COD(:,oi,ti,mi,:)));
            nan_log = isnan(pc{2});
            pc{1}(nan_log) = [];
            pc{2}(nan_log) = [];
            patch(ax(oi,ti),pc{:},method_color(mi,:),'FaceAlpha',0.2,'DisplayName',disp_names(mi),'EdgeColor','none');
            plot(ax(oi,ti),x,y_mean,'color',method_color(mi,:),'LineWidth',1,'displayname',disp_names(mi));
            xlabel(ax(oi,ti),"Time (d)")
        end
    end
end
title(ax(1,1),["Codensity of","Proliferating"])
title(ax(1,2),["Codensity of","Quiescent"])
ylabel(ax(1,1),["Relative to","Proliferating"],"FontWeight","bold");
ylabel(ax(2,1),["Relative to","Quiescent"],"FontWeight","bold");

set(ax,'FontSize',font_size*res_factor)
f(end).Units = "inches";
f(end).Position = two_x_two_panel_fig_pos*res_factor;

end