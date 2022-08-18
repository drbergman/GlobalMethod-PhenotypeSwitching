function scatterPlot(compartment_grid,ms,out,compartment_size,TME_size,make_svgs,delta_svg_num,cell_labels,num_timepoints_to_plot,method,cutaway)

% allow for plotting drug concentration if present

t = out.t;
nc = out.num_cells;
nc_per_c = out.cell_types_per_compartment;
ntypes = size(nc,2);
if ~exist('make_svgs','var')
    make_svgs = false;
end
f = figure;
f.Units = "normalized";
f.Position(1:2) = 0;
f.Position(3:4) = [0.3 0.4];
ax_sc = subplot(4,3,[1 9]);
hold on;
ax_ts = subplot(4,3,[10 12]);
hold on;
axis(ax_sc,'equal')

colors = lines(ntypes);
compartment_start_x = (0:compartment_grid(1)-1)*ms(1);
compartment_start_y = (0:compartment_grid(2)-1)*ms(2);
if length(TME_size)==2
    [Cx,Cy] = ndgrid(compartment_start_x,compartment_start_y);
elseif length(TME_size)==3
    compartment_start_z = (0:compartment_grid(3)-1)*ms(3);
    [Cx,Cy,Cz] = ndgrid(compartment_start_x,compartment_start_y,compartment_start_z);
end

if isfield(out,'num_agents_per_site')
    num_agents_per_site = out.num_agents_per_site;
else
    num_agents_per_site = ones(ntypes,1);
end

if isfield(out,'same_lattice')
    same_lattice = out.same_lattice;
else
    same_lattice = true;
end

if ~isempty(cutaway)
    switch cutaway
        case "none"
            show_these = true(TME_size);
        case "octant"
            show_these = true(TME_size);
            show_these(1:floor(TME_size(1)/2),1:floor(TME_size(2)/2),ceil(TME_size(3)/2):TME_size(3)) = false;
    end
else
    show_these = true(TME_size);
end

L_now = squeeze(sliceof(nc_per_c,1,1));
for type_ind = ntypes:-1:1
    if length(TME_size)==2
        p(:,:,type_ind) = (L_now==type_ind)/(compartment_size*num_agents_per_site(type_ind));
    elseif length(TME_size)==3
        p(:,:,:,type_ind) = (L_now==type_ind)/(compartment_size*num_agents_per_site(type_ind));
    end
end

if ~same_lattice
    p = p / ntypes;
end

for type_ind = 1:ntypes
    if nc(1,type_ind)>0
        temp_ind = show_these & sliceof(p,ndims(p),type_ind)>0;
        sc_pl(type_ind) = scatter3(ax_sc,Cx(temp_ind),Cy(temp_ind),Cz(temp_ind),200,colors(type_ind,:),'filled','Marker','o','MarkerEdgeColor','k','LineWidth',1);
    else
        sc_pl(type_ind) = scatter3(ax_sc,[],[],[],200,colors(type_ind,:),'filled','Marker','o','MarkerEdgeColor','k','LineWidth',1);
    end
end
az = -45;
el = 25;
view(ax_sc,az,el)

for i = ntypes:-1:1
    ts_pl(i) = plot(ax_ts,t,nc(:,i),'Color',colors(i,:),'LineWidth',2);
    ts_m(i) = scatter(ax_ts,t(1),nc(1,i),20,colors(i,:),'filled');
end
fig_title = title(ax_sc,sprintf('t = %3.2f',t(1)));
ax_sc.FontSize = 20;
ax_ts.FontSize = 20;
axis(ax_sc,[0 TME_size(1) 0 TME_size(2)])

ax_ts.XLim = [0 t(end)];
ax_ts.XAxis.Label.String = "Time (min)";
ax_ts.YAxis.Label.String = "Cell Count";
ax_ts.YLim(1) = 0;
legend(ax_ts,ts_pl,cell_labels,'Location','best')
drawnow;

if make_svgs
    svg_num = 0;
    print(sprintf('output/snapshot_%s_%08d',method,svg_num),'-dsvg')
    svg_num = 1;
    next_svg_num = delta_svg_num;
end

tinds = round(linspace(1,length(t),num_timepoints_to_plot));
tinds = unique(tinds);

for ti = 2:length(tinds)
    L_now = squeeze(sliceof(nc_per_c,1,tinds(ti)));
    for type_ind = ntypes:-1:1
        if length(TME_size)==2
            p(:,:,type_ind) = (L_now==type_ind)/(compartment_size*num_agents_per_site(type_ind));
        elseif length(TME_size)==3
            p(:,:,:,type_ind) = (L_now==type_ind)/(compartment_size*num_agents_per_site(type_ind));
        end
    end

    if ~same_lattice
        p = p / ntypes;
    end
%     delete(ax_sc.Children)

    for type_ind = 1:ntypes
        if nc(tinds(ti),type_ind)>0
            temp_ind = show_these & sliceof(p,ndims(p),type_ind)>0;
            sc_pl(type_ind).XData = Cx(temp_ind);
            sc_pl(type_ind).YData = Cy(temp_ind);
            sc_pl(type_ind).ZData = Cz(temp_ind);
        else
            sc_pl(type_ind).XData = [];
            sc_pl(type_ind).YData = [];
            sc_pl(type_ind).ZData = [];
        end

    end


    for i = 1:ntypes
        ts_m(i).XData = t(tinds(ti));
        ts_m(i).YData = nc(tinds(ti),i);
    end
    fig_title.String = sprintf('t = %3.2f',t(tinds(ti)));
    drawnow;
    if make_svgs && tinds(ti) >= next_svg_num+1
        print(sprintf('output/snapshot_%s_%08d',method,svg_num),'-dsvg')
        svg_num = svg_num+1;
        next_svg_num = next_svg_num + delta_svg_num;
    end

end