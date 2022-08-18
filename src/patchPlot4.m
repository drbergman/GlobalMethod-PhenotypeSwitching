function patchPlot4(compartment_grid,ms,out,compartment_size,TME_size,make_svgs,delta_svg_num,cell_labels,num_timepoints_to_plot,method,ntypes,substrate_inds,substrates)

% allow for plotting drug concentration if present

warning("instead of passing in method (from the first substrate usually), should just pass in the desired file name format")

t = out.t;
nc = out.num_cells;
nc_per_c = out.cell_types_per_compartment;
ntypes = size(nc,2);
nsubstrates = length(substrate_inds);
if nsubstrates~=0
    substrate_names = regexprep(arrayify(substrates,"name"),'_',' ');
end

for si = length(substrate_inds):-1:1
    if substrates(substrate_inds(si)).method~="local" || ~isempty(substrates(substrate_inds(si)).solver.regions)
        Regions{si} = substrates(substrate_inds(si)).solver.regions;
    else
        Regions{si} = reshape(1:prod(TME_size),TME_size);
    end
end

if ~exist('make_svgs','var')
    make_svgs = false;
end
f = figure;
if ~isfield(out,'C') || isempty(substrate_inds)
    hm_plot = false;
    ax_sc = subplot(4,3,[1 9]);
    hold on;
    ax_ts = subplot(4,3,[10 12]);
    hold on;
    axis(ax_sc,'equal')
else
    hm_plot = true;
    nrows = 4;
    ncols = 3*(nsubstrates+1);
    ax_sc = subplot(nrows,ncols,[1 r2c(nrows,ncols,[3,3])]);
    hold on;
    ax_hm = gobjects(nsubstrates,1);
    for si = 1:nsubstrates
        ax_hm(si) = subplot(nrows,ncols,[4 + 3*(si-1), r2c(nrows,ncols,[3,3*(si+1)])]);
    end
    ax_ts = subplot(nrows,ncols,[3*(3*(nsubstrates+1))+1 nrows*ncols]);
    hold on;
    axis(ax_sc,'equal')
    axis(ax_hm,'equal')

end
colors = lines(ntypes);
compartment_start_x = (0:compartment_grid(1)-1)*ms(1);
compartment_start_y = (0:compartment_grid(2)-1)*ms(2);
if length(TME_size)==2
    [Cx,Cy] = ndgrid(compartment_start_x,compartment_start_y);
elseif length(TME_size)==3
    error('not set up for 3d yet')
    compartment_start_z = (0:compartment_grid(3)-1)*mz;
    [Cx,Cy,Cz] = ndgrid(compartment_start_x,compartment_start_y,compartment_start_z);
end

% if isfield(out,'num_ctl_per_site')
%     num_ctl_per_site = out.num_ctl_per_site;
% else
%     num_ctl_per_site = 1;
% end

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
L_now = squeeze(sliceof(nc_per_c,1,1));
for type_ind = ntypes:-1:1
%     p(:,:,type_ind) = sliceof(L_now,ndims(L_now),type_ind)/(compartment_size*num_agents_per_site(type_ind));
    p(:,:,type_ind) = (L_now==type_ind)/(compartment_size*num_agents_per_site(type_ind));
end
% p{1} = L_now(:,:,1)/compartment_size;
% p{2} = L_now(:,:,2)/(compartment_size*num_ctl_per_site);
% if ~same_lattice
%     p{1} = 0.5*p{1};
%     p{2} = 0.5*p{2};
% end

if ~same_lattice
    p = p / ntypes;
end

nonempty_compartments = any(L_now>0,3);
% if nc(1,2)>0
%     type2_ind = nonempty_compartments & p{2}>0;
%     rect_pl(2) = patch(ax_sc,Cx(type2_ind)'+ms(1)*(1+p{2}(type2_ind)'.*[-1;0;0;-1]),Cy(type2_ind)'+[ms(2)*ones(2,sum(type2_ind,'all'));zeros(2,sum(type2_ind,'all'))  ],colors(2,:),'EdgeColor','none');
% end
% if nc(1,1)>0
%     type1_ind = nonempty_compartments & p{1}>0;
%     rect_pl(1) = patch(ax_sc,Cx(type1_ind)'+ms(1)*(p{1}(type1_ind)'.*[ 0;1;1; 0]),Cy(type1_ind)'+[zeros(2,sum(type1_ind,'all'))  ;ms(2)*ones(2,sum(type1_ind,'all'))],colors(1,:),'EdgeColor','none');
% end

box_left_relative_start = zeros(size(nonempty_compartments));
% rect_pl = gobjects(ntypes,1);
for type_ind = 1:ntypes
    if nc(1,type_ind)>0
        this_type_ind = nonempty_compartments & p(:,:,type_ind)>0;
        p_this_type = p(:,:,type_ind);
        xcoords = (Cx(this_type_ind)'+box_left_relative_start(this_type_ind)')+ms(1)*p_this_type(this_type_ind)'.*[0;1;1;0];
        ycoords = Cy(this_type_ind)'+[ms(2)*ones(2,sum(this_type_ind,'all'));zeros(2,sum(this_type_ind,'all'))  ];
%         rect_pl(type_ind) = patch(ax_sc,xcoords,ycoords,colors(type_ind,:),'EdgeColor','none');
        patch(ax_sc,xcoords,ycoords,colors(type_ind,:),'EdgeColor','none');
        box_left_relative_start(this_type_ind) = box_left_relative_start(this_type_ind) + p_this_type(this_type_ind);
    else
%         rect_pl(type_ind) = patch(ax_sc,[],[],colors(type_ind,:),'EdgeColor','none');
    end
end

if hm_plot
    C_new = zeros([length(t),TME_size,nsubstrates]);
    hm_colors = colormap('bone');
    colormap(ax_sc,flipud(hm_colors))
    c = gobjects(nsubstrates,1);
    hm_alone = gobjects(nsubstrates,1);
    for si = 1:nsubstrates
        for ti = 1:length(t)
            C_new(ti,:,:,si) = reshape(out.C{substrate_inds(si)}(ti,Regions{si}),[1,TME_size]);
        end
        out.C{substrate_inds(si)} = C_new(:,:,:,si);

        hm_alone(si) = imagesc(ax_hm(si),((1:compartment_grid(1))-0.5)*ms(1),((1:compartment_grid(2))-0.5)*ms(2),permute(out.C{substrate_inds(si)}(1,:,:),[3,2,1]),'AlphaData',0.9);

        colormap(ax_hm(si),flipud(hm_colors))
        c(si) = colorbar(ax_hm(si));
%         c(si).Label.String = substrates(substrate_inds(si)).name;
        c(si).Label.String = substrate_names(si);
        if nsubstrates==1
            c(si).Location = 'east';
            c(si).Label.VerticalAlignment = 'baseline';
            c(si).Label.Rotation = 270;
            c(si).Position(1) = c(si).Position(1) + c(si).Position(3);
            ax_hm(si).Position(1) = ax_hm(si).Position(1) - c(si).Position(3);
            c(si).AxisLocation = "out";
        else
            c(si).Location = "northoutside";
            c(si).AxisLocation = "out";

        end
        %     caxis(ax_hm(si),[0 max(out.C{substrate_inds(si)},[],'all')])
        %     caxis(ax_sc,[0 max(out.C{substrate_inds(si)},[],'all')])
        max_c = max(removeslice(out.C{substrate_inds(si)},1,1),[],'all');
        if max_c==0
            max_c = 1;
        end
        caxis(ax_hm(si),[0 max_c])
        xlim(ax_hm(si),[0 TME_size(1)])
        ylim(ax_hm(si),[0 TME_size(2)])
        ax_hm(si).YDir = 'normal';
        ax_hm(si).FontSize = 20;
        ax_hm(si).XTick = [];
        ax_hm(si).YTick = [];

    end
    if nsubstrates==1
        ax_sc.Position(1) = ax_sc.Position(1) - c(si).Position(3);
        ax_ts.Position(1) = ax_ts.Position(1) - c(si).Position(3);
    end
    max_c = max(removeslice(out.C{substrate_inds(si)},1,1),[],'all');
    if max_c==0
        max_c = 1;
    end
    caxis(ax_sc,[0 max_c])
    set(ax_sc.Children,'FaceAlpha',0.6);
    if nsubstrates==1
        imagesc(ax_sc,((1:compartment_grid(1))-0.5)*ms(1),((1:compartment_grid(2))-0.5)*ms(2),permute(out.C{substrate_inds(si)}(1,:,:),[3,2,1]),'AlphaData',0.9);
        ax_sc.Children = flip(ax_sc.Children);
    end
    ax_sc.XTick = [];
    ax_sc.YTick = [];

    f.Position(3:4) = 2*f.Position(3:4);
    f.Position(1:2) = 0;

    axis(ax_sc,'equal')
    axis(ax_hm,'equal')
%     axis(ax_hm,'square')
end

for i = ntypes:-1:1
    ts_pl(i) = plot(ax_ts,t,nc(:,i),'Color',colors(i,:),'LineWidth',2);
    ts_m(i) = scatter(ax_ts,t(1),nc(1,i),20,colors(i,:),'filled');
end
fig_title = title(ax_sc,sprintf('t = %3.2f',t(1)));
ax_sc.FontSize = 20;
ax_ts.FontSize = 20;
axis(ax_sc,[0 TME_size(1) 0 TME_size(2)])
if hm_plot
    drawnow; % make sure resizing the ax_sc axes above takes effect before setting this
    for si = 1:nsubstrates
        ax_hm(si).PlotBoxAspectRatio = ax_sc.PlotBoxAspectRatio;
        if nsubstrates>1
            ax_hm(si).Position(2) = .355;
        end
    end
end
ax_ts.XLim = [0 t(end)];
ax_ts.XAxis.Label.String = "Time (min)";
ax_ts.YAxis.Label.String = "Cell Count";
ax_ts.YLim(1) = 0;
legend(ax_ts,ts_pl,regexprep(cell_labels,"_"," "),'Location','best')
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
        %     p(:,:,type_ind) = sliceof(L_now,ndims(L_now),type_ind)/(compartment_size*num_agents_per_site(type_ind));
        p(:,:,type_ind) = (L_now==type_ind)/(compartment_size*num_agents_per_site(type_ind));
    end

%     p{1} = L_now(:,:,1)/compartment_size;
%     p{2} = L_now(:,:,2)/(compartment_size*num_ctl_per_site);
%     if ~same_lattice
%         p{1} = 0.5*p{1};
%         p{2} = 0.5*p{2};
%     end
%     for type_ind = ntypes:-1:1
%         p(:,:,type_ind) = L_now(:,:,type_ind)/(compartment_size*num_agents_per_site(type_ind));
%     end
    if ~same_lattice
        p = p / ntypes;
    end
    nonempty_compartments = any(L_now>0,3);
    delete(ax_sc.Children)
%     if nc(tinds(ti),2)>0
%         type2_ind = nonempty_compartments & p{2}>0;
%         rect_pl(2) = patch(ax_sc,Cx(type2_ind)'+ms(1)*(1+p{2}(type2_ind)'.*[-1;0;0;-1]),Cy(type2_ind)'+[ms(2)*ones(2,sum(type2_ind,'all'));zeros(2,sum(type2_ind,'all'))  ],colors(2,:),'EdgeColor','none');
%     end
%     if nc(tinds(ti),1)>0
%         type1_ind = nonempty_compartments & p{1}>0;
%         rect_pl(1) = patch(ax_sc,Cx(type1_ind)'+ms(1)*(  p{1}(type1_ind)'.*[ 0;1;1; 0]),Cy(type1_ind)'+[zeros(2,sum(type1_ind,'all'))  ;ms(2)*ones(2,sum(type1_ind,'all'))],colors(1,:),'EdgeColor','none');
%     end

    box_left_relative_start = zeros(size(nonempty_compartments));
    for type_ind = 1:ntypes
        if nc(tinds(ti),type_ind)>0
            this_type_ind = nonempty_compartments & p(:,:,type_ind)>0;
            p_this_type = p(:,:,type_ind);
            xcoords = (Cx(this_type_ind)'+box_left_relative_start(this_type_ind)')+ms(1)*p_this_type(this_type_ind)'.*[0;1;1;0];
            ycoords = Cy(this_type_ind)'+[ms(2)*ones(2,sum(this_type_ind,'all'));zeros(2,sum(this_type_ind,'all'))  ];
%             rect_pl(type_ind) = patch(ax_sc,xcoords,ycoords,colors(type_ind,:),'EdgeColor','none');
            patch(ax_sc,xcoords,ycoords,colors(type_ind,:),'EdgeColor','none');
            box_left_relative_start(this_type_ind) = box_left_relative_start(this_type_ind) + p_this_type(this_type_ind);
        end
    end

    if hm_plot
        for si = 1:nsubstrates
            hm_alone(si).CData = permute(out.C{substrate_inds(si)}(tinds(ti),:,:),[3,2,1]);
        end
        set(ax_sc.Children,'FaceAlpha',0.6);
        if nsubstrates==1
            imagesc(ax_sc,((1:compartment_grid(1))-0.5)*ms(1),((1:compartment_grid(2))-0.5)*ms(2),permute(out.C{substrate_inds(si)}(tinds(ti),:,:),[3,2,1]),'AlphaData',0.9);
            ax_sc.Children = flip(ax_sc.Children);
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
        print(sprintf('output/snapshot_%s_%08d',method,svg_num),'-djpeg')
        svg_num = svg_num+1;
        next_svg_num = next_svg_num + delta_svg_num;
    end

end