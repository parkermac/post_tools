function compare_roms_obs_z_plot(Obs, Roms, Forcing)

n = length(Obs); %number of stations

%-----%-----%-----%-----%-----%-----%-----%-----%-----%-----%-----%-----%-----%-----
%----- plot T/S comparisons -------------------------------------------------------- 
%%% plot
 figure; set(gcf,'position',[177   207   809   731]);
 for i = 1:n
    subplot(n,2,2*i-1)  %plot salt 
        plot(Obs(i).td, Obs(i).salt, 'bo-', 'linewidth', 1.5)
        hold on; grid on;
        plot(Roms(i).td, Roms(i).salt, 'kv-', 'linewidth', 1.5)
        ax = axis; axis([Obs(i).td(1)-1 Obs(i).td(end)+1 ax(3) ax(4)]);
        set(gca,'xtick',Obs(i).td(1):7:Obs(i).td(end));
        datetick('x',6,'keeplimits','keepticks');
        title([num2str(Obs(i).tdepth) 'm Sal at : ' num2str(Obs(i).pos(1),'%5.2f') 'W, ' ...
            num2str(Obs(i).pos(2),'%5.2f') 'N'])
        ylabel('Salinity');
    subplot(n,2,2*i)  %plot temp 
        plot(Obs(i).td, Obs(i).temp, 'ro-', 'linewidth', 1.5)
        hold on; grid on;
        plot(Roms(i).td, Roms(i).temp, 'kv-', 'linewidth', 1.5)
        ax2 = axis; axis([Obs(i).td(1)-1 Obs(i).td(end)+1 ax2(3) ax2(4)]);
        set(gca,'xtick',Obs(i).td(1):7:Obs(i).td(end));
        datetick('x',6,'keeplimits','keepticks');
        title([num2str(Obs(i).tdepth) 'm Temp at : ' num2str(Obs(i).pos(1),'%5.2f') 'W, ' ...
            num2str(Obs(i).pos(2),'%5.2f') 'N'])
        ylabel('temperature (\circC)');
 end
%-----%-----%-----%-----%-----%-----%-----%-----%-----%-----%-----%-----%-----%-----

%-----%-----%-----%-----%-----%-----%-----%-----%-----%-----%-----%-----%-----%-----
%----- plot forcing ----------------------------------------------------
colors=['b-';'r-';'g-';'m-';'c-';'k-'];
 figure; set(gcf,'position',[214   501   892   412]);
    subplot(211)  %plot river forcing
        plot(Forcing.rivertime, Forcing.riverQ, 'k', 'linewidth', 2)
        hold on; grid on;
        set(gca,'xtick',Forcing.rivertime(1):7:Forcing.rivertime(end));
        datetick('x',6,'keeplimits','keepticks');
        ax = axis; axis(ax);
        box = [Obs(1).td(1) Obs(1).td(1) Obs(1).td(end) Obs(1).td(end) Obs(1).td(1)];
        plot(box, [ax(3) ax(4) ax(4) ax(3) ax(3)], 'b','linewidth',1.5)
        ylabel('Q (m^3 s^{-1})'); 
        title([Forcing.rivername ' River discharge'],'fontsize',14);
    subplot(212)  %plot N/S wind stress at pos 2 above
        plot(Forcing.time, Forcing.vstress, 'k', 'linewidth', 2.5)
        hold on; grid on;
        plot([Forcing.time(1) Forcing.time(end)],[0 0],'k','linewidth', .75);
        set(gca,'xtick',Forcing.rivertime(1):7:Forcing.rivertime(end));
        datetick('x',6,'keeplimits','keepticks');
        ax = axis; ax(3) = min(Forcing.vstress)-.03; ax(4) = max(Forcing.vstress)+.03;
        axis(ax);
        box = [Obs(1).td(1) Obs(1).td(1) Obs(1).td(end) Obs(1).td(end) Obs(1).td(1)];
        plot(box, [ax(3) ax(4) ax(4) ax(3) ax(3)], 'b','linewidth',1.5)
        ylabel('N-S wind stress (N m^{-2})'); xlabel(['Date, ' datestr(Obs(1).td(1),10)]);
%-----%-----%-----%-----%-----%-----%-----%-----%-----%-----%-----%-----%-----%-----

%-----%-----%-----%-----%-----%-----%-----%-----%-----%-----%-----%-----%-----%-----
%----- map figure ------------------------------------------------------------------
filename = [Roms(1).dirname,roms_filename(Roms(1).basename,1)];
[G,S,T] = Z_get_basic_info(filename);
h = nc_varget(filename,'h');
zbot = -h;
zbot(~G.mask_rho)=5;
coast = 'regional';
if(strcmp(Forcing.rivername,'Columbia')); ax = [-126 -123 46 49];
    elseif(strcmp(Forcing.rivername,'Fraser')); ax = [-125 -122.2 47.5 50];
    else ax = [-123.4 -122.1 47 48.5]; 
    coast = 'detailed';
end;    
figure; set(gcf,'position',[465   442   610   503]);
    plot_WAcoast(coast, 'linewidth', 1.5);
    hold on; axis image;
    axis(ax)
    clevs = -[1000 500 200 100];
    [cc,hh] = contour(G.lon_rho,G.lat_rho,zbot,clevs,'k-');
    set(hh(:),'linewidth',1.5,'color',[.5 .5 .5]);
    clabel(cc,hh,'fontsize',10,'color',[.5 .5 .5],'labelspacing',288)
    plot(Forcing.pos(1), Forcing.pos(:,2), 'mp', 'markersize', 8)
    for i=1:n
            plot(Roms(i).pos(1),Roms(i).pos(2),[colors(i,1) 's'])
            hold on; grid on;
            plot(Obs(i).pos(1),Obs(i).pos(2),[colors(i,1) 'o'],'markerfacecolor',colors(i,1))
    end
    
    
%-----%-----%-----%-----%-----%-----%-----%-----%-----%-----%-----%-----%-----%-----    