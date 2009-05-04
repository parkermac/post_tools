function compare_roms_obs_ctd_plot(Obs, Roms, Forcing)

n = length(Obs); %number of stations

%-----%-----%-----%-----%-----%-----%-----%-----%-----%-----%-----%-----%-----%-----
%----- plot T/S comparisons -------------------------------------------------------- 
%%% plot
 figure; set(gcf,'position',[177   207   809   731]);
 for i = 1:n
    subplot(n,2,2*i-1)  %plot salt 
        plot(Obs(i).salinity, -Obs(i).depth, 'bo-', 'linewidth', 2)
        hold on; grid on;
        plot(Roms(i).salt, Roms(i).coords.zm, 'kv-', 'linewidth', 2)
        ax = axis; ax(4)=0; axis(ax);
        title(['Sal at : ' num2str(Obs(i).pos(1),'%5.2f') 'W, ' ...
            num2str(Obs(i).pos(2),'%5.2f') 'N, ' datestr(Obs(i).td,1)],'fontsize',16)
        ylabel('Depth (m)','fontsize',16)
        set(gca,'fontsize',14);
    subplot(n,2,2*i)  %plot temp 
        plot(Obs(i).temperature, -Obs(i).depth, 'ro-', 'linewidth', 2)
        hold on; grid on;
        plot(Roms(i).temp, Roms(i).coords.zm, 'kv-', 'linewidth', 2)
        ax2 = axis; axis([ax2(1:2) ax(3:4)]);
        title(['Temp at : ' num2str(Obs(i).pos(1),'%5.2f') 'W, ' ...
            num2str(Obs(i).pos(2),'%5.2f') 'N, ' datestr(Obs(i).td,1)],'fontsize',16)
        ylabel('Depth (m)','fontsize',16)
        set(gca,'fontsize',14);
 end
%-----%-----%-----%-----%-----%-----%-----%-----%-----%-----%-----%-----%-----%-----

%-----%-----%-----%-----%-----%-----%-----%-----%-----%-----%-----%-----%-----%-----
%----- plot forcing ----------------------------------------------------
colors=['b-';'r-';'g-';'m-';'c-';'k-'];
 figure; set(gcf,'position',[214   501   892   412]);
    subplot(211)  %plot river forcing
        plot(Forcing.rivertime, Forcing.riverQ, 'k', 'linewidth', 2)
        hold on; grid on;
        set(gca,'xtick',Forcing.rivertime(1):7:Forcing.rivertime(end),'fontsize',14);
        datetick('x',6,'keeplimits','keepticks');
        ax = axis; axis(ax);
        for i=1:n
            plot([Obs(i).td Obs(i).td], [ax(3) ax(4)], colors(i,:),'linewidth',1.5)
        end
        ylabel('Q (m^3 s^{-1})','fontsize',16); 
        title([Forcing.rivername ' River discharge'],'fontsize',16);
    subplot(212)  %plot N/S wind stress at pos 2 above
        plot(Forcing.time, Forcing.vstress, 'k', 'linewidth', 2.5)
        hold on; grid on;
        plot([Forcing.time(1) Forcing.time(end)],[0 0],'k','linewidth', .75);
        set(gca,'xtick',Forcing.rivertime(1):7:Forcing.rivertime(end),'fontsize',14);
        datetick('x',6,'keeplimits','keepticks');
        ax = axis; ax(3) = min(Forcing.vstress)-.03; ax(4) = max(Forcing.vstress)+.03;
        axis(ax);
        for i=1:n
            plot([Obs(i).td Obs(i).td], [ax(3) ax(4)], colors(i,:),'linewidth',1.5)
        end
        ylabel('N-S wind stress (N m^{-2})','fontsize',16); 
        xlabel(['Date, ' datestr(Obs(1).td,10)],'fontsize',16);
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
    set(gca,'fontsize',14);
    ylabel('Latitude (\circN)','fontsize',16); 
    xlabel('Longitude (\circW)','fontsize',16);
    
%-----%-----%-----%-----%-----%-----%-----%-----%-----%-----%-----%-----%-----%-----    