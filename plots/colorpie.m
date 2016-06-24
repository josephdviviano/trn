%% color pie

%%%

interval = 0.1;
increment = 1800;
colorflip = 0;
printtitle = 0;

%%%

figure('Color',[1 1 1]);

x = ones(1,1800);
h = pie(x);
colormap(hsv(3600));

for hc=1:2:length(h) % segments
    set(h(hc),'EdgeColor','none');
end

% texts
tickx = ['\pi', '2\pi'];

for hc=2:2:length(h) 
    if mod(hc,increment) ~= 0;
        % delete some texts
        delete(h(hc)); 
    else
        % set new texts
        set(h(hc),'string',num2str(hc/10));
        %hc = format_tick(h, tickx);
    end
end

% label chart
if printtitle == 1; title('Rotation in Degrees'); end;

% make all text in the figure to size 14 and bold
figureHandle = gcf;
set(findall(figureHandle,'type','text'),'fontSize',24,'fontWeight','bold');

% reverse order of pie chart (and colormap?)
set(gca,'XDir','reverse');
if colorflip == 1; colormap(flipud(colormap)); end

set(gcf, 'Position', [0 0 1024 768]);
%export_fig(['fig_color_pie.tiff'], '-tif', '-m2');


%%

% syms x
% f = taylor(log(1+x));
% ezplot(f)
% hold on
% title(['$' latex(f) '$'],'interpreter','latex')
% hold off