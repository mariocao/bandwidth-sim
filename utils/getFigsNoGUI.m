function [ figures ] = getFigsNoGUI( )
%GETFIGSNOGUI Get figures but not GUI
%   Get list of all figures excepting the GUI figure

figures=[];
fig_h = findobj( 0, 'Type', 'Figure' );
for fh = 1:length(fig_h)
    uih = findobj( fig_h(fh), 'Type', 'uicontrol' );
    if isempty( uih )
        figures(end+1) = fig_h(fh);
    end
end

end

