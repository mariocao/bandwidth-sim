function [  ] = saveFigures( filePath, figures, gui_on )
%SAVEFIGURES Saves displayed figures
%   Saves displayed figures in different file formats

GUI_ON = 0;
if (nargin < 2 || nargin > 3 )
    error('usage: saveFigures( filePath, hfigs, [gui_on] )')
elseif nargin == 3
    GUI_ON=gui_on;
end

% Save list of figures
if(GUI_ON)
    h = waitbar(0,'Saving figures...','Name','Bandwitdh Allocation Simulation');
    total_figs = length(figures);
    for m = 1:length(figures)
        waitbar(m/total_figs,h);
        if (ishandle(figures(m)))
            %figure(hfigs(m))                                %Bring Figure to foreground
            %filename = input('Filename? (0 to skip)\n', 's')%Prompt user
            filename=strcat(filePath,'\\figure',num2str(figures(m)));
            saveas(figures(m), [filename '.fig']); %Matlab .FIG file
            saveas(figures(m), [filename '.emf']); %Windows Enhanced Meta-File (best for powerpoints)
            saveas(figures(m), [filename '.png']); %Standard PNG graphics file (best for web)
            saveas(figures(m),filename,'epsc2');
            %eval(['print -depsc2 ' filename])   %Enhanced Postscript (Level 2 color) (Best for LaTeX documents)
        end
    end
    waitbar(1,h);
    delete(h);
else %NO GUI
    for m = 1:length(figures)
        if (ishandle(figures(m)))
            %figure(hfigs(m))                                %Bring Figure to foreground
            %filename = input('Filename? (0 to skip)\n', 's')%Prompt user
            filename=strcat(filePath,'\\figure',num2str(figures(m)));
            saveas(figures(m), [filename '.fig']); %Matlab .FIG file
            saveas(figures(m), [filename '.emf']); %Windows Enhanced Meta-File (best for powerpoints)
            saveas(figures(m), [filename '.png']); %Standard PNG graphics file (best for web)
            saveas(figures(m),filename,'epsc2');
            %eval(['print -depsc2 ' filename])   %Enhanced Postscript (Level 2 color) (Best for LaTeX documents)
        end
    end
end

end

