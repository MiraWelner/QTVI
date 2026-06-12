function ShowDbgPlots(figs)

    for i = 1:length(figs)
        fig = figure(figs{i});
        set(fig, 'DefaultFigureVisible', 'on');
        set(fig, 'Position', get(0, 'Screensize'));
        
        while true
            disp('Press Space to continue');
            try
                w = waitforbuttonpress;
            catch
                return
            end
            switch w
                case 1 % (keyboard press)
                key = get(fig, 'currentcharacter');
                switch key
                    case 32 % space
                        break
                    case 27 % escape key
                        return % break out of the while loop
                    case 13 % return key
                        break
                    otherwise
                        % Wait for a different command.
                end
            end
        end
        close(fig);
    end
end

