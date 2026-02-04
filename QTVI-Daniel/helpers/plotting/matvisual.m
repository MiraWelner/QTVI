function matvisual(A,axis)

    % check the input
    if ~isreal(A) || isempty(A) || ischar(A) || ndims(A) > 3
        errordlg('The data array is not suitable for visualization!', ...
                 'Error!', 'modal')
        return
    end

    % determine the matrix size
    [M, N, P] = size(A);

    % loop through the matrix pages
    if ~exist('axis','var')
        figure;
    else
        hold(axis, 'on');
        xlim([0.5,N+.5]);
        ylim([0.5,M+.5]);
    end
    for p = 1:P

        % visualize the matrix page by page
        himg = imagesc(A(:, :, p));
        colormap(flipud(jet));

        % annotation
        if ~exist('axis','var')
            set(gca, 'FontName', 'Times New Roman', 'FontSize', 12)
        else
            set(gca, 'FontName', 'Times New Roman', 'FontSize', 6)
        end
        xlabel('Column number')
        ylabel('Row number')
        if P > 1, title(['Matrix page ' num2str(p)]), end

        set(gca, 'XTick', 1:N)
        set(gca, 'YTick', 1:M)

        hclb = colorbar;
        hclb.Label.String = 'Value';
        hclb.Label.FontName = 'Times New Roman';
        if ~exist('axis','var')
            hclb.Label.FontSize = 12;
        else
            hclb.Label.FontSize = 6;
        end
        
        for m = 1:M
            for n = 1:N
                text(n, m, num2str(A(m, n, p), 3), ...
                    'FontName', 'Times New Roman', ...
                    'FontSize', round(6 + 25./sqrt(M.*N)), ...
                    'HorizontalAlignment', 'center', ...
                    'Rotation', 45)
            end
        end

        % set the datatip UpdateFcn
        cursorMode = datacursormode(gcf);
        set(cursorMode, 'UpdateFcn', {@datatiptxt, himg})

    end

end

function text_to_display = datatiptxt(~, hDatatip, himg)

% determine the current datatip position
pos = get(hDatatip, 'Position');

% form the datatip label
text_to_display = {['Row: ' num2str(pos(2))], ...
                   ['Column: ' num2str(pos(1))], ...
                   ['Value: ' num2str(himg.CData(pos(2), pos(1)))]};

end