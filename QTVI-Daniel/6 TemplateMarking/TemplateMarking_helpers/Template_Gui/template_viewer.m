function varargout = template_viewer(varargin)
    % TEMPLATE_VIEWER MATLAB code for template_viewer.fig
    %      TEMPLATE_VIEWER, by itself, creates a new TEMPLATE_VIEWER or raises the existing
    %      singleton*.
    %
    %      H = TEMPLATE_VIEWER returns the handle to a new TEMPLATE_VIEWER or the handle to
    %      the existing singleton*.
    %
    %      TEMPLATE_VIEWER('CALLBACK',hObject,eventData,handles,...) calls the local
    %      function named CALLBACK in TEMPLATE_VIEWER.M with the given input arguments.
    %
    %      TEMPLATE_VIEWER('Property','Value',...) creates a new TEMPLATE_VIEWER or raises the
    %      existing singleton*.  Starting from the left, property value pairs are
    %      applied to the GUI before template_viewer_OpeningFcn gets called.  An
    %      unrecognized property name or invalid value makes property application
    %      stop.  All inputs are passed to template_viewer_OpeningFcn via varargin.
    %
    %      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
    %      instance to run (singleton)".
    %
    % See also: GUIDE, GUIDATA, GUIHANDLES

    % Edit the above text to modify the response to help template_viewer

    % Last Modified by GUIDE v2.5 11-Mar-2019 17:48:34

    % Begin initialization code - DO NOT EDIT
    gui_Singleton = 1;
    gui_State = struct('gui_Name', mfilename, ...
        'gui_Singleton', gui_Singleton, ...
        'gui_OpeningFcn', @template_viewer_OpeningFcn, ...
        'gui_OutputFcn', @template_viewer_OutputFcn, ...
        'gui_LayoutFcn', [], ...
        'gui_Callback', []);

    if nargin && ischar(varargin{1})
        gui_State.gui_Callback = str2func(varargin{1});
    end

    if nargout
        [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
    else
        gui_mainfcn(gui_State, varargin{:});
    end

end

% End initialization code - DO NOT EDIT

% --- Executes just before template_viewer is made visible.
function template_viewer_OpeningFcn(hObject, eventdata, handles, varargin)
    % This function has no output args, see OutputFcn.
    % hObject    handle to figure
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    % varargin   command line arguments to template_viewer (see VARARGIN)

    % Choose default command line output for template_viewer
    handles.output = hObject;

    templates = varargin{1};
    title_str = varargin{2};
    num_plots = numSubplots(size(templates, 1));
    handles.outputLoc = varargin{3};

    handles.bad_r_templates = false(length(templates),1);
    handles.bad_ppg_templates = false(length(templates),1);
    handles.images = {};
    handles.closearg = 0;
    handles.name = title_str;
    handles.title.String = title_str;
    handles.ppg_samplingrates=cellfun(@(x) x.ppgSamplingRate,templates);

%     handles.ppg_amp_axes = subplot(num_plots(1), num_plots(2), 1);
    handles.markingMode.String = "Mode: Move Subsequent";
    subplot = @(m,n,p,args) subtightplot(m,n,p,[0.06 0.03], [0.01 0.05], [0.01 0.05], args{:});
    handles.dnotches = ones(size(templates,1),1);

    for i = 1:size(templates, 1)
        ax = subplot(num_plots(1), num_plots(2), i, {'ButtonDownFcn', {@plot_ButtonDownFcn, handles}});
        l = template_plot(ax, templates{i}.ppgTemplate, templates{i}.ppgSamplingRate, templates{i}.ecgTemplate,templates{i}.ecgSamplingRate, templates{i}.alignment_point);
        handles.dnotches(i) = l;
        set(ax, 'Tag', ['Bin ' num2str(i)]);

        title(['Bin # ' num2str(templates{i}.index)]);
    end
    updateDnotchs(handles);
%     sgtitle(title_str, 'Interpreter', 'none');
    
%     set(gcf,'WindowState','fullscreen');
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    % Update handles structure
    guidata(hObject, handles);
    
    
    uiwait(handles.figure1);
end

% UIWAIT makes template_viewer wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = template_viewer_OutputFcn(hObject, eventdata, handles)
    % varargout  cell array for returning output args (see VARARGOUT);
    % hObject    handle to figure
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Get default command line output from handles structure

    %     tmp.End = findobj('Tag','End','Type','line');
    %     tmp.Dicrotic = findobj('Tag','Dicrotic','Type','line');
    %     tmp.Peak = findobj('Tag','Peak','Type','line');
    %     tmp.Onset = findobj('Tag','Onset','Type','line');

    
    axes = findobj('-regexp','Tag','Bin');
    outputs = cell(length(axes), 1);

    for i = 1:length(axes)
        ax = axes(i);
        ppg.TemplateBad = handles.bad_ppg_templates(length(axes)-i+1);
        ppg.bad_r_templates = handles.bad_r_templates(length(axes)-i+1);
        ppg.bad_ppg_templates = handles.bad_ppg_templates(length(axes)-i+1);
        if isempty(findobj(ax.Children, 'Tag', 'Dicrotic', 'Type', 'line'))
            ppg.Onset = nan;
            ppg.Peak = nan;
            ppg.Dicrotic = nan;
            ppg.End = nan;
            outputs{i} = ppg;
            continue
        end

        %t = findobj(ax.Children, 'Tag', 'Onset', 'Type', 'line');
        ppg.Onset = nan;
        %t = findobj(ax.Children, 'Tag', 'Peak', 'Type', 'line');
        ppg.Peak = nan;
        t = findobj(ax.Children, 'Tag', 'Dicrotic', 'Type', 'line');
%         ppg.Dicrotic = t.XData(1)*handles.ppg_samplingrates(i);
        ppg.Dicrotic = t.XData(1);
        %t = findobj(ax.Children, 'Tag', 'End', 'Type', 'line');
        ppg.End = nan;

        outputs{i} = ppg;
    end

    varargout{1} = flipud(outputs);
    varargout{2} = handles.closearg;
    if  handles.closearg == 1
        %disp('Saving image of Templates...');

        %export_fig(gcf, fullfile(handles.outputLoc, [handles.name '.png']));
    end

    delete(handles.output);
end

% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, ~, ~)
    % hObject    handle to figure1 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hint: delete(hObject) closes the figure
    if isequal(get(hObject, 'waitstatus'), 'waiting')
        % The GUI is still in UIWAIT, us UIRESUME
        uiresume(hObject);
    else
        delete(hObject);
    end

end

function plot_ButtonDownFcn(hObject, ~, ~)
    handles = guidata(hObject);

    if (strcmp(handles.figure1.SelectionType, 'alt'))
        press = 'r';
    elseif (strcmp(handles.figure1.SelectionType, 'normal'))
        press = 'l';
    else
        press = 'n';
    end

    bin = str2num(hObject.Tag(4:end));
    if press == 'l' && handles.bad_ppg_templates(bin) == 0
        pos = get(hObject, 'CurrentPoint');
        XPOS = pos(1, 1);
        x = inf(4, 1);

        if numel(hObject.Children) == 2

            for i = 1:1
                x(i) = abs(hObject.Children(i).XData(1) - XPOS);
            end

            if ~isempty(x(x < 200))
                [~, idx] = min(x);

                move_all = 0;

                switch hObject.Children(idx).Tag
                    case 'End'

                        if strcmpi(handles.markers_menu.Children(1).Checked, 'on')
                            move_all = 1;
                        end

                    case 'Dicrotic'

                        if strcmpi(handles.markers_menu.Children(2).Checked, 'on')
                            move_all = 1;
                        end

                    case 'Peak'

                        if strcmpi(handles.markers_menu.Children(3).Checked, 'on')
                            move_all = 1;
                        end

                    case 'Onset'

                        if strcmpi(handles.markers_menu.Children(4).Checked, 'on')
                            move_all = 1;
                        end

                end

                if move_all == 1
                    MoveAll(hObject.Children(idx).Tag, hObject.Children(idx).Parent.Tag, handles, XPOS)
                else
                    curr = numel(handles.figure1.Children) + 1 - bin;
                    ax = handles.figure1.Children(curr);
                    if XPOS <= ax.XLim(2) && XPOS >= ax.XLim(1)
                        hObject.Children(idx).XData = [round(XPOS) round(XPOS)];
                    end
                end

                set(ancestor(hObject, 'figure'), 'windowbuttonupfcn', @StopDragging)
                set(ancestor(hObject, 'figure'), 'windowbuttonmotionfcn', {@DragLine, handles, hObject.Children(idx)})
            end

        end
    elseif press == 'n'
        
        if handles.bad_r_templates(bin)
            handles.bad_r_templates(bin) = 0;
            delete(handles.images{bin});
            lines = findobj(hObject.Children,'Type','line','-not','Tag','Dicrotic');
            try
                xmin = 0;
                xmax = -inf;
                ymin = inf;
                ymax = -inf;
                for x = 1:length(lines)
                    xmax = max(max(lines(x).XData),xmax);
                    ymax = max(max(lines(x).YData),ymax);
                    ymin = min(min(lines(x).YData),ymin);
                end
            

                xlim([xmin xmax]);
                ylim([ymin ymax]);
            catch
                xlim('auto');
                ylim('auto');
            end

        else
            warning('off');
            handles.bad_r_templates(bin) = 1;
            %ax4b = get(hObject, 'ButtonDownFcn');    
            hold on;
            i = imread('bad_r.png');
            handles.images{bin}=image(i,'Parent',hObject, 'HitTest','off','PickableParts','none');
            %set(hObject, 'ButtonDownFcn', ax4b{1});
            %set(handles.images{bin}, 'ButtonDownFcn', ax4b{1});
            xlim([1 handles.images{bin}.XData(2)]);
            ylim([1 handles.images{bin}.YData(2)]);

            warning('on');

        end
        guidata(hObject, handles);
    elseif press == 'r'

        if handles.bad_ppg_templates(bin)
            handles.bad_ppg_templates(bin) = 0;
            delete(handles.images{bin});
            lines = findobj(hObject.Children,'Type','line','-not','Tag','Dicrotic');
            try
                xmin = 0;
                xmax = -inf;
                ymin = inf;
                ymax = -inf;
                for x = 1:length(lines)
                    xmax = max(max(lines(x).XData),xmax);
                    ymax = max(max(lines(x).YData),ymax);
                    ymin = min(min(lines(x).YData),ymin);
                end
            

                xlim([xmin xmax]);
                ylim([ymin ymax]);
            catch
                xlim('auto');
                ylim('auto');
            end

        else
            warning('off');
            handles.bad_ppg_templates(bin) = 1;
            %ax4b = get(hObject, 'ButtonDownFcn');    
            hold on;
            i = imread('red_x.png');
            handles.images{bin}=image(i,'Parent',hObject, 'HitTest','off','PickableParts','none');
            %set(hObject, 'ButtonDownFcn', ax4b{1});
            %set(handles.images{bin}, 'ButtonDownFcn', ax4b{1});
            xlim([1 handles.images{bin}.XData(2)]);
            ylim([1 handles.images{bin}.YData(2)]);

            warning('on');

        end
        guidata(hObject, handles);

    end
end

function DragLine(hObject, eventdata, handles, line)
    pos = get(gca, 'CurrentPoint');
    XPOS = pos(1, 1);

    move_all = 0;

    switch line.Tag
        case 'End'

            if strcmpi(handles.markers_menu.Children(1).Checked, 'on')
                move_all = 1;
            end

        case 'Dicrotic'

            if strcmpi(handles.markers_menu.Children(2).Checked, 'on')
                move_all = 1;
            end

        case 'Peak'

            if strcmpi(handles.markers_menu.Children(3).Checked, 'on')
                move_all = 1;
            end

        case 'Onset'

            if strcmpi(handles.markers_menu.Children(4).Checked, 'on')
                move_all = 1;
            end

    end

    if move_all == 1
        MoveAll(line.Tag, line.Parent.Tag, handles, XPOS)
    else
        curr = numel(handles.figure1.Children) + 1 - str2num(line.Parent.Tag);
        ax = handles.figure1.Children(curr);
        if XPOS <= ax.XLim(2) && XPOS >= ax.XLim(1)
            line.XData = [round(XPOS) round(XPOS)];
        end
    end

end

function updateDnotchs(handles)
    warning('off');
    dnotches = handles.dnotches;
    curr = numel(handles.figure1.Children) + 1 - 1;
    outliers = isoutlier(dnotches);
    XPOS = nanmedian(dnotches(~outliers));
    idx = length(outliers);
%     XPOS
%     dnotches
    for i = 5:curr
        if(outliers(idx))
            ax = handles.figure1.Children(i);
            set(ax,'color',[252/255, 254/255, 242/255]);
            if numel(ax.Children) > 1

                if XPOS <= ax.XLim(2) && XPOS >= ax.XLim(1)

                    line = findobj(ax.Children,'Tag','Dicrotic');
                    if ~isempty(line)
                        line.XData = [round(XPOS) round(XPOS)];
                    end

                end

            end
        end
            
        idx = idx-1;

    end

    warning('on');
end


function MoveAll(type, num, handles, XPOS)
    warning('off');

    curr = numel(handles.figure1.Children) + 1 - str2num(num);

    for i = 5:curr
        idx = 0;

%         switch type
%             case 'End'
%                 idx = 1;
%             case 'Dicrotic'
%                 idx = 2;
%             case 'Peak'
%                 idx = 3;
%             case 'Onset'
%                 idx = 4;
%         end

        ax = handles.figure1.Children(i);

        if numel(ax.Children) > 1

            if XPOS <= ax.XLim(2) && XPOS >= ax.XLim(1)
                switch type
                    case 'End'
                        idx = 1;
                    case 'Dicrotic'
                        line = findobj(ax.Children,'Tag','Dicrotic');
                        if ~isempty(line)
                            line.XData = [round(XPOS) round(XPOS)];
                        end
                    case 'Peak'
                        idx = 3;
                    case 'Onset'
                        idx = 4;
                end
            end

        end

    end

    warning('on');
end

function StopDragging(hObject, ~, ~)
    set(hObject, 'windowbuttonmotionfcn', '');
    set(hObject, 'windowbuttonupfcn', '');
end

% --------------------------------------------------------------------
function markers_menu_Callback(hObject, eventdata, handles)
    % hObject    handle to markers_menu (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
end

% --------------------------------------------------------------------
function PPG_Onset_bool_Callback(hObject, eventdata, handles)
    % hObject    handle to PPG_Onset_bool (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    if strcmpi(hObject.Checked, 'on')
        handles.PPG_Onset_bool = 0;
        hObject.Checked = 'off';
    else
        handles.PPG_Onset_bool = 1;
        hObject.Checked = 'on';
    end

    guidata(hObject, handles);
end

% --------------------------------------------------------------------
function PPG_peak_bool_Callback(hObject, eventdata, handles)
    % hObject    handle to PPG_peak_bool (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    if strcmpi(hObject.Checked, 'on')
        handles.PPG_peak_bool = 0;
        hObject.Checked = 'off';
    else
        handles.PPG_peak_bool = 1;
        hObject.Checked = 'on';
    end

    guidata(hObject, handles);
end

% --------------------------------------------------------------------
function PPG_Dicrotic_Notch_bool_Callback(hObject, eventdata, handles)
    % hObject    handle to PPG_Dicrotic_Notch_bool (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    if strcmpi(hObject.Checked, 'on')
        handles.PPG_Dicrotic_Notch_bool = 0;
        hObject.Checked = 'off';
        handles.markingMode.String = "Mode: Move Individual";
    else
        handles.PPG_Dicrotic_Notch_bool = 1;
        hObject.Checked = 'on';
        handles.markingMode.String = "Mode: Move Subsequent";
    end

    guidata(hObject, handles);
end

% --------------------------------------------------------------------
function PPG_end_bool_Callback(hObject, eventdata, handles)
    % hObject    handle to PPG_end_bool (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    if strcmpi(hObject.Checked, 'on')
        handles.PPG_end_bool = 0;
        hObject.Checked = 'off';
    else
        handles.PPG_end_bool = 1;
        hObject.Checked = 'on';
    end

    guidata(hObject, handles);
end

% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over finish_button.
function finish_button_ButtonDownFcn(hObject, eventdata, handles)
    % hObject    handle to finalize_btn (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
%     resp = questdlg("Are you sure?", "Finalize Selections");
% 
%     if (strcmp(resp, "Yes"))
%         handles.closearg = 1;
%         guidata(hObject, handles);
% 
%         close(handles.output);
%     end

    handles.closearg = 1;
    guidata(hObject, handles);

    close(handles.output);
end
