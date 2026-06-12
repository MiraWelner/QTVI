function varargout = ReviewData(varargin)
    % REVIEWDATA MATLAB code for ReviewData.fig
    %      REVIEWDATA, by itself, creates a new REVIEWDATA or raises the existing
    %      singleton*.
    %
    %      H = REVIEWDATA returns the handle to a new REVIEWDATA or the handle to
    %      the existing singleton*.
    %
    %      REVIEWDATA('CALLBACK',hObject,eventData,handles,...) calls the local
    %      function named CALLBACK in REVIEWDATA.M with the given input arguments.
    %
    %      REVIEWDATA('Property','Value',...) creates a new REVIEWDATA or raises the
    %      existing singleton*.  Starting from the left, property value pairs are
    %      applied to the GUI before ReviewData_OpeningFcn gets called.  An
    %      unrecognized property name or invalid value makes property application
    %      stop.  All inputs are passed to ReviewData_OpeningFcn via varargin.
    %
    %      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
    %      instance to run (singleton)".
    %
    % See also: GUIDE, GUIDATA, GUIHANDLES

    % Edit the above text to modify the response to help ReviewData

    % Last Modified by GUIDE v2.5 28-Dec-2020 14:33:29

    % Begin initialization code - DO NOT EDIT
    gui_Singleton = 1;
    gui_State = struct('gui_Name', mfilename, ...
        'gui_Singleton', gui_Singleton, ...
        'gui_OpeningFcn', @ReviewData_OpeningFcn, ...
        'gui_OutputFcn', @ReviewData_OutputFcn, ...
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

% --- Executes just before ReviewData is made visible.
function ReviewData_OpeningFcn(hObject, ~, handles, varargin)
    % This function has no output args, see OutputFcn.
    % hObject    handle to figure
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    % varargin   command line arguments to ReviewData (see VARARGIN)

    % Choose default command line output for ReviewData
    handles.output = hObject;

    %  setting up variables
    handles.ColorScheme = [173, 73, 255;
                            72, 101, 232;
                            75, 202, 255;
                            73, 232, 130;
                            79, 255, 113];

    handles.ColorScheme = (handles.ColorScheme - 25) / 256;
    handles.ctrl = 0;
    
    inputs = varargin{1};

    handles.edf_label.String = join(["Currently processing:" strcat(string(inputs.name), ".edf")]);
    handles.ModeText.String = "Mode: Mark Individual";

    handles.prev_key_press = datetime('now');
    
    handles.reviewing = inputs.reviewing;

    handles.ppg = inputs.ppg;
    handles.ppgTime_sec = inputs.ppgTime_sec;
    handles.ppgSamplingRate = inputs.ppgSamplingRate;
    handles.ppgAmpogram = inputs.ppgAmpogram;
    handles.ppgAmpogram_timeIndex = inputs.ppgAmpogram_timeIndex;
    handles.ppg_std = std(handles.ppg);
    handles.ppg_median = median(handles.ppg);
    handles.ppg_mean = mean(handles.ppg);

    handles.ecg = inputs.ecg;
    handles.ecgTime_sec = inputs.ecgTime_sec;
    handles.ecgSampleingRate = inputs.ecgSamplingRate;
    handles.ecgAmpogram = inputs.ecgAmpogram;
    handles.ecgAmpogram_timeIndex = inputs.ecgAmpogram_timeIndex;
    handles.ecg_std = std(handles.ecg);
    handles.ecg_median = median(handles.ecg);
    handles.ecg_mean = mean(handles.ecg);
    % [~, ~, R_amp, ~, ~, ~] = rpeakdetect(inputs.ecg', inputs.ecgSamplingRate);
    R_amp = handles.ecg(inputs.rIndex);
    handles.ecg_r = R_amp;
    handles.ecg_peak_mean = mean(R_amp);
    handles.ecg_peak_std = std(R_amp);
    handles.ecg_peak_median = median(R_amp);

    handles.sleepStageTime_hrs = inputs.sleepStageTime_hrs;
    handles.sleepStages = inputs.sleepStages;
    handles.scoredEvents = inputs.scoredEvents;
    handles.sleepStageSampleRate = inputs.sleepStageSampleRate;

    max_length = max([numel(inputs.ppgTime_sec) max([numel(inputs.ecgTime_sec) numel(inputs.sleepStageTime_hrs)])]);
    max_samplerate = max([inputs.ecgSamplingRate inputs.ppgSamplingRate inputs.sleepStageSampleRate]);
    handles.time_sec = (0:max_length - 1) / max_samplerate;
    handles.sleepStages_timeIndex = 1:(1 / handles.sleepStageSampleRate * max_samplerate):length(handles.time_sec);
    handles.sleepStages_timeIndex = handles.sleepStages_timeIndex(1:length(handles.sleepStages));

    handles.closeResponse = 0;
    handles.viewWidth = 30;
    handles.CurrentPosition = 30; % seconds
    handles.startSeg = 0;
    handles.endSeg = 30;
    handles.GenExc.ind = [];
    handles.GenExc.noiseExc = [];
    handles.GenExc.ind = [];
    handles.GenExc.last = [];
    handles.GenExc.noisecnt = 0;
    
    
    
    if isfield(inputs,'noiseInfo') == 1
        
        for x = 1:size(inputs.noiseInfo.noise_markings, 1)
            handles.GenExc.noiseExc(end + 1, 1:2) = [inputs.noiseInfo.noise_markings(x,1), inputs.noiseInfo.noise_markings(x,2)];
            handles.GenExc.ind(end + 1, 1:2) = [inputs.noiseInfo.noise_markings(x,1), inputs.noiseInfo.noise_markings(x,2)];
            handles.GenExc.last(end + 1) = inputs.noiseInfo.noise_markings(x,5);
            handles.GenExc.noisecnt = handles.GenExc.noisecnt + 1;
        end

    end
%     handles.GenExc.noiseExc(end + 1, 1:2) = [min([handles.GenExc.Lbegin, Lend]), max([handles.GenExc.Lbegin, Lend])];
%     handles.GenExc.ind(end + 1, 1:2) = [min([handles.GenExc.Lbegin, Lend]), max([handles.GenExc.Lbegin, Lend])];
%     handles.GenExc.last(end + 1) = press;
%     handles.GenExc.Lbegin = 0;
%     handles.GenExc.noisecnt = handles.GenExc.noisecnt + 1;
    
    handles.GenExc.Lbegin = 0;
    handles.GenExc.Rbegin = 0;
    handles.setup = 0;
    handles.scrollMod = .5;
    handles.data_y_scale = 1;
    handles = updatePlots(handles);
    handles = updateMetadataPlots(handles, 0);
    if handles.GenExc.noisecnt > 0
        handles = updateMetadataPlots(handles, 1);
    end
    handles = updateCursor(handles);
    handles = updateScroll(handles);

    set(handles.sleep_state_axes, 'YTick', []);
    set(handles.events_axes, 'YTick', []);
    set(handles.ekg_amp_axes, 'YTick', []);
    set(handles.ppg_amp_axes, 'YTick', []);

    set(handles.ppg_axes, 'YColor', 'k');

    grid(handles.sleep_state_axes, 'on');
    grid(handles.events_axes, 'on');
    grid(handles.ekg_amp_axes, 'on');
    grid(handles.ppg_amp_axes, 'on');
    grid(handles.ekg_axes, 'on');

    handles.sleep_state_axes.YGrid = 'off';
    handles.sleep_state_axes.XGrid = 'on';
    handles.events_axes.YGrid = 'off';
    handles.events_axes.XGrid = 'on';
    handles.ekg_amp_axes.YGrid = 'off';
    handles.ekg_amp_axes.XGrid = 'on';
    handles.ppg_amp_axes.YGrid = 'off';
    handles.ppg_amp_axes.XGrid = 'on';
    handles.ekg_axes.YGrid = 'off';
    handles.ekg_axes.XGrid = 'on';

    set(handles.sleep_state_axes, 'xticklabel', [])
    set(handles.events_axes, 'xticklabel', [])
    set(handles.ekg_amp_axes, 'xticklabel', [])
    set(handles.ekg_axes, 'xticklabel', [])

    % Update handles structure
    guidata(hObject, handles);
    uiwait(handles.figure1);

end

% UIWAIT makes ReviewData wait for user response (see UIRESUME)
% --- Outputs from this function are returned to the command line.

function varargout = ReviewData_OutputFcn(~, ~, handles)
    % varargout  cell array for returning output args (see VARARGOUT);
    % hObject    handle to figure
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Get default command line output from handles structure
    varargout{1} = handles;
    delete(handles.output);
end

% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(hObject, eventdata, handles)
    % hObject    handle to figure1 (see GCBO)
    % eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
    %     Key: name of the key that was pressed, in lower case
    %     Character: character interpretation of the key(s) that was pressed
    %     Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
    % handles    structure with handles and user data (see GUIDATA)
    %disp("");
    keypress = eventdata.Key;
    
    if (strcmp(keypress, 'rightarrow') || strcmp(keypress, 'leftarrow')) || (strcmp(keypress, 'a') || strcmp(keypress, 'd'))
        if (strcmp(keypress, 'rightarrow')) || strcmp(keypress, 'd')
            handles.CurrentPosition = min(floor(handles.ppgTime_sec(end) - handles.viewWidth / 2), handles.CurrentPosition + handles.viewWidth * handles.scrollMod);
            handles.startSeg = handles.CurrentPosition - handles.viewWidth / 2;
            handles.endSeg = handles.CurrentPosition + handles.viewWidth / 2;
        else
            handles.CurrentPosition = max(0 + handles.viewWidth / 2, handles.CurrentPosition - handles.viewWidth * handles.scrollMod);
            handles.startSeg = handles.CurrentPosition - handles.viewWidth / 2;
            handles.endSeg = handles.CurrentPosition + handles.viewWidth / 2;
        end
    elseif strcmp(keypress, 'space')
        if handles.ctrl == 1
            handles.ctrl = 0;
            handles.ModeText.String = "Mode: Mark Individual";
        else
            handles.ctrl = 1;
            handles.ModeText.String = "Mode: Mark Both";
        end

    end

    handles = updatePlots(handles);
    handles = updateCursor(handles);
    guidata(hObject, handles);
end

% --- Executes on key release with focus on figure1 and none of its controls.
function figure1_KeyReleaseFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
%	Key: name of the key that was released, in lower case
%	Character: character interpretation of the key(s) that was released
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) released
% handles    structure with handles and user data (see GUIDATA)
%     keypress = eventdata.Key;
% %     disp('release');
% %     disp(keypress);
%     if strcmp(keypress, 'space')
%         handles.ctrl = 0;
%         handles.ModeText.String = "Mode: Mark Individual";
%     end
    guidata(hObject, handles);
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

% --- Executes on button press in undo_btn.
function undo_btn_Callback(hObject, ~, handles)
    % hObject    handle to undo_btn (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    handles.GenExc.ind(end, :) = [];

    handles.GenExc.noisecnt = handles.GenExc.noisecnt - 1;
    handles.GenExc.noiseExc(end, :) = [];

    handles.GenExc.last(end) = [];
    handles = updateMetadataPlots(handles, 1);
    handles = updatePlots(handles);
    handles = updateCursor(handles);
    guidata(hObject, handles);
    set(handles.undo_btn, 'Enable', 'off');
    drawnow;
    set(handles.undo_btn, 'Enable', 'on');
end

% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over undo_btn.
function undo_btn_ButtonDownFcn(~, ~, ~)
    % hObject    handle to undo_btn (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
end

% --- Executes on button press in clear_btn.
function clear_btn_Callback(hObject, ~, handles)
    % hObject    handle to clear_btn (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    resp = questdlg("Are you sure? this will clear all data selections made", "Clear All Selections", 'No');

    if (strcmp(resp, "Yes"))
        handles.GenExc.ind(:, :) = [];
        handles.GenExc.noisecnt = 0;
        handles.GenExc.noiseExc(:, :) = [];
        handles.GenExc.last(:) = [];
        handles = updateMetadataPlots(handles, 0);
        handles = updatePlots(handles);
        handles = updateCursor(handles);
        guidata(hObject, handles);
    end

end

% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over clear_btn.
function clear_btn_ButtonDownFcn(~, ~, ~)
    % hObject    handle to clear_btn (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
end

% --- Executes on button press in finalize_btn.
function finalize_btn_Callback(hObject, ~, handles)
    % hObject    handle to finalize_btn (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    resp = questdlg("Are you sure?", "Finalize Selections");

    if (strcmp(resp, "Yes"))
        handles.closeResponse = 1;
        guidata(hObject, handles);

        close(handles.output);
    end

end

% --- Executes on mouse press over axes background.
function ekg_axes_ButtonDownFcn(hObject, ~, handles)
    % hObject    handle to ppg_axes (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    ppg_axes_ButtonDownFcn(hObject, 1, handles);
end

% --- Executes on mouse press over axes background.
function ppg_axes_ButtonDownFcn(hObject, ~, handles)
    % hObject    handle to ppg_axes (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    pos = get(hObject, 'CurrentPoint');
    XPOS = pos(1, 1);

    if (strcmp(handles.figure1.SelectionType, 'alt'))
        press = 'r';
    elseif (strcmp(handles.figure1.SelectionType, 'normal'))
        press = 'l';
    else
        press = 'n';
    end

    if ((press == 'l' || press == 'r') && handles.ctrl == 1)
        if ((press == 'l' || press == 'l') && handles.GenExc.Rbegin == 0)

            if (handles.GenExc.Lbegin == 0)
                handles.GenExc.Lbegin = XPOS;
                handles.GenExc.vline = vline(seconds(XPOS), 'r');

            else
                Lend = XPOS;
                handles.GenExc.noiseExc(end + 1, 1:2) = [min([handles.GenExc.Lbegin, Lend]), max([handles.GenExc.Lbegin, Lend])];
                handles.GenExc.ind(end + 1, 1:2) = [min([handles.GenExc.Lbegin, Lend]), max([handles.GenExc.Lbegin, Lend])];
                handles.GenExc.last(end + 1) = 'c';
                handles.GenExc.Lbegin = 0;
                handles.GenExc.noisecnt = handles.GenExc.noisecnt + 1;
                handles = updateMetadataPlots(handles, 2);
                handles = updatePlots(handles);
                delete(handles.GenExc.vline)
            end
        end
        guidata(hObject, handles);
    elseif (press == 'l' || press == 'r')

        if (press == 'l' && handles.GenExc.Rbegin == 0)

            if (handles.GenExc.Lbegin == 0)
                handles.GenExc.Lbegin = XPOS;
                handles.GenExc.vline = vline(seconds(XPOS), 'b');

            else
                Lend = XPOS;
                handles.GenExc.noiseExc(end + 1, 1:2) = [min([handles.GenExc.Lbegin, Lend]), max([handles.GenExc.Lbegin, Lend])];
                handles.GenExc.ind(end + 1, 1:2) = [min([handles.GenExc.Lbegin, Lend]), max([handles.GenExc.Lbegin, Lend])];
                handles.GenExc.last(end + 1) = press;
                handles.GenExc.Lbegin = 0;
                handles.GenExc.noisecnt = handles.GenExc.noisecnt + 1;
                handles = updateMetadataPlots(handles, 2);
                handles = updatePlots(handles);
                delete(handles.GenExc.vline)
            end

        elseif (handles.GenExc.Lbegin == 0)

            if (handles.GenExc.Rbegin == 0)
                handles.GenExc.Rbegin = XPOS;
                handles.GenExc.vline = vline(seconds(XPOS), 'g');

            else
                Rend = XPOS;
                handles.GenExc.noiseExc(end + 1, 1:2) = [min([handles.GenExc.Rbegin, Rend]), max([handles.GenExc.Rbegin, Rend])];
                handles.GenExc.ind(end + 1, 1:2) = [min([handles.GenExc.Rbegin, Rend]), max([handles.GenExc.Rbegin, Rend])];
                handles.GenExc.last(end + 1) = press;
                handles.GenExc.Rbegin = 0;
                handles.GenExc.noisecnt = handles.GenExc.noisecnt + 1;
                handles = updateMetadataPlots(handles, 2);
                handles = updatePlots(handles);
                delete(handles.GenExc.vline)

            end

        end

        guidata(hObject, handles);
    end

end

% --- Executes on mouse press over axes background.
function ppg_amp_axes_ButtonDownFcn(hObject, ~, handles)
    pos = get(hObject, 'CurrentPoint');
    XPOS = pos(1, 1);
    handles.CurrentPosition = round(XPOS * 60 * 60);
    handles.startSeg = handles.CurrentPosition - handles.viewWidth / 2;
    handles.endSeg = handles.CurrentPosition + handles.viewWidth / 2;
    handles = updatePlots(handles);
    handles = updateCursor(handles);
    guidata(hObject, handles);
end

% --- Executes on mouse press over axes background.
function ekg_amp_axes_ButtonDownFcn(hObject, ~, handles)
    pos = get(hObject, 'CurrentPoint');
    XPOS = pos(1, 1);
    handles.CurrentPosition = round(XPOS * 60 * 60);
    handles.startSeg = handles.CurrentPosition - handles.viewWidth / 2;
    handles.endSeg = handles.CurrentPosition + handles.viewWidth / 2;
    handles = updatePlots(handles);
    handles = updateCursor(handles);
    guidata(hObject, handles);
end

% --- Executes during object creation, after setting all properties.
function window_time_selection_CreateFcn(hObject, ~, ~)
    % hObject    handle to window_time_selection (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    empty - handles not created until after all CreateFcns called

    % Hint: listbox controls usually have a white background on Windows.
    %       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject, 'BackgroundColor'), get(0, 'defaultUicontrolBackgroundColor'))
        set(hObject, 'BackgroundColor', 'white');
    end

end

% --- Executes on selection change in window_time_selection.
function window_time_selection_Callback(hObject, ~, handles)
    % hObject    handle to window_time_selection (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hints: contents = cellstr(get(hObject,'String')) returns window_time_selection contents as cell array
    %        contents{get(hObject,'Value')} returns selected item from window_time_selection
    handles = updateWidth(handles);
    handles = updatePlots(handles);
    handles = updateCursor(handles);
    guidata(hObject, handles);
end

% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over window_time_selection.
function window_time_selection_ButtonDownFcn(hObject, ~, handles)
    % hObject    handle to window_time_selection (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    handles = updateWidth(handles);
    handles = updatePlots(handles);
    handles = updateCursor(handles);
    guidata(hObject, handles);
end

% --- Executes on slider movement.
function scroll_length_slider_Callback(hObject, ~, handles)
    % hObject    handle to scroll_length_slider (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    handles = updateScroll(handles);
    guidata(hObject, handles);
end

% --- Executes during object creation, after setting all properties.
function scroll_length_slider_CreateFcn(hObject, ~, ~)
    % hObject    handle to scroll_length_slider (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    empty - handles not created until after all CreateFcns called

    % Hint: slider controls usually have a light gray background.
    if isequal(get(hObject, 'BackgroundColor'), get(0, 'defaultUicontrolBackgroundColor'))
        set(hObject, 'BackgroundColor', [.9 .9 .9]);
    end

end

% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over scroll_length_slider.
function scroll_length_slider_ButtonDownFcn(hObject, ~, handles)
    % hObject    handle to scroll_length_slider (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    handles = updateScroll(handles);
    guidata(hObject, handles);
end

% --- Executes during object creation, after setting all properties.
function scroll_length_txt_CreateFcn(hObject, ~, ~)
    % hObject    handle to scroll_length_txt (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    empty - handles not created until after all CreateFcns called

    % Hint: edit controls usually have a white background on Windows.
    %       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject, 'BackgroundColor'), get(0, 'defaultUicontrolBackgroundColor'))
        set(hObject, 'BackgroundColor', 'white');
    end

end

function scroll_length_txt_Callback(~, ~, ~)
    % hObject    handle to scroll_length_txt (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hints: get(hObject,'String') returns contents of scroll_length_txt as text
    %        str2double(get(hObject,'String')) returns contents of scroll_length_txt as a double
end

% --------------------------------------------------------------------
function Options_menu_Callback(hObject, eventdata, handles)
    % hObject    handle to Options_menu (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
end

% --------------------------------------------------------------------
function scale_window_to_data_min_max_option_Callback(hObject, eventdata, handles)
    % hObject    handle to scale_window_to_data_min_max_option (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    handles.fix_window_to_data_option.Checked = 'off';
    handles.scale_window_to_data_min_max_option.Checked = 'on';
    handles.data_y_scale = 1;
    handles = updatePlots(handles);
    guidata(hObject, handles);
end

% --------------------------------------------------------------------
function fix_window_to_data_option_Callback(hObject, eventdata, handles)
    % hObject    handle to fix_window_to_data_option (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    handles.fix_window_to_data_option.Checked = 'on';
    handles.scale_window_to_data_min_max_option.Checked = 'off';
    handles.data_y_scale = 0;
    handles = updatePlots(handles);
    guidata(hObject, handles);
end



% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
activity_legend();
end



% --- Executes on mouse press over axes background.
function events_axes_ButtonDownFcn(hObject, eventdata, handles)
    pos = get(hObject, 'CurrentPoint');
    XPOS = pos(1, 1);
    handles.CurrentPosition = round(XPOS * 60 * 60);
    handles.startSeg = handles.CurrentPosition - handles.viewWidth / 2;
    handles.endSeg = handles.CurrentPosition + handles.viewWidth / 2;
    handles = updatePlots(handles);
    handles = updateCursor(handles);
    guidata(hObject, handles);
end
function clear_btn_KeyPressFcn(hObject, eventdata, handles)

end
% --- Executes on mouse press over axes background.
function sleep_state_axes_ButtonDownFcn(hObject, eventdata, handles)
    pos = get(hObject, 'CurrentPoint');
    XPOS = pos(1, 1);
    handles.CurrentPosition = round(XPOS * 60 * 60);
    handles.startSeg = handles.CurrentPosition - handles.viewWidth / 2;
    handles.endSeg = handles.CurrentPosition + handles.viewWidth / 2;
    handles = updatePlots(handles);
    handles = updateCursor(handles);
    guidata(hObject, handles);
end
