function handles = updateScroll(handles)
    scollMod = get(handles.scroll_length_slider, 'Value');
    reallength = scollMod * handles.viewWidth;
    set(handles.scroll_length_txt, 'String', num2str(reallength));
    handles.scrollMod = scollMod;
end
