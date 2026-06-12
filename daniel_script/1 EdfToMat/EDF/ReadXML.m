function [length_sec, sleepStages, scoredEvents] = ReadXML(path)
    xDoc = xmlread(path);
    
    length_sec = str2double(xDoc.getElementsByTagName('EpochLength').item(0).getTextContent); 
    
    xmlSleepStages = xDoc.getElementsByTagName('SleepStage');
    sleepStages = nan(xmlSleepStages.getLength,1);
    for k = 0:xmlSleepStages.getLength-1
        sleepStages(k+1) = str2double(xmlSleepStages.item(k).item(0).getTextContent);
    end
    
    xmlScoredEvents = xDoc.getElementsByTagName('ScoredEvent');
    scoredEvents = cell(xmlScoredEvents.getLength,1);
    %tmp = [];
    for k = 0:xmlScoredEvents.getLength-1
        event.Name = xmlScoredEvents.item(k).getElementsByTagName('Name').item(0).getTextContent;
        %tmp = [tmp;event.Name];

        event.Start_sec = str2double(xmlScoredEvents.item(k).getElementsByTagName('Start').item(0).getTextContent);
        event.Duration_sec = str2double(xmlScoredEvents.item(k).getElementsByTagName('Duration').item(0).getTextContent);
        scoredEvents{k+1} = event;
    end
end