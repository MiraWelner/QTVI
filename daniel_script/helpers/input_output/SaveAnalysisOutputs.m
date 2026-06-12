function SaveAnalysisOutputs(output, filename)
    output = table(output.RR_interval(:), output.RtoBP(:), output.RtoSyst(:), output.MPS_PreSystolic_d(:), output.MNS_PostSystolic_d(:), ...
        output.MNS_PostDicrotic_d(:), output.DicroticNotch_d(:), output.DicroticPeak_d(:), output.SystolicPeak_d(:), output.Foot_d(:), ...
        output.MPSlope(:), output.BPmaxAMP(:), output.PulseIntegral(:), ...
        'VariableNames', {'RR_interval''RtoBP''RtoSyst''MPS_PreSystolic_d'...
        'MNS_PostSystolic_d''MNS_PostDicrotic_d''DicroticNotch_d''DicroticPeak_d''SystolicPeak_d''Foot_d''MPSlope''BPmaxAMP''PulseIntegral'});

    filename = strcat(filename, '_PulseAnalysis_');
    filename = strcat(filename, date);
    filename = strcat(filename, '.csv');
    writetable(output, filename);

    disp('Saved As:');
    disp(filename);
end
