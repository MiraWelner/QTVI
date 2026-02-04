clear;
close all;


% dictionary contains the name in beats_flattened as key and the output
% name as value
% dict = containers.Map({'abs_amp_foot' 'abs_amp_peak' 'adjusted_sleep_state' 'amp_baselined_diastolic_peaks' 'amp_baselined_dicrotic_notches' 'amp_baselined_feets' 'amp_baselined_neg_slopes_after_dnotch' 'amp_baselined_neg_slopes_pre_dnotch' 'amp_baselined_pos_slopes' 'amp_baselined_systolic_peaks' 'amp_baselined_tP_20' 'amp_baselined_tP_20_inv' 'amp_baselined_tP_50' 'amp_baselined_tP_50_inv' 'amp_baselined_tP_80' 'amp_baselined_tP_80_inv' 'amp_baselined_tR_20' 'amp_baselined_tR_20_inv' 'amp_baselined_tR_50' 'amp_baselined_tR_50_inv' 'amp_baselined_tR_80' 'amp_baselined_tR_80_inv' 'amp_delta_systolic' 'amp_raw_diastolic_peaks' 'amp_raw_dicrotic_notches' 'amp_raw_feets' 'amp_raw_neg_slopes_after_dnotch' 'amp_raw_neg_slopes_pre_dnotch' 'amp_raw_pos_slopes' 'amp_raw_systolic_peaks' 'amp_raw_tP_20' 'amp_raw_tP_20_inv' 'amp_raw_tP_50' 'amp_raw_tP_50_inv' 'amp_raw_tP_80' 'amp_raw_tP_80_inv' 'amp_raw_tR_20' 'amp_raw_tR_20_inv' 'amp_raw_tR_50' 'amp_raw_tR_50_inv' 'amp_raw_tR_80' 'amp_raw_tR_80_inv' 'amp_raw_vallies' 'area' 'area_baselined' 'correct_idx_begin' 'corrected_time_sec' 'edge_beat_mask' 'error_ppg_segmentation' 'idx_begin' 'idx_diastolic' 'idx_dnotch' 'idx_end' 'idx_foot' 'idx_neg_slope_after' 'idx_neg_slope_b4' 'idx_pos_slope' 'idx_systolic' 'msec_beat_length' 'msec_R_2_diastolic_peak' 'msec_R_2_dicrotic_notch' 'msec_R_2_first_valley' 'msec_R_2_foot' 'msec_R_2_negslopes_post_dnotch' 'msec_R_2_negslopes_pre_dnotch' 'msec_R_2_pos_slope' 'msec_R_2_second_valley' 'msec_R_2_systolic_peak' 'msec_R_2_tP_20' 'msec_R_2_tP_20_inv' 'msec_R_2_tP_50' 'msec_R_2_tP_50_inv' 'msec_R_2_tP_80' 'msec_R_2_tP_80_inv' 'msec_R_2_tR_20' 'msec_R_2_tR_20_inv' 'msec_R_2_tR_50' 'msec_R_2_tR_50_inv' 'msec_R_2_tR_80' 'msec_R_2_tR_80_inv' 'msec_total_duration_20' 'msec_total_duration_50' 'msec_total_duration_80' 'msec_total_duration_tR_20' 'msec_total_duration_tR_50' 'msec_total_duration_tR_80' 'msec_tP_50_2_diastolic_peak' 'msec_tP_50_2_dicrotic_notch' 'msec_tP_50_2_first_valley' 'msec_tP_50_2_foot' 'msec_tP_50_2_negslopes_post_dnotch' 'msec_tP_50_2_negslopes_pre_dnotch' 'msec_tP_50_2_pos_slope' 'msec_tP_50_2_second_valley' 'msec_tP_50_2_systolic_peak' 'msec_tP_50_2_tP_20' 'msec_tP_50_2_tP_20_inv' 'msec_tP_50_2_tP_80' 'msec_tP_50_2_tP_80_inv' 'msec_tP_50_2_tR_20' 'msec_tP_50_2_tR_20_inv' 'msec_tP_50_2_tR_50' 'msec_tP_50_2_tR_50_inv' 'msec_tP_50_2_tR_80' 'msec_tP_50_2_tR_80_inv' 'ppg_flat_time_msec' 'ppg_wout_noise' 'proportional_pulse_amp' 'review_bad_ppg_template' 'review_bad_r_template' 'sec_diastolic_2_diastolic' 'sec_dnotch_2_dnotch' 'sec_foot_2_foot' 'sec_from_last_onset_of_sleep' 'sec_neg_slope_after_2_neg_slope_after' 'sec_neg_slope_b4_2_neg_slope_b4' 'sec_pos_slope_2_pos_slope' 'sec_systolic_2_systolic' 'sec_to_first_onset_of_sleep' 'sec_tP_20_2_tP_20' 'sec_tP_20_inv_2_tP_20_inv' 'sec_tP_50_2_tP_50' 'sec_tP_50_inv_2_tP_50_inv' 'sec_tP_80_2_tP_80' 'sec_tP_80_inv_2_tP_80_inv' 'sec_tR_20_2_tR_20' 'sec_tR_20_inv_2_tR_20_inv' 'sec_tR_50_2_tR_50' 'sec_tR_50_inv_2_tR_50_inv' 'sec_tR_80_2_tR_80' 'sec_tR_80_inv_2_tR_80_inv' 'sec_valley_2_valley' 'sleep_stages' 'sqi_corrcoff_direct' 'sqi_corrcoff_interp' 'sqi_dtw' 'sqi_frechet' 'sqi_mean_corr_dtw' 'tP_20_x' 'tP_20_x_inv' 'tP_50_x' 'tP_50_x_inv' 'tP_80_x' 'tP_80_x_inv' 'tR_20_x' 'tR_20_x_inv' 'tR_50_x' 'tR_50_x_inv' 'tR_80_x' 'tR_80_x_inv'},{'aamp_foot' 'aamp_peak' 'slee_sleepstate' 'bamp_diapeak' 'bamp_dnotch' 'bamp_feet' 'bamp_madn' 'bamp_mbdn' 'bamp_tpm' 'bamp_systpeak' 'bamp_tp20' 'bamp_tp20inv' 'bamp_tp50' 'bamp_tp50inv' 'bamp_tp80' 'bamp_tp80inv' 'bamp_tr20' 'bamp_tr20inv' 'bamp_tr50' 'bamp_tr50inv' 'bamp_tr80' 'bamp_tr80inv' 'bamp_delsystol' 'ramp_diapeak' 'ramp_dnotch' 'ramp_feet' 'ramp_mdn' 'ramp_mbdn' 'ramp_tpm' 'ramp_systpeak' 'ramp_tp20' 'ramp_tp20inv' 'ramp_tp50' 'ramp_tp50inv' 'ramp_tp80' 'ramp_tp80inv' 'ramp_tr20' 'ramp_tr20inv' 'ramp_tr50' 'ramp_tr50inv' 'ramp_tr80' 'ramp_tr80inv' 'ramp_valley' 'area_area' 'barea_area' 'indx_begin' 'seco_time' 'bool_edgebeat' 'bool_errorppg' 'indx_intrabegin' 'indx_intradiastolic' 'indx_intradnotch' 'indx_intraend' 'indx_intrafoot' 'indx_intramadn' 'indx_intramadn' 'indx_intratpm' 'indx_intrasystolic' 'msec_length' 'msec_r2diapeak' 'msec_r2dnotch' 'msec_r2firval' 'msec_r2foot' 'msec_r2madn' 'msec_r2mbdn' 'msec_r2tpm' 'msec_r2secval' 'msec_r2syspeak' 'msec_r2tp20' 'msec_r2tp20inv' 'msec_r2tp50' 'msec_r2tp50inv' 'msec_r2tp80' 'msec_r2tp80inv' 'msec_r2tr20' 'msec_r2tr20inv' 'msec_r2tr50' 'msec_r2tr50inv' 'msec_r2tr80' 'msec_r2tr80inv' 'msec_totdur20' 'msec_totdur50' 'msec_totdur80' 'msec_totdur20inv' 'msec_totdur50inv' 'msec_totdur80inv' 'msec_tp502diepe' 'msec_tp502dnotc' 'msec_tp502firva' 'msec_tp502foot' 'msec_tp502sadn' 'msec_tp502sbdn' 'msec_tp502tpm' 'msec_tp502secva' 'msec_tp502speak' 'msec_tp502tp50' 'msec_tp502tp50inv' 'msec_tp502tp80' 'msec_tp502tp80inv' 'msec_tp502tr20' 'msec_tp502tr20inv' 'msec_tp502tr50' 'msec_tp502tr50inv' 'msec_tp502tr80' 'msec_tp502tr50inv' 'msec_ppgflattime' 'samp_ppgnoiseless' 'rati_proppulse' 'bool_badppgtemplate' 'bool_badecgtemplate' 'seco_dia2dia' 'seco_dn2dn' 'seco_foot2foot' 'seco_time2wake' 'seco_madn2madn' 'seco_mbdn2mbdn' 'seco_tpslo2slo' 'seco_sys2sys' 'seco_time2sleep' 'seco_tp202tp20' 'seco_tp202tp20inv' 'seco_tp202tp50' 'seco_tp202tp50inv' 'seco_tp202tp80' 'seco_tp202tp80inv' 'seco_tr202tr20' 'seco_tr202tr20inv' 'seco_tr502tr50' 'seco_tr502tr50inv' 'seco_tr802tr80' 'seco_tr802tr80inv' 'seco_val2val' 'slee_unadjustedsleepstage' 'psqi_corrdir' 'psqi_coorint' 'psqi_dtw' 'psqi_frechet' 'psqi_mcoordtw' 'indx_intratp20' 'indx_intratp20inv' 'indx_intratp50' 'indx_intratp50inv' 'indx_tp80' 'indx_intratp80inv' 'indx_intratr20' 'indx_intratr20inv' 'indx_intratr50' 'indx_intratr50inv' 'indx_intratr80' 'indx_intratr80inv'});


props = readProps('config.txt');
% metadata_path = props('Save_mesa_ids_matching');
features_path = props('Save_input');
outputLoc = props('Save_output');

FID = -1;
headers = [];

first = 0;
[analysisFiles] = SaveSetup(features_path, outputLoc);
time = 0;
parfor i = 1:size(analysisFiles,1)
%     tStart = tic;

%     avg_time = time/i;
%     disp(['Avg Time (s): ' num2str(avg_time)]);

%     disp(['Est finish (min): ' num2str((avg_time*(size(analysisFiles,1)-i))/60)]);
    disp(join(['Saving analysis of ' analysisFiles{i, 1}]));
    try
        beats_flattened = load(analysisFiles{i, 2});
    catch
        continue
    end
    beats_flattened = beats_flattened.beats_flattened;
%     beats_flattened = rmfield(beats_flattened,'ppg_wout_noise');
%     beats_flattened = rmfield(beats_flattened,'ppg_flat_time_msec');

    mkdir(fullfile(outputLoc,analysisFiles{i, 1}))
    headers = fieldnames(beats_flattened);

    for x = 1:length(headers)
%         if ~contains(headers{x},'idx') && %&& ~contains(headers{x},'sleep') % this does not output sleep variables need to export them also for proper analysis
%          if isKey(dict, headers{x})
              
        output_file = fullfile(fullfile(outputLoc, analysisFiles{i, 1}), [analysisFiles{i, 1} '-'  headers{x} '.mat']);
        if isfile(output_file)
            continue
        end
        tmp = beats_flattened.(headers{x});
        parsave(output_file,tmp);
%         end
    end
    

%     disp(['____________________________________________________________________________________________________' newline]);
%     time = time+toc(tStart);
%     toc(tStart);
end

function structOut = renamefield(structIn, oldField, newField)                         
    for i = 1:length(structIn)          
        structIn = setfield(structIn,{i},newField,getfield(structIn(i),oldField));                
    end         
    structOut = rmfield(structIn,oldField);                    
end

function parsave(location, tmp)
    save(location, 'tmp');
end