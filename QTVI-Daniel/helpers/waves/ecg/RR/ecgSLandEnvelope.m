function listcycok = ecgSLandEnvelope(rwave,sampling,EnvelopeBasesMsec_MP,startingAvgUnitsIndexes,startingStdMsec,MinMaxRRintMsecAllowed )

% 2/17/13 - need more outlier removal other than SL and L and VERYshort
%               - see Ron Berger viz Larisa emails, Feb 2013, see data file++ 76_DUR_ESM.ecg
% 2/17/13 - REM - this is one place that outlier detection is done,
%           but maybe MORE IMPORTANT is :
%                   - thruRsB.m in SHORT FILE ANALYSIS
%                   - baseqt_gissipheno.m section in HOLTER ANALYSIS
% 2/24/13 - only generate <=-3.99, DO NOT generate a -4 in either this program _OR_ ecgSL.m


if nargin<5
    error('ecgSL\\andEnvelope\\.m needs at least 5 inputs : - quitting.');
end

if exist('MinMaxRRintMsecAllowed','var')~=1
    MinMaxRRintMsecAllowed = [ 0     10000 ] ;
end
if numel( MinMaxRRintMsecAllowed )~=2
    MinMaxRRintMsecAllowed = [ 0     10000 ] ;
end


% this code/key won't formally be linked to the calling program,
% but the code/key should be universal in all programs

shortofSLsetmark = -1 ;  % using 1-RRlimperc & 1+RRlimperc
longofSLsetmark  = -2 ;  % -1,-2 always a pair
justlongmark     = -3 ;  % using RRlimbig
shortPVCenvelmark= -3.25 ; % added 2/24/13
longPVCenvelmark = -3.50 ; % added 2/24/13
%------------- ^ only these are PVC's ^
%QTshortmark     = -4 ;  % added "big new thing" 02/17/13 , REMOVED 2/24/13
%QTsfORmark      = -5 ;  % using QTsffromRM
%QThiNssvmark    = -6 ;  % using Nssv
%cycflatmark     = -7 ;  % if submitted cycle is flat



% this too is universal :
RRlimperc = 0.12; % meant for when an (early) PVC is detected as an R wave - 10%, 12%, . . .
RRlimbig  = 1.60; % meant when a PVC is not detected as an R wave - - i should go very close to 2xRRint as a length
%%%don't do : PH1 = 1.20; PH2 = .80;

% put code here that runs thru R waves, finds PVCs based on S-L, justL criteria
Nbackcycavg = 15 ; %%% 10 ; % 2/17/13 - increase to 15, to get a better std value

everyLseg     = diff(rwave);
everyLsegMSEC = 1000*everyLseg/sampling ;
lenrwave      = numel(rwave);
if numel(everyLseg) < Nbackcycavg
    Nbackcycavg = numel(everyLseg)-2;
end

uptorwave_i = lenrwave-2 ;
rwaverng    = 1:uptorwave_i ;
listcycok   = ones(1,uptorwave_i); % pre-make

% speed is an issue with q2off_cd, so do this :
prevavgRRint = startingAvgUnitsIndexes ; % <--- better ; %%% PREV : median(everyLseg(1:Nbackcycavg));	% MEDIAN, NOT AVERAGE - this hasn't been investiagted yet, may have PVCs in here
SUMlast10L   = Nbackcycavg*( prevavgRRint ) ;
BLOCKlast10L = prevavgRRint ;
BLOCKlast10L = BLOCKlast10L( ones(1,Nbackcycavg) ) ;
%
prevavgRRintMSEC = 1000 * prevavgRRint / sampling ;
prevSTDrrintMSEC = startingStdMsec ; % from an immed. prev execution of "old" ecgSL.m in the caller program
MultFactorOnStd = 2.0 ;

%%%%%%%%%%%%%%%%%%%InDevelopOnly_StdVector=[];


%last10L = ones(1,Nbackcycavg) * ( median(everyLseg(1:Nbackcycavg)) ) ; 	      prevavgRRint = mean(last10L); % will use rem2 to sub-in new values in circular indexing, just taking avg, order not important, etc
%last10H = ones(1,Nbackcycavg) * ( median(datasignal(rwave(1:Nbackcycavg))) ) ; 	prevavgRhite = mean(last10H); % will use rem2 to sub-in new values in circular indexing, just taking avg, order not important, etc

nextcyc_ok='y';
n_ok=0;


%-------------------------
for perrwave = rwaverng
    thiscyc_ok = nextcyc_ok;    
    nextcyc_ok='y';  % all of the indexes inside
    
    Lseg     = everyLseg(perrwave);    
    Lnextseg = everyLseg(perrwave+1);
    LsegMSEC = everyLsegMSEC(perrwave);
    
    try_these_indexes = perrwave+1  :  perrwave+10 ;
    try_these_indexes = try_these_indexes(try_these_indexes<=numel(LsegMSEC)) ;
    
    thesefutureRRintMsec  = LsegMSEC(try_these_indexes) ;
    if numel(thesefutureRRintMsec)==0
        thesefutureRRintMsec = prevavgRRintMSEC  ; 
    end  % or = LsegMSEC
    
    thisMinusEnvelopeMsec = abs(  EnvelopeBasesMsec_MP(1)  +  MultFactorOnStd*prevSTDrrintMSEC  ) ;
    thisPlusEnvelopeMsec  = abs(  EnvelopeBasesMsec_MP(2)  +  MultFactorOnStd*prevSTDrrintMSEC  ) ;
    
    %&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    % L of SL sequence :
    if thiscyc_ok=='n', 									listcycok(perrwave) = longofSLsetmark ;
        
        % S of SL sequence :
    elseif Lseg <= prevavgRRint*(1-RRlimperc) && Lnextseg >= prevavgRRint*(1+RRlimperc)
        listcycok(perrwave) = shortofSLsetmark;		% but don't add to Npvc's
        nextcyc_ok='n'; % <- only time this can be set to 'n'
        
        % just Long, infer must contain a S-L set
    elseif Lseg >= prevavgRRint*(RRlimbig),                 listcycok(perrwave) = justlongmark ;
        
        
        % new 2/17/13 : just Short (via Envelope)
    elseif LsegMSEC                <(prevavgRRintMSEC-thisMinusEnvelopeMsec)   ...
            && any(thesefutureRRintMsec>(prevavgRRintMSEC-thisMinusEnvelopeMsec/2))
        listcycok(perrwave) = shortPVCenvelmark ;
        
        % new 2/17/13 : just Long-via-Envelope
    elseif LsegMSEC                >(prevavgRRintMSEC+thisPlusEnvelopeMsec)   ...
            && any(thesefutureRRintMsec<(prevavgRRintMSEC+thisPlusEnvelopeMsec/2))
        listcycok(perrwave) = longPVCenvelmark ;
        
        
        % new 3/3/13 :
    elseif LsegMSEC < MinMaxRRintMsecAllowed(1)
        listcycok(perrwave) = shortPVCenvelmark ;
        
        % new 3/3/13 :
    elseif LsegMSEC > MinMaxRRintMsecAllowed(2)
        listcycok(perrwave) = longPVCenvelmark ;
        
        
    end% elseif's
    %&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    
    
    % still inside of the perrwave loop, for this 1 cycle :
    
    
    % an if has to be inside the loop, results here may [or may not] affect the next iteration of for
    
    
    %|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    if listcycok(perrwave)==1
        
        n_ok=n_ok+1;
        nowi_in_block = rem2(n_ok,Nbackcycavg);
        
        SUMlast10L   = SUMlast10L  - BLOCKlast10L(nowi_in_block)  +   Lseg  ;
        prevavgRRint = SUMlast10L / Nbackcycavg ;
        
        prevavgRRintMSEC = 1000 * prevavgRRint / sampling ;
        
        
        %%%last10L( rem2(perrwave,Nbackcycavg) ) = Lseg ;  								prevavgRRint = mean(last10L);
        %%%last10H( rem2(perrwave,Nbackcycavg) ) = datasignal(rwave(perrwave)) ;   prevavgRhite = mean(last10H);
        
        % and replace in block
        BLOCKlast10L(nowi_in_block) = Lseg ;   % the next time, when i do n_ok=n_ok+1, now's will be at i-1, and thereby be the oldest
        
        % only compute/update STD if a full set of the block is "real" and good :
        if n_ok>=Nbackcycavg
            prevSTDrrintMSEC = 1000 * std(BLOCKlast10L) / sampling ;
            %%%%%%%%%%%%%%%%%%%InDevelopOnly_StdVector = [ InDevelopOnly_StdVector   prevSTDrrintMSEC ] ;
        end % if n_ok>=Nbackcycavg,
        
        
        % but more importantly :
        %%% only for validation of "SW" SUM/BLOCK avg'ing, remove when done :
        %%%HOLDprevavgRRint(perrwave) = prevavgRRint ;
        % yes, works by plot(everyLseg,'b') ; plot(HOLDprevavgRRint,'r') , *** -> changes in RRint are tracked by average
        
        
    end  % if listcycok(perrwave)==1,
    %|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    
    
    
    
    
    
end% for perrwave
%-------------------------



%%%%%%%%%%%%%%%%%%%save InDevelop.mat ; error('crash inside of ecgSLandEnvelope.m');

% NOTE : numel(listcycok) = numel(rwave) - 2


% 10/14/15 : found on 2 Sleep cases : 107456, 107676 :
% many many PVC_Long's in a row - "must" be the case of an instantaneous jump to a new normal of longer RRints
% do not have this one ex post facto "fix" ruin everything for all cases except these, do a light touch :
% won't fix all of them, but will fix a lot of them : will miss the first and last few, is OK :
% have a big barrier to this - must be a full 10 in a row :
% !!!! must do -10 in the for declare, else the +9 will exceed limits !

nnnn = min([ numel(listcycok)  numel(everyLsegMSEC) ]) ;
if nnnn>30
    TENONES = ones(1,10) ; % make once not 4 times x ______
    for perchunk = 1:10:nnnn-10
        if  all(listcycok(perchunk:perchunk+9)==     justlongmark*TENONES)==1 ...
                ||  all(listcycok(perchunk:perchunk+9)==shortPVCenvelmark*TENONES)==1 ...
                ||  all(listcycok(perchunk:perchunk+9)== longPVCenvelmark*TENONES)==1
            
            % orig : all 1 shot :  listcycok(perchunk:perchunk+9)=              (+1)*TENONES ;
            
            % additional test before making good :
            for perperadd=0:9
                if everyLsegMSEC(perchunk+perperadd) > MinMaxRRintMsecAllowed(1) ...
                        && everyLsegMSEC(perchunk+perperadd) < MinMaxRRintMsecAllowed(2) 
                    listcycok(perchunk+perperadd)=  +1 ;
                end 
            end 
           
        end 
    end 
end 


