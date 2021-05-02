function this = FitTrapezoidsForMeanMS2Traces(this, SetIndex, NC, APindex, TraceType)

alpha = this.alpha;

[MeanTrace, StdErrorTrace, NCTimes] = GetSpecifiedTracesForFitting(this, SetIndex, NC, APindex, TraceType);

if length(MeanTrace) >= this.MinimumTimePoints
    pos_slope = NaN;
    t_elongation = NaN;
    pos_yintercept = NaN;
    pos_xintercept = NaN;
    t_off = NaN;
    t_peak = NaN;
    span = ceil(1/this.time_delta)*2+1;
    SmoothedTrace = smooth(MeanTrace, span);%, 'loess');
    PositiveFitIncludedIndices = GetIndicesForPositiveSlopeFit(SmoothedTrace);
    
    if (length(PositiveFitIncludedIndices)>= this.MinimumFittingPoints) 
        [pos_fitresult,pos_gof1,pos_out1] = fit(NCTimes(PositiveFitIncludedIndices).',MeanTrace(PositiveFitIncludedIndices).','poly1');
        pos_slope = pos_fitresult.p1;
        ci = confint(pos_fitresult, alpha);
        ci_slope = ci(:,1).';
        t = tinv((1+alpha)/2, 1);
        se_pos_slope = (ci_slope(2)-ci_slope(1)) ./ (2*t); % Standard Error
        pos_yintercept = pos_fitresult.p2;
        ci_yintercept = ci(:,2).';
        t = tinv((1+alpha)/2, 1);
        se_pos_yintercept = (ci_yintercept(2)-ci_yintercept(1)) ./ (2*t); % Standard Error
        %         this = SetSpecifiedTrapezoidFitParameters(this, pos_slope, se_pos_slope, ci_slope,...
        %             SetIndex, NC, APindex, TraceType, 'MeanInitiationRates');
        
        [pos_fitresult2,pos_gof2,pos_out2] = fit(MeanTrace(PositiveFitIncludedIndices).',NCTimes(PositiveFitIncludedIndices).','poly1');
        pos_xintercept = pos_fitresult2.p2;
        ci2 = confint(pos_fitresult2, alpha);
        ci_xintercept = ci2(:,2).';
        t = tinv((1+alpha)/2, 1);
        se_pos_xintercept = (ci_xintercept(2)-ci_xintercept(1)) ./ (2*t); % Standard Error
        %         this = SetSpecifiedTrapezoidFitParameters(this, pos_xintercept, se_pos_xintercept, ci_xintercept,...
        %             SetIndex, NC, APindex, TraceType, 'TimeOns');
        
        
    else
        pos_slope = NaN;
        se_pos_slope = NaN;
        pos_yintercept = NaN;
        se_pos_yintercept = NaN;
        pos_xintercept = NaN;
        se_pos_xintercept = NaN;
       
    end
    
    
    NegativeFitIncludedIndices = GetIndicesForNegativeSlopeFit(SmoothedTrace);
    if NC < 14
        if (length(NegativeFitIncludedIndices)>= this.MinimumFittingPoints)
            [neg_fitresult,neg_gof1,neg_out1] = fit(NCTimes(NegativeFitIncludedIndices).',MeanTrace(NegativeFitIncludedIndices).','poly1');
            neg_slope = neg_fitresult.p1;
            ci = confint(neg_fitresult, alpha);
            ci_neg_slope = ci(:,1).';
            t = tinv((1+alpha)/2, 1);
            se_neg_slope = (ci_neg_slope(2)-ci_neg_slope(1)) ./ (2*t); % Standard Error
            neg_yintercept = neg_fitresult.p2;
            ci_neg_yintercept = ci(:,2).';
            t = tinv((1+alpha)/2, 1);
            se_neg_yintercept = (ci_neg_yintercept(2)-ci_neg_yintercept(1)) ./ (2*t); % Standard Error
            %             this = SetSpecifiedTrapezoidFitParameters(this, neg_slope, se_neg_slope, ci_neg_slope,...
            %                 SetIndex, NC, APindex, TraceType, 'UnloadingRates');
        else
            neg_slope = NaN;
            se_neg_slope = NaN;
            neg_yintercept = NaN;
            se_neg_yintercept = NaN;
        end
        
        
        
        
        
        
        %plot(NCTimes, SmoothedTrace)
        
        
        
        
        if ~isempty(PositiveFitIncludedIndices) & ~isempty(NegativeFitIncludedIndices)
            PlateauFitIncludedIndices = max(PositiveFitIncludedIndices):min(NegativeFitIncludedIndices);
            if length(PlateauFitIncludedIndices) >= 2
                [plateau_height, plateau_S, plateau_mu] = polyfit(NCTimes(PlateauFitIncludedIndices).',MeanTrace(PlateauFitIncludedIndices).',0);
                ci_plateau = polyparci(plateau_height, plateau_S, alpha).';
                t = tinv((1+alpha)/2, 1);
                se_plateau = (ci_plateau(2)-ci_plateau(1)) ./ (2*t); % Standard Error
            elseif length(PlateauFitIncludedIndices) == 1
                plateau_height = MeanTrace(PlateauFitIncludedIndices);
                se_plateau = NaN;
            else
                plateau_height = NaN;
                se_plateau = NaN;
            end
            if all(~isnan([pos_slope, pos_yintercept, plateau_height]))
                [t_peak, ~] = calculateIntersection(pos_slope, pos_yintercept, 0, plateau_height);
                t_elongation = t_peak-pos_xintercept;
                %                 this = SetSpecifiedTrapezoidFitParameters(this, t_elongation, [], [],...
                %                     SetIndex, NC, APindex, TraceType, 'ElongationTimes');
            else
                t_peak = NaN;
                t_elongation = NaN;
            end
            if all(~isnan([neg_slope, neg_yintercept, plateau_height]))
                [t_off, ~] = calculateIntersection(neg_slope, neg_yintercept, 0, plateau_height);
                
            else
                t_off = NaN;
            end
        end
    else
        neg_slope = NaN;
        if ~isempty(PositiveFitIncludedIndices) & ~isempty(NegativeFitIncludedIndices)
            PlateauFitIncludedIndices = max(PositiveFitIncludedIndices):min(NegativeFitIncludedIndices);
            if length(PlateauFitIncludedIndices) >= 2
                [plateau_height, plateau_S, plateau_mu] = polyfit(NCTimes(PlateauFitIncludedIndices).',MeanTrace(PlateauFitIncludedIndices).',0);
                ci_plateau = polyparci(plateau_height, plateau_S, alpha).';
                t = tinv((1+alpha)/2, 1);
                se_plateau = (ci_plateau(2)-ci_plateau(1)) ./ (2*t); % Standard Error
            elseif length(PlateauFitIncludedIndices) == 1
                plateau_height = MeanTrace(PlateauFitIncludedIndices);
                se_plateau = NaN;
            else
                plateau_height = NaN;
                se_plateau = NaN;
            end
            if all(~isnan([pos_slope, pos_yintercept, plateau_height]))
                [t_peak, ~] = calculateIntersection(pos_slope, pos_yintercept, 0, plateau_height);
                t_elongation = t_peak-pos_xintercept;
                %                 this = SetSpecifiedTrapezoidFitParameters(this, t_elongation, [], [],...
                %                     SetIndex, NC, APindex, TraceType, 'ElongationTimes');
            else
                t_peak = NaN;
                t_elongation = NaN;
            end
            
            t_off = NCTimes(min(NegativeFitIncludedIndices));
            
        end
    end
    %
   
    
    
    if ~isempty(PositiveFitIncludedIndices)
        StartIndex = PositiveFitIncludedIndices(1);
    else
        StartIndex = 1;
    end
    if ~isempty(NegativeFitIncludedIndices) & NC < 14
        EndIndex = NegativeFitIncludedIndices(end);
    elseif ~isempty(NegativeFitIncludedIndices) & NC == 14
        EndIndex = NegativeFitIncludedIndices(1);
    else
        EndIndex = length(NCTimes);
    end
    
    if ~isempty(PositiveFitIncludedIndices) & ~isempty(NegativeFitIncludedIndices) & (NC < 14)
        ft = fittype( 'trapezoidFitFunction(x, a, b, c, t1, t2 )' );
        
        if EndIndex-StartIndex+1 >= this.MinimumTimePoints
             StartingPoints = zeros(1, 5);
             if ~isnan(pos_slope)
                 StartingPoints(1) = pos_slope;
             else
                 StartingPoints(1) = 1;
             end
             if ~isnan(pos_slope) & ~isnan(pos_xintercept)
                 StartingPoints(2) = -pos_slope*pos_xintercept;
             end
             if ~isnan(neg_slope)
                 StartingPoints(3) = neg_slope;
             end
             if ~isnan(t_elongation)& ~isnan(pos_xintercept)
                 StartingPoints(4) = t_elongation+pos_xintercept;
             end
             if ~isnan(t_off)
                 StartingPoints(5) = t_off;
             end
            [f, gof] = fit( NCTimes(StartIndex:EndIndex).', MeanTrace(StartIndex:EndIndex).', ft,...
                'StartPoint', StartingPoints,...
                'Lower', [0 -Inf -max(MeanTrace(StartIndex:EndIndex))  0 0],...
                'Upper', [max(MeanTrace(StartIndex:EndIndex)) 0 0 max(NCTimes) max(NCTimes) ]);
            
            
            pos_slope = f.a;
            try
                ci = confint(f);
            catch
                ci = NaN(2,5);
            end
            ci_slope = ci(:,1).';
            t = tinv((1+alpha)/2, 1);
            se_pos_slope = (ci_slope(2)-ci_slope(1)) ./ (2*t); % Standard Error
            pos_yintercept =f.b;
            ci_yintercept = ci(:,2).';
            t = tinv((1+alpha)/2, 1);
            se_pos_yintercept = (ci_yintercept(2)-ci_yintercept(1)) ./ (2*t); % Standard Error
            
            
            pos_xintercept = -pos_yintercept/pos_slope;
            
            se_pos_xintercept = sqrt((1/pos_slope)^2*se_pos_yintercept^2+(pos_yintercept/(pos_slope^2))^2*se_pos_slope^2);
            
            neg_slope = f.c;
            
            ci_neg_slope = ci(:,3).';
            t = tinv((1+alpha)/2, 1);
            se_neg_slope = (ci_neg_slope(2)-ci_neg_slope(1)) ./ (2*t); % Standard Error
            %     neg_yintercept = neg_fitresult.p2;
            %     ci_neg_yintercept = ci(:,2).';
            %     t = tinv((1+alpha)/2, 1);
            %     se_neg_yintercept = (ci_neg_yintercept(2)-ci_neg_yintercept(1)) ./ (2*t); % Standard Error
            
            t_off = f.t2;
            ci_t_off = ci(:,5).';
            t = tinv((1+alpha)/2, 1);
            se_t_off = (ci_t_off(2)-ci_t_off(1)) ./ (2*t); % Standard Error
            
            
            t_peak = f.t1;
            ci_t_peak = ci(:,4).';
            t = tinv((1+alpha)/2, 1);
            se_t_peak = (ci_t_peak(2)-ci_t_peak(1)) ./ (2*t); % Standard Error
            t_elongation = t_peak-pos_xintercept;
            se_t_elongation = sqrt(se_t_peak^2 + se_pos_xintercept^2);
            
            this = SetSpecifiedTrapezoidFitParameters(this, f, [], [],...
                SetIndex, NC, APindex, TraceType, 'Fits');
            this = SetSpecifiedTrapezoidFitParameters(this, pos_slope, se_pos_slope, ci_slope,...
                SetIndex, NC, APindex, TraceType, 'MeanInitiationRates');
            this = SetSpecifiedTrapezoidFitParameters(this, pos_xintercept, se_pos_xintercept, NaN(1,2),...
                SetIndex, NC, APindex, TraceType, 'TimeOns');
            
            this = SetSpecifiedTrapezoidFitParameters(this, neg_slope, se_neg_slope, ci_neg_slope,...
                SetIndex, NC, APindex, TraceType, 'UnloadingRates');
            this = SetSpecifiedTrapezoidFitParameters(this, t_elongation, se_t_elongation, NaN(1,2),...
                SetIndex, NC, APindex, TraceType, 'ElongationTimes');
            this = SetSpecifiedTrapezoidFitParameters(this, t_off, se_t_off, ci_t_off,...
                SetIndex, NC, APindex, TraceType, 'TimeOffs');
            this = SetSpecifiedTrapezoidFitParameters(this, gof.rsquare, [], [],...
                SetIndex, NC, APindex, TraceType, 'R2s');
        end
        %end
    elseif ~isempty(PositiveFitIncludedIndices)
        ft2 = fittype( 'leftHalfTrapezoidFitFunction(x, a, b, t1)' );
        
        if EndIndex-StartIndex+1 >= this.MinimumTimePoints
             StartingPoints = zeros(1, 3);
             if ~isnan(pos_slope)
                 StartingPoints(1) = pos_slope;
             else
                 StartingPoints(1) = 1;
             end
             if ~isnan(pos_slope) & ~isnan(pos_xintercept)
                 StartingPoints(2) = -pos_slope*pos_xintercept;
             end
             if ~isnan(t_elongation)& ~isnan(pos_xintercept)
                 StartingPoints(3) = t_elongation+pos_xintercept;
             end
       
            [f, gof] = fit( NCTimes(StartIndex:EndIndex).', MeanTrace(StartIndex:EndIndex).', ft2,...
                'StartPoint', StartingPoints,...
                'Lower', [0 -Inf 0],...
                'Upper', [max(MeanTrace(StartIndex:EndIndex)) 0 max(NCTimes)]);
            
            
            pos_slope = f.a;
            try
                ci = confint(f);
            catch
                ci = NaN(2,3);
            end
            ci_slope = ci(:,1).';
            t = tinv((1+alpha)/2, 1);
            se_pos_slope = (ci_slope(2)-ci_slope(1)) ./ (2*t); % Standard Error
            pos_yintercept =f.b;
            ci_yintercept = ci(:,2).';
            t = tinv((1+alpha)/2, 1);
            se_pos_yintercept = (ci_yintercept(2)-ci_yintercept(1)) ./ (2*t); % Standard Error
            
            
            pos_xintercept = -pos_yintercept/pos_slope;
            
            se_pos_xintercept = sqrt((1/pos_slope)^2*se_pos_yintercept^2+(pos_yintercept/(pos_slope^2))^2*se_pos_slope^2);
            
          
            
            
            t_peak = f.t1;
            ci_t_peak = ci(:,3).';
            t = tinv((1+alpha)/2, 1);
            se_t_peak = (ci_t_peak(2)-ci_t_peak(1)) ./ (2*t); % Standard Error
            t_elongation = t_peak-pos_xintercept;
            se_t_elongation = sqrt(se_t_peak^2 + se_pos_xintercept^2);
            
            
            this = SetSpecifiedTrapezoidFitParameters(this, f, [], [],...
                SetIndex, NC, APindex, TraceType, 'Fits');
            this = SetSpecifiedTrapezoidFitParameters(this, pos_slope, se_pos_slope, ci_slope,...
                SetIndex, NC, APindex, TraceType, 'MeanInitiationRates');
            this = SetSpecifiedTrapezoidFitParameters(this, pos_xintercept, se_pos_xintercept, NaN(1,2),...
                SetIndex, NC, APindex, TraceType, 'TimeOns');
         
            this = SetSpecifiedTrapezoidFitParameters(this, t_elongation, se_t_elongation, NaN(1,2),...
                SetIndex, NC, APindex, TraceType, 'ElongationTimes');

            this = SetSpecifiedTrapezoidFitParameters(this, gof.rsquare, [], [],...
                SetIndex, NC, APindex, TraceType, 'R2s');
            %end
        end
        
    end
    
    
    
end
