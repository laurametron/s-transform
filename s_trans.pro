; The S transform function
; written by bob stockwell, 1993
;
; Returns a matrix if succesful
; a structure  if a sampling rate is given
;    The structure is {st:, Freq:, Time:}
;	 where st is the S Transform, freq is a vector containing the
;    frequencies
;
;    Optional Paramters
;   WIDTH   size of time resolution
;			if not set, width = 1
;
;
;    Keywords
;   \HELP				explains all the keywords and parameters
;   \VERBOSE			flags errors and size
;	\SAMPLINGRATE   	if set returns array of frequency
;	\MAXFREQ
;	\MINFREQ
;	\FREQSAMPLINGRATE
;	\PIECEWISENUMBER    divides the time series, and passes back array
; 	\POWER				 returns the power spectrum
; 	\AMPLITUDE				 returns the absolute value spectrum
;	\REMOVEEDGE		 removes the edge with a 5% taper, and takes
;							 out least-sqares fit parabola
;	\NO_ANALYTIC      suppresses the analytic signal that is performed on a real time series
;	\NEGATIVE_FREQUENCIES calculates and returns the negative frequency components.

; added a "mask out edges" keyword,
; which STs a line (ie st of edges) and thresholds the returned st matrix.
; The value of masked edges is the percent at which to  make the threshold
; default = 5 %.

; added an EXAMPLE keyword, will display a time series and the
; amplitude of the ST to the current graphics device
; WARNING, will call !P.multi=0

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; modified hilbert transform ;;;;;;;;;;;;;;;;;;;;;
;+
; NAME:
;       HILBERT
;
; PURPOSE:
;       Return a series that has all periodic terms shifted by 90 degrees.
;

;^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
; REVISION HISTORY:
;       JUNE, 1985,     Written, Leonard Kramer, IPST (U. of Maryland) on site
;                       contractor to NASA(Goddard Sp. Flgt. Cntr.)
;-
;================================================================================
;  MODIFIED:   APRIL 20 1994 by RGS
;                  - returns real function only
;                  - added keyword     /ANALYTIC
;                    to return the analytic signal
;^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
; the S Transform function
FUNCTION s_trans, ts, factor, HELP = HELP,VERBOSE=VERBOSE,SAMPLINGRATE=SAMPLINGRATE $
    ,MAXFREQ=MAXFREQ,MINFREQ=MINFREQ  $
    ,FREQSAMPLINGRATE=FREQSAMPLINGRATE     $
    ,POWER=POWER,AMPLITUDE=AMPLITUDE,REMOVEEDGE=REMOVEEDGE   $
    , maskedges=maskedges,EXAMPLE=EXAMPLE $
    ,NO_ANALYTIC=NO_ANALYTIC,NEGATIVE_FREQUENCIES=NEGATIVE_FREQUENCIES,abs=abs

; BACkWARDS COMPATIBILTY, i changed the name of the keyword
if keyword_set(abs) then begin
	AMPLITUDE = abs
	print,'keyword (abs) obsolete, use AMPLITUDE'
endif


if keyword_set(EXAMPLE) then begin  ; show example of a cos of chirp function
	ex_len = 512
	ex_time = findgen(ex_len)
	ex_freq = 5;ex_len/16
	ex_ts =  cos(2*!Pi*ex_freq*ex_time/ex_len)
	ex_ts =  cos(2*!Pi*(ex_len/5+2*ex_ts)*ex_time/ex_len)
	; crossed chirp example commented out
		;ex_ts =  cos(2*!Pi*ex_freq*ex_time*(1+2*ex_time/ex_len)/ex_len)
		;ex_ts = ex_ts + reverse(ex_ts)
;	!P.multi=[0,1,2]
;	plot,ex_ts,xtitle='Time (units)',title='Time Series [h(t) = cos(cos(wt))]'
	s = s_trans(ex_ts,/samp, /AMPLITUDE,verbose=verbose)  ; returns structure, amps only returned
	nlevels = 14
	levels = findgen(nlevels)/(nlevels-1)*1.5

;	contour,s.st,ex_time,s.freq,levels=levels,/fill,xtitle='Time (units)', $
;		ytitle='Frequency (1/unit)',title='Amplitude of S-Transform'
	return,0
	!P.multi=0
endif

if keyword_set(HELP) then begin
	Print,"S_TRANS()  HELP COMMAND"
	Print,"S_trans() returns a matrix if succesful or a structure if required"
	Print,"S_trans() returns  - 1 or an error message if it fails"
	Print,"USAGE::    localspectra = s_trans(timeseries)
	Print," "
	Print,"Optional Parameters"
	Print,"WIDTH  -size of time resolution"
	Print,"       -if not set, default WIDTH = 1"
	Print," "
	Print,"Keywords:
	Print,"\HELP            -explains all the keywords and parameters
	Print,"\VERBOSE         -flags errors and size, tells time left etc.
	Print,"\SAMPLINGRATE    -if set returns array of frequency
	Print,"\MAXFREQ
	Print,"\MINFREQ
	Print,"\FREQSAMPLINGRATE
	Print,"\POWER           -returns the power spectrum
	Print,"\AMPLITUDE             -returns the absolute value spectrum
	Print,"\REMOVEEDGE      -removes the edge with a 5% taper, and takes
	Print,"                 -out least-squares fit parabola
	Print,'\NO_ANALYTIC     - does not perform analytic signal of real-valued time series.'
	return, -1
endif

;print,'paramters',N_params()
;Check number of arguments.
CASE N_PARAMS() OF
	1: begin
		if n_elements(ts) eq 0 then MESSAGE, 'Invalid timeseries (check your spelling).'
		factor = 1
	end
	2: if n_elements(factor) eq 1 then begin ;Two-argument case.
			; NOW factor as an array is accepted
		 	; factor = 1                          ;Make sure factor is a number.
	   		; if keyword_set(VERBOSE)then print,'Error in second parameter. Using default values.'
	   endif
ELSE: MESSAGE, 'Incorrect number of arguments'
ENDCASE


if n_elements(ts) eq 0 then MESSAGE, 'Invalid timeseries (check your spelling).'
time_series = ts ; don't change the original dataset passed to this function

; check to see if it is a vector, not a 1 x N matrix
sz = size(time_series)
if sz(0) ne 1 then begin
    if sz(1) eq 1 and sz(2) gt 1 then begin
       time_series = reform(time_series)  ; a column vector, change it
       if keyword_set(VERBOSE)then print,'Reforming timeseries'
	endif else MESSAGE, 'Must enter an array of data'
endif

if keyword_set(VERBOSE) then print
if keyword_set(VERBOSE) then print,'Performing S transform:'

if keyword_set(REMOVEEDGE)  then begin
	 if keyword_set(VERBOSE) then  print,'Removing edges, LSF parabola removed, 5% taper on ends applied.'
 	 ind = findgen(n_elements(time_series))
 	 tstype = size(time_series,/type)
 	 if (tstype eq 6) or (tstype eq 9) then begin ; complex input
		rr = poly_fit(ind,double(time_series),2,rfit,yband,sigma,AM)
		ri = poly_fit(ind,imaginary(time_series),2,ifit,yband,sigma,AM)
		time_series = time_series - complex(rfit,ifit)
 	 endif else begin
 		 r = poly_fit(ind,time_series,2,fit,yband,sigma,AM)
		 time_series = time_series - fit
	 endelse
 	 ;ts_power = sqrt(total(abs(time_series)^2)/n_elements(time_series))
 	 ;if keyword_set(VERBOSE) then  print,'Sigma is:',sigma,'power',ts_power,'ratio',sigma/ts_power
 	 ;if keyword_set(VERBOSE) then print, 'total error',total(yband)/n_elements(yband)
 	 sh_len = n_elements(time_series)/10 ; 5% window on the result
	 if sh_len gt 1 then begin
	 	wn = hanning(sh_len)
 	 	time_series(0:sh_len/2-1) = time_series(0:sh_len/2-1)*wn(0:sh_len/2-1)
 	 	time_series(n_elements(time_series)-sh_len/2:*) = time_series(n_elements(time_series)-sh_len/2:*)*wn(sh_len/2:*)
 	endif
endif


; here its dimension is one
sz = size(time_series)
if (sz(2) ne 6) and (sz(2) ne 9) then begin
	if keyword_set(VERBOSE) then print,'Not complex data, finding analytic signal.'
	;take hilbert transfrom
	; do NOT do hilbert if we want neg freqs, or else they will all be  = zero
	if keyword_set(NEGATIVE_FREQUENCIES) then NO_ANALYTIC = 1 
	if not(keyword_set(NO_ANALYTIC)) then begin
		time_series = (hilbert(time_series,/analytic))
	endif
endif


length = n_elements(time_series)
spe_length = length/2
b = dcomplexarr(length)
gw = dcomplexarr(length)
h = fft(time_series,-1)

; do the different sampling cases here:
if (keyword_set(MAXFREQ))  then begin
	if maxfreq lt 1 then begin
		if keyword_set(SAMPLINGRATE) then begin  ; make structure
			maxfreq = long(length*samplingrate*maxfreq)
		endif else begin
			maxfreq = long(length*maxfreq)
		endelse
	endif
endif
if not(keyword_set(MINFREQ))  then begin
		if keyword_set(VERBOSE) then print,'Minimum Frequency is 0.'
		MINFREQ = 0  ; loop starts at 0
endif else begin
		if MINFREQ gt spe_length then begin
				MINFREQ = spe_length
				print,'MINFREQ too large, using default value'
		endif
		if minfreq lt 1 then begin
			if keyword_set(SAMPLINGRATE) then begin  ; make structure
				minfreq = long(length*samplingrate*minfreq)
			endif else begin
				minfreq = long(length*minfreq)
			endelse
		endif
		if keyword_set(VERBOSE) then print,strcompress('Minimum Frequency is '+string(MINFREQ)+'.')
endelse
if not(keyword_set(MAXFREQ))  then begin
		if keyword_set(VERBOSE) then print,strcompress('Maximum Frequency is '+string(spe_length)+'.')
		MAXFREQ = spe_length
endif else begin
		if MAXFREQ gt spe_length then begin
				MAXFREQ = spe_length
				print,'MAXFREQ too large, using default value'
		endif
		if keyword_set(VERBOSE) then print,strcompress('Maximum Frequency is '+string(MAXFREQ)+'.')
endelse
if not(keyword_set(FREQSAMPLINGRATE))  then begin
		if keyword_set(VERBOSE) then print,'Frequency sampling rate is 1.'
		FREQSAMPLINGRATE = 1
endif else if keyword_set(VERBOSE) then print,strcompress('Frequency sampling rate is '+string(FREQSAMPLINGRATE)+'.')

if FREQSAMPLINGRATE eq 0 then FREQSAMPLINGRATE = 1     ; if zero then use default

; check for errors in frequency parameters
if MAXFREQ lt MINFREQ then begin  ; if min > max, switch them
	temp = MAXFREQ
	MAXFREQ = MINFREQ
    MINFREQ= temp
    temp = 0
    print,'Switching frequency limits.'+strcompress(' Now, (MINFREQ = '+string(MINFREQ) + ') and (MAXFREQ ='+string(MAXFREQ)+').')
endif

if MAXFREQ ne MINFREQ then begin
	if FREQSAMPLINGRATE gt (MAXFREQ - MINFREQ)   then  begin
		print,strcompress('FreqSamplingRate='+string(FREQSAMPLINGRATE)+' too big, using default = 1.')
		FREQSAMPLINGRATE = 1 ; if too big then use default
	endif
endif else FREQSAMPLINGRATE = 1 ; if there is only one frequency

spe_nelements = floor((MAXFREQ - MINFREQ )/FREQSAMPLINGRATE)+1



if keyword_set(AMPLITUDE) and keyword_set(POWER) then begin
	print,'Invalid Keyword Pattern! Defaulting to Local Amplitude Spectra calculation'
	Power = 0
endif

;********************************
; Calculate the ST of the data
;********************************
; make sure that factor is an array of approriate size
if n_elements(factor) eq 1 then begin
	factor = dblarr(spe_nelements)+factor
endif
; make sure array passes is of right length
if n_elements(factor) ne spe_nelements then begin
	print,strcompress('Invalid Factor Array. Length is: '+string(n_elements(factor))+' and should be equal to '+string(spe_nelements))
	print,strcompress('!!! Replacing factor array with its first value !!!')
	factor = dblarr(spe_nelements)+factor(0)
endif
; make sure array does not have a zero
wfactorzero = where(factor eq 0, wfactorzerocount)
if wfactorzerocount gt 0 then begin
	Message,strcompress('Invalid Factor Array. Factor has a value of zero!')
endif




; Calculate the POSITIVE frequencies from minfreq to maxfreq
;********************************
if keyword_set(AMPLITUDE) or keyword_set(POWER) then begin
    loc = dblarr(length,spe_nelements)
	if keyword_set(AMPLITUDE)  then if keyword_set(VERBOSE) then print,'Calculating Local Amplitude Spectra.'  $
		else if keyword_set(VERBOSE) then print,'Calculating Local Power Spectra.'
	h = shift(h,-MINFREQ)
	if MINFREQ eq 0 then  begin
		gw = dblarr(length)
		gw(0) = 1
		loc(*,0) = abs(fft(h*gw,1))
	endif else begin
		f = double(MINFREQ)
    	width = factor(0) * length/f
		gw = gaussian_window(length,width)
		b = h * gw
    	loc(*,0) = abs(fft(b,1))
	endelse
	for index = 1d,spe_nelements-1 do begin
 		f = double(MINFREQ) + index*FREQSAMPLINGRATE
    	width = factor(index) * length/f
		gw = gaussian_window(length,width)
		h = shift(h,-FREQSAMPLINGRATE)
    	b = h * gw
    	loc(*,index) = abs(fft(b,1))
 	endfor
	if keyword_set(POWER) then loc = loc^2
endif else begin  ; calculate complex ST
	if keyword_set(VERBOSE) then print,'Calculating Local Complex Spectra'
	loc = dcomplexarr(length,spe_nelements)
	h = shift(h,-MINFREQ)
	if MINFREQ eq 0 then begin
		gw = dblarr(length)
		gw(0) = 1
		loc(*,0) = fft(h*gw,1)     ; 0 freq. equal to DC level
	endif else begin
		f = double(MINFREQ)
    	width = factor(0) * length/f
		gw = gaussian_window(length,width)
		b = h * gw
    	loc(*,0) = fft(b,1)
	endelse
 	for index = 1d,spe_nelements-1  do begin
 		f = float(MINFREQ) + index*FREQSAMPLINGRATE
		width = factor(index) * length/f
    	gw = gaussian_window(length,width)
	    h = shift(h,-FREQSAMPLINGRATE)
    	b = h * gw
    	loc(*,index) = fft(b,1)
    endfor
endelse


; now calculate negative frequencies from  - minfreq - maxfreq
;********************************
; NEGATIVE FREQUENCIES  -  calc the ST again, this time on the neg freqs
;********************************
if keyword_set(NEGATIVE_FREQUENCIES) then begin
	;create a new spectrum
	h = fft(time_series,-1)
	h = reverse(h); reverse spectrum so negative frequencies are first
	h = shift(h,1) ; shift DC to index 0
	if MINFREQ eq 0 then begin
			posminfreq = MINFREQ  ; save original value if its needed for frequency calcs
			MINFREQ = 1 ; step past Dc and start at first (negative) harmonic
	endif
	ori_spe_nelements = spe_nelements
	spe_nelements = floor((MAXFREQ - MINFREQ )/FREQSAMPLINGRATE)+1 ; recalc spe_nelements (in case minfreq changed)
	if keyword_set(AMPLITUDE) or keyword_set(POWER) then begin
    	loc_n = dblarr(length,spe_nelements)
		if keyword_set(AMPLITUDE)  then if keyword_set(VERBOSE) then begin
			print,'Calculating Local Amplitude Spectra (negative frequencies).'  
		endif else begin
			if keyword_set(VERBOSE) then print,'Calculating Local Power Spectra (negative frequencies).'
		endelse
		h = shift(h,-MINFREQ)
		f = double(MINFREQ)
   		width = factor(0) * length/f
		gw = gaussian_window(length,width)
		b = h * gw
    	loc_n(*,0) = abs(fft(b,1))
		for index = 1d,spe_nelements-1 do begin
 			f = double(MINFREQ) + index*FREQSAMPLINGRATE
	    	width = factor(index) * length/f
			gw = gaussian_window(length,width)
			h = shift(h,-FREQSAMPLINGRATE)
	    	b = h * gw
    		loc_n(*,index) = abs(fft(b,1))
 		endfor
		if keyword_set(POWER) then loc_n = loc_n^2
	endif else begin  ; calculate complex ST
		if keyword_set(VERBOSE) then print,'Calculating Local Complex Spectra (negative frequencies)'
		loc_n = dcomplexarr(length,spe_nelements)
		h = shift(h,-MINFREQ)
		f = double(MINFREQ)
    	width = factor(0) * length/f
		gw  = gaussian_window(length,width)
		b = h * gw
    	loc_n(*,0) = fft(b,1)
 		for index = 1d,spe_nelements-1  do begin
			f = double(MINFREQ) + index*FREQSAMPLINGRATE
 			width = factor(index) * length/f
    		gw = gaussian_window(length,width)
	    	h = shift(h,-FREQSAMPLINGRATE)
    		b = h * gw
    		loc_n(*,index) = fft(b,1)
    	endfor
	endelse
	; if sampled negative nyquist freq component (and n mod 2 = 0) then drop that redundant voice
	if maxfreq eq spe_length then begin
		if (length+1) mod 2 then begin
			loc_n = loc_n(*,0:spe_nelements-2)
		endif
	endif

	loc_n=(reverse(loc_n,2)) ; flip columns to normal FFT frequency indexing (0, pos, nyq ; neg nyqu, neg)
	loc_n=(reverse(loc_n,1)) ; undo previous reverse() which causes time reversal
	loc_n=(shift(loc_n,1))	 ; ditto
	loc = transpose([transpose(loc),transpose(loc_n)])
	loc_n = 0
endif

nfreq_elements = (size(loc))(2)
if keyword_set(VERBOSE) then print,strcompress('The number of frequency elements is'+string(nfreq_elements)+'.')


if keyword_set(maskedges) then begin
	if maskedges eq 1 then maskthreshold=0.05
	if maskedges gt 0 and maskedges lt 1 then maskthreshold=maskedges
	if maskedges gt 1 and maskedges le 100 then maskthreshold=float(maskedges)/100. $
		else  maskthreshold=0.05
	edgets = dindgen(length)/length
	st = s_trans(edgets,/abs)
	mask=where(st gt maskthreshold,maskcount) ; 5 % is good = 0.05 based on snooping around
	loc(mask) = 0
endif

if keyword_set(SAMPLINGRATE) then begin  ; make structure
	; MINFREQ was changed in negative frequency calculations to avoid the DC
	if n_elements(posminfreq) gt 0 then MINFREQ = posminfreq 
	if n_elements(ori_spe_nelements) gt 0 then  spe_nelements = ori_spe_nelements
	frequencies = (MINFREQ + dindgen(spe_nelements)*FREQSAMPLINGRATE)/(SAMPLINGRATE*length)
	negfreqs =  -reverse(frequencies)
	if keyword_set(NEGATIVE_FREQUENCIES) then begin ; get rid of reduntaly sampled nyquist point
		if maxfreq eq spe_length then begin
			if (length+1) mod 2 then begin
				negfreqs = negfreqs(1:*)
			endif
		endif
		; drop the DC component from negative freqs
		if negfreqs(N_elements(negfreqs)-1) eq 0 then negfreqs = negfreqs(0:n_elements(negfreqs)-2)
		frequencies = [frequencies,negfreqs]
	endif
	time = dindgen(length)*SAMPLINGRATE
	a = {st: loc, time: time,freq: frequencies}
	return,a
endif else return, loc



end ; end of function

