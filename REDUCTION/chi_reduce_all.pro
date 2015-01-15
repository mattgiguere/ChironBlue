;+
;
;  NAME: 
;     chi_reduce_all
;
;  PURPOSE: 
;   A procedure designed to run Andrei's sorting hat routine for all modes
;
;  CATEGORY:
;      CHIRON
;
;  CALLING SEQUENCE:
;
;      chi_reduce_all
;
;  KEYWORD PARAMETERS:
;    
;  EXAMPLE:
;      chi_reduce_all
;
;  MODIFICATION HISTORY:
;        c. Matt Giguere 2011.11.18 11:05:47 AM
;
;-
pro chi_reduce_all, $
help = help, $
date = date, $
doppler = doppler, $
skipbary = skipbary, $
skipdistrib = skipdistrib, $
skipqc = skipqc

print, 'Started @: ', systime()
spawn, 'hostname', host

rawdir = '/raw/mir7/'
lfn = '/tous/mir7/logsheets/20'+strmid(date, 0, 2)+'/'+strt(date)+'.log'

;This part gets the image prefix:
spawn, 'ls -1 '+rawdir+date+'/', filearr
nel = n_elements(filearr)
nfa = strarr(nel)
for i=0, nel-1 do nfa[i] = strmid(filearr[i], 0, strlen(filearr[i])-9)
uniqprefs =  nfa(uniq(nfa))
pref = 'junk'
ii=-1
repeat begin
  ii++
  pref = uniqprefs[ii]
  print, 'pref is: ', uniqprefs[ii]
  print, 'pref 1st 2 are: ', strmid(uniqprefs[ii], 0,2)
endrep until ( ((strmid(uniqprefs[ii],0,2) eq 'qa') or $
        (strmid(uniqprefs[ii],0,2) eq 'ch')) and $
        (strmid(uniqprefs[ii], 0, 4) ne 'chir') and $
        (strmid(uniqprefs[ii], 0, 1, /reverse) eq '.'))

print, 'pref is: ', pref
;stop, 'pref is: ', pref

modearr = [$
'narrow', $ 
'slicer', $
'slit', $
'fiber']

for i=0, 3 do begin
  print, '*************************************************'
  print, ' NOW ON TO THE ', MODEARR[I], ' MODE...'
  print, '*************************************************'
  sorting_hat,date,run=pref,mode=modearr[i],/reduce,/getthid,/iod2fits
endfor

print, 'Finished @: ', systime()

end;chi_reduce_all.pro