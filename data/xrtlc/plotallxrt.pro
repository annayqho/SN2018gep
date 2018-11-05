function readxrtlc, filename

  if file_test(filename) eq 0 then begin
     print, 'Cannot find ', filename
     return, -1
  endif

  n = countlines(filename)

  lc = replicate({t:0., f:0., dt:0., ef:0.}, n)
  
  openr, 1, filename
  inline = ''
  i = 0
  for l = 0, n-1 do begin
    readf, 1, inline
    if strlen(inline) lt 20 then continue
    inline = clip(repstr(inline,'!#',''))
    char1 = strmid(inline,0,1)
    if (char1 lt '0' or char1 gt '9') and char1 ne '-' then continue
    indata = strsplit(inline, /extract)
    if n_elements(indata) lt 5 then continue
    lc[i].t = indata[0] 
    ;if strmid(clip(indata[2]),0,1) eq '*' then  print, filename, inline else $
    lc[i].dt = abs(indata[1]) + abs(indata[2])
    lc[i].f = indata[3]
    lc[i].ef = indata[4]
    i += 1
  endfor

  close, 1
  lc = lc[0:i-1]
  return, lc
end


pro plotallxrt

   za = 0.338
   zeff = 1.0

   psopen, 'allxrt.eps', xsize=11, ysize=5.5, /inches, /encaps
   !p.font = 0
   device, /helvetica, font_index = 17
   device, /color
   colors = transpose([[0,0,0], [200,0,0], [0,125,0], [0,0,200], $
                                [255,125,125], [125,255,125], [155,155,255], $
                                [200,200,200]])
   tvlct, colors
   !p.multi = [0,2,1]

; 060124, 100902A, 110801A - X-ray flares
; 060218, 111209A          - mid-time due to lonhg lasting activity
; 100615                   - late-time due to +perr and short exposure
; 100621, 060729           - late-time
; 090417B                  - flare, late time

   readcol, 'grbz.txt', grb, grbz, format='a,f', comment='#', /silent
   grb = reverse(grb)
   grbz = reverse(grbz)
   w = where(grb eq '090417B' or grb eq '130427A')
   grbz = grbz[w]
   grb = grb[w]
   ;grbz = grbz[where(grb eq '060124' or grb eq '060218' or grb eq '100902A' or grb eq '110801A' or grb eq '111209A'  or grb eq '130427A')]
   ;grb = grb[where(grb eq '060124' or grb eq '060218' or  grb eq '100902A' or grb eq '110801A' or grb eq '111209A'  or grb eq '130427A')]
   color = 1+indgen(n_elements(grb)) mod 6
   
   good = intarr(n_elements(grb))
   for i = 0, n_elements(grb)-1 do begin
      filename = './'+strlowcase(grb[i])+'_xrt_bin.txt'
      if file_test(filename) eq 0 then continue
      lc = readxrtlc(filename)
      if n_elements(lc.t) lt 3 then continue
      t0 = min(lc.t)
      if t0 gt 1000. then continue
      good[i] = 1
   endfor

   grb = grb[where(good)]
   grbz = grbz[where(good)]
   color[where(clip(grb) eq '130427A')] = 0


   !p.position = [0.055, 0.105, 0.45,0.975]

   xrange = [3e-4, 100]
   yrange = [1e-14, 1e-5]
   plot, [0],[0], /xlog, /ylog, xrange=xrange, yrange=yrange, /xsty, ysty=1, title='!17', xtitle='t (days)', ytitle='XRT flux (0.3-10 keV)', xtickformat='expaxis'


   for i = 0, n_elements(grb)-1 do begin
      filename = './'+strlowcase(grb[i])+'_xrt_bin.txt'
      lc = readxrtlc(filename)
      if clip(grb[i]) eq '130427A' then thick=5 else thick=1
      oplot, lc.t/(24.*3600), lc.f, color=color[i], thick=thick
   endfor


   !p.position = [0.555, 0.105, 0.95,0.97]
   plot, [0],[0], /xlog, /ylog, xrange=xrange, yrange=yrange, xsty=1+8, ysty=1+8, $
         title='!17', xtitle='!17t (days, observer frame if z=1)', ytitle='XRT flux (0.3-10 keV), z='+clipzero(zeff), xtickformat='expaxis'
   lumrange=yrange * 4 * !pi * (3.0857d24*lumdist(zeff,/s))^2 / (1.+zeff)   ; last term is the default k correction for f_nu
   axis, yaxis=1, yrange=lumrange, /ysty, /ylog, ytitle='0.6-20 keV Luminosity (erg/s)'
   axis, xaxis=1, /xlog, xrange=xrange/2., xstyle=1, xtitle='t (days, rest frame)', xtickformat='expaxis'

   for i = 0, n_elements(grb)-1 do begin
      lc = readxrtlc('./'+strlowcase(grb[i])+'_xrt_bin.txt')
      z = grbz[i]
      if clip(grb[i]) eq '130427A' then thick=5 else thick=1
      kredshift, 50., lc.f, z, 50., fzeff, zeff, index=-1. ; note using a uniform index
      ;print, clip(grb[i],7), clip(grb[i]) eq '130427A', grbz[i], fzeff[0]/lc[0].f, color[i], thick
      oplot, lc.t*(1.+zeff)/(1.+z)/(24.*3600.), fzeff>1e-20, color=color[i], thick=thick

      w = where(lc.t/(24.*3600.) gt 0.5 and lc.f*lc.t/(24.*3600.) gt 1d-11*1., ct)
      if ct gt 1 then print, grb[i]

   endfor


   psclose


end

