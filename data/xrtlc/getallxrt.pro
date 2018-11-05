pro grabxrtlc, grb
   trigger = getbattrigger(grb)
   outfile = './'+strlowcase(grb)+'_xrt.txt'
   if file_test(outfile) eq 0 then begin
      cmd = 'wget http://www.swift.ac.uk/burst_analyser/' + string(long(trigger),'(I08)') + '/xrt/xrt_flux_XRTBAND.qdp.gz -O '+outfile
      spawn, cmd
   endif
   binxrt, './'+strlowcase(grb)+'_xrt.txt',   './'+strlowcase(grb)+'_xrt_bin.txt'
end


pro getallxrt
   readcol, '~/progs/idl/research/triggerindex.txt', grb, trig, comment='#', format='a,a'
   for i = 0, n_elements(grb)-1 do begin
       print, grb[i]
       grabxrtlc, grb[i]
   endfor

end
