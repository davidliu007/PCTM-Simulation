;*************************************************************************
; This program aims to provide a quick and convenient way to check and   *
; verify the concerned total global flux surface mass.                   *
;========================================================================*
; Arguments as inputs
;
; Arguments as outputs
;
; NOTES:
;
; Revision History
;
;*************************************************************************
@ /home/ypliu/tools/idl/cvs/idl4pctm/cellarea_multigrid.pro
;GOTO,plt
;GOTO,outputs
;
;========================================================================*
;goto,pp
;Parameters for specific run
      yr=2016 & syr=strtrim(yr,2)
      ydays=365
;      trc='fuel'
;      trc='fire_emiss_daily_'
;      trc='fire_emiss_3hrly'
;      trc='woodfuel_n_wildfire_3hrly'
;      trc='m2NEE'
      trc='NPP'
      trc='resp'
;      trc='flux.e632.co2.3hrly_576x361_' & res='H'
;      trc='GPP' & res='H'
;      trc='flux_newcasa_3hrly_288x181_'  & res='L'

;      dir0='/misc/srk_raid05/ypliu/CASA_GFED3_2014plus/'
;      dir0='/misc/srk_raid05/ypliu/base_fossil_ocn/'
      dir0='/misc/srk_raid05/ypliu/m2CO2flux/'
;      dir0='/misc/srk_raid05/ypliu/casa_flux_cntrl/'

      fname0=dir0+trc+syr+'.txt'
;      fname0=dir0+trc+syr+'merra2.txt' ;For raw daily fire
;      fname1=dir0+trc+syr    ; for flux

      fname1=dir0+trc+syr+'_576x361.dat' 
;      fname1=dir0+trc+syr+'_1.25x1_e614.dat'
;      fname1=dir0+trc+syr+'_3hrly_576x361.dat' 
;      fname1=dir0+trc+'_288x181.'+syr ;wood fuel and wild fire
;      fname1=dir0+trc+syr+'_3hrly_1.25x1_e614.dat'

;      fname1=dir0+trc+syr+'_3hrly_1.25x1_e629.dat'
;      fname1='/misc/srk_raid05/ypliu/CASA_GFED2012_e614/NEP2012_3hrly_1x1.25_180.dat'
;
;      fname0=dir0+'fire_emiss_daily_2012.txt'
;      fname1=dir0+'fire_emiss_daily_2009_1x1.25_180.dat'
;      fname1=dir0+'Taka09.1x1.25.3hr.data'  ;unit: kg/m2/sec
;      fname1=dir0+trc+syr        ;for flux, unit: kg/m2/sec
;       fname1=dir0+'Taka02.1x1.25.3hr.data' & res='L'
;       fname1=dir0+'Taka09.576x361.3hr.data' & res='H'
;       fname1=dir0+'Taka02.576x361.3hr.data' & res='H'
;       fname1=dir0+'taka09_kgcm2s_576x361.dat' & res='H'
;       fname1=dir0+'Taka02.576x361.data' & res='H'
;       fname1=dir0+'D_float.2005' & res='L'
;       fname1=dir0+'D_float576x361.2005' & res='H'

      print,syr, ' ',trc ,' Month 1-12'
      print, '------------------------------------'
       print,'Raw data file name',fname0
       print,'Regrided data file name',fname1
      dmon=[0,31,28,31,30,31,30,31,31,30,31,30,31]
      if((yr mod 4) eq 0) then begin
         dmon[2]=29
         ydays=366
      endif
;
      output_ps=0         ;flag of output format
      raw=0 & tinc0='m'   ;flag of if process raw casa data and if it is daily data 
      npd=1

      if(raw eq 1) then goto,raw0

;      nflx=1                        ;for one layer data
      tinc='m' & nflx=12            ;for monthly data
;      tinc='d' & npd=8 & nflx=ydays*npd           ;for daily data and 3hourly data
      res='H'

      if(res eq 'L') then begin
        xdim=288
        ydim=181
      endif else begin
        xdim=576      
        ydim=361
      endelse
;      inc=long(xdim)*long(ydim)*4L       ;increase of each map size
;      loc=inc*(ndoff*long(npd)+td)       ;data read in start point      
;------------------------------------------------------------------------
      dx = 360.0/xdim
      dy = 180.0/(ydim-1)
      longitudes = findgen(xdim)*dx      ;for [0:360] D_float.yyyy
      longitudes = findgen(xdim)*dx-180  ;for [-180:180] in longitude
      latitudes = findgen(ydim)*dy-90.0
;
;Open and read regrided file 
;
      flux1=fltarr(xdim,ydim,nflx)

      openr, lun, fname1, /get_lun;,/swap_endian  ;for ocean, fossil, nee and flux
      readu, lun, flux1

      free_lun,lun
;      stop,'mmmm'
;     goto,plt
goto,outputs
pp:
for m=3,8 do begin
ncols=40 & zmin=min(flux1[*,*,m])  &  zmax=max(flux1[*,*,m])  
step=(zmax-zmin)/ncols
levs=Indgen(ncols)*step+zmin

;Setting for screen outputs, good for quick check.
ctit=fname1
device,decomposed=0
window,/free,title=ctit+string(m+1,format='(I2.2)')
loadct,33

print,fname1
lat0=0  &  lon0=0  &  rot=0
cdat1=Wrap(flux1(*,*,m),longitudes,lon1=lon1)     ;for monthly display
lon1(xdim)=180.000  ;for [-180:180]
Map_set,lat0,lon0,rot,/nob,/adv,tit=ctit,col=10,ymarg=[3,3],xmarg=[3,3],/mer;,scale=35e6    
Contour,cdat1,lon1,latitudes,lev=levs,/fill, /over
Map_continents  &  Map_grid,/label,color=150
endfor

      stop,'mmmm'

outputs:
;cellarea288x181 is in cellarea_multigrid.pro, it also works well with resolution of 576x361.
;     cellarea=cellarea288x181(longitudes,latitudes) 
     cellarea=cellarea_ctr2(longitudes,latitudes) 
;      print, '------------------------------------'
;      print,'total regrided earth surface area: ',total(cellarea)
      dmon=[0,31,28,31,30,31,30,31,31,30,31,30,31]
      if((yr mod 4) eq 0) then dmon[2]=29
;      sum1=0.0D
      sum1=[]
      sum0=[]
      toff=0
      flx_mon=fltarr(xdim,ydim,12)

      if( tinc eq 'm') then begin
        for m=0,nflx-1 do begin ; month loop
	  subm=0.0D
          for i=0,ydim-1 do begin
            for j=0,xdim-1 do begin
              subm=flux1(j,i,m)*cellarea(j,i)+subm
            endfor
          endfor
;          print,'regrided surface flux sum(gC/Month)=: ',subm
	  sum1=[sum1,subm]
        endfor
;	print,'regrided surface flux sum(gC/yr)=: ',sum1
      endif else begin
        for m=0,11 do begin ; month loop
        subm=0.0
        toff=toff+dmon(m)*npd
;        print,m+1, toff, toff+dmon(m+1)*npd
        for nt=toff,toff+dmon(m+1)*npd-1 do begin
;          for i=0,ydim-1 do begin
;            for j=0,xdim-1 do begin
              flx_mon[*,*,m]=flx_mon[*,*,m]+flux1(*,*,nt)*cellarea
              sumhr1=total(flux1(*,*,nt)*cellarea)
              subm=sumhr1+subm
	      sum0=[sum0,sumhr1]
;            endfor
;          endfor
        endfor
        print,'3hrly regrided surface flux sum(gC/Month)=: ',subm*3*3600*1000
;        print,'monthly regrided surface flux sum(gC/Month)=: ',subm
;        sum1=subm+sum1
        sum1=[sum1,subm]
        endfor ; end month loop
      print,'3 hrly regrided surface flux sum(gC/yr)=: ',total(sum1)*3*3600*1000
;      print,'Annual regrided surface flux sum(gC/yr)=: ',sum1
      sum0=sum0*3*3600*1000
      sum1=sum1*3*3600*1000
      flx_mon=flx_mon*3*3600*1000
      endelse

goto,raw0
;stop,'look up'
;Plot regrided data
plt:

for m=221*8,221*8+3 do begin
ncols=40 & zmin=min(flux1[*,*,m])  &  zmax=max(flux1[*,*,m])  
;ncols=40 & zmin=min(flx_mon[*,*,m])  &  zmax=max(flx_mon[*,*,m])  
step=(zmax-zmin)/ncols
levs=Indgen(ncols)*step+zmin

;Setting for screen outputs, good for quick check.
ctit=fname1
device,decomposed=0
window,/free,title=ctit+string(m+1,format='(I2.2)')
loadct,33

print,fname1
lat0=0  &  lon0=0  &  rot=0
cdat1=Wrap(flux1(*,*,m),longitudes,lon1=lon1)     ;for monthly display
;cdat1=Wrap(flx_mon(*,*,m),longitudes,lon1=lon1)     ;for monthly display
lon1(xdim)=180.000  ;for [-180:180]
Map_set,lat0,lon0,rot,/nob,/adv,tit=ctit,col=10,ymarg=[3,3],xmarg=[3,3],/mer;,scale=35e6    
Contour,cdat1,lon1,latitudes,lev=levs,/fill, /over
Map_continents  &  Map_grid,/label,color=150
endfor
stop,'Check regrided Results'
;========================================================================*
;Following processes the raw gridcell (720latx360lon)
raw0:
raw=1
sumday=[]
    IF(raw gt 0) THEN BEGIN

;      FORM='(720(1x,f8.3))'   ;MOD15 NPP, resp and Fuel
      FORM='(720(1x,f9.3))'   ;e614 NPP, resp and Fuel
;      FORM='(720(1x,f9.4))'   ;e614 NPP, resp and Fuel
      nflx0=12

      IF(tinc0 eq 'd') THEN BEGIN
        nflx0=ydays
        FORM='(720(1x,f8.3))'    ;MOD15 Daily Fire
      ENDIF

      ydim0=360
      xdim0=720
      dx0 = 360.0/xdim0
      dy0 = 180.0/ydim0
;      longitudes0 = findgen(xdim0)*dx0      ;for [0:360]
      longitudes0 = findgen(xdim0)*dx0-180  ;for [-180:180] in longitude
      latitudes0 = findgen(ydim0)*dy0-90.0

      radius=6.37122e6
      pi=3.14159

      celledge0=fltarr(ydim0+1)
      cellarea0=fltarr(xdim0,ydim0)
      dphi0=180./ydim0                 ;increment at lat direction
      celledge0(0)=-90.0                   ;edge at south pole
      for j=1,ydim0 do celledge0(j)=celledge0(j-1)+dphi0
      celledge0=celledge0*pi/180.0               ;convert to rad from degree
;compute areas of each cellnee_3hrly_m2gpp_nep4lfeng.driver
      fact0=2.*pi*radius*radius/xdim0    ;half global surface divided into number 'long' of equal parts
      for j=1,ydim0 do begin
        bandarea=(sin(celledge0(j))-sin(celledge0(j-1)))*fact0
        for i=0,xdim0-1 do cellarea0(i,j-1)=bandarea        ;for east-west band all cells have same area2.18599e+15
      endfor
;
;Read in raw 720x360 data
      flux0=fltarr(xdim0,ydim0,nflx0)
      flux2=fltarr(xdim0,ydim0,12)
      tmp0=fltarr(xdim0)
      openr, lun, fname0, /get_lun
;      point_lun,lun,loc
      n=0L
      WHILE NOT EOF(lun) DO BEGIN
        i=n/ydim0
	j=n mod ydim0
;        readf,lun,tmp0, form=FORM    ;for mod15, mod18
        readf,lun,tmp0                ;for e614 raw NPP and daily fire files
        flux0(*,ydim0-1-j,i)=tmp0
        n=n+1
      ENDWHILE

      free_lun,lun
;goto,rplt      
;
;compute 720x360 resolution data
;
;goto,day
;      print, '------------------------------------'
;      print,'total raw earth surface area: ',total(cellarea0)
      sum0=0
      for nt=0,nflx0-1 do begin
        summ=0.0
        for i=0,ydim0-1 do begin
          for j=0,xdim0-1 do begin
            summ=flux0(j,i,nt)*cellarea0(j,i)+summ
          endfor
        endfor
        sumday=[sumday,summ]
        sum0=summ+sum0
        endfor
;      print,'raw surface flux sum(gC/yr)=: ',sum0
      print, '------------------------------------'
      print, 'MON     RAW       REGRIDED      DIFFERENCE[raw-reg] PERCENTAGE[(raw-reg)x100/raw]'
      for k=0,11 do begin
        print,string(k+1,'(I2.2)'),sumday[k], sum1[k], sumday[k]-sum1[k],(sumday[k]-sum1[k])/sumday[k]*100.
      endfor
      print, 'Year Sum up:'
      print,total(sumday), total(sum1), total(sumday-sum1),total(sumday-sum1)/total(sumday)*100.

stop,'MMM'
day:
      print,'print monthly wild fire sum:'
      sum0=[]
      sum1=[]
      toff=0
      for m=0,11 do begin ; month loop
        subm=0.0
        tmp=fltarr(xdim0,ydim0)
        toff=toff+dmon(m)*npd
        print,m+1, toff, toff+dmon(m+1)*npd
        for nt=toff,toff+dmon(m+1)*npd-1 do begin
	  sumhr=total(flux0[*,*,nt]*cellarea0)
	  subm=sumhr+subm
	  sum0=[sum0,sumhr]
	  tmp=flux0[*,*,nt]+tmp
        endfor
;        print,'monthly regrided surface flux sum(gC/Month)=: ',subm
;       sum1=subm+sum1
        sum1=[sum1,subm]
	flux2[*,*,m]=tmp/(dmon[m+1]*npd)
        endfor ; end month loop


    ENDIF  ;End raw data processing

stop,'end raw fire computing'
;Plot raw data
rplt:
m=6
ncols=40 & zmin=min(flux0[*,*,m])  &  zmax=max(flux0[*,*,m])  
step=(zmax-zmin)/ncols
levs=Indgen(ncols)*step+zmin

;Setting for screen outputs, good for quick check.
device,decomposed=0
window,/free,title=ctit+string(m+1,format='(I2.2)')
loadct,33

ctit=fname0
lat0=0  &  lon0=0  &  rot=0
cdat0=Wrap(flux0(*,*,m),longitudes0,lon1=lon10)     ;for monthly display
lon10(xdim0)=180.000  ;for [-180:180]
Map_set,lat0,lon0,rot,/nob,/adv,tit=ctit,col=10,ymarg=[3,3],xmarg=[3,3],/mer;,scale=35e6    
Contour,cdat0,lon10,latitudes0,lev=levs,/fill, /over
Map_continents  &  Map_grid,/label,color=150

stop,'*** Check raw Results ***'

;#######################################################################
plts:
  
      rg_plt=1
      ndoff=195L                           ;days offset/1st day from begining of the year
      td=6L                              ;begining time of the day (0,1,2,3,4,5,6,7)
      days=2                             ;days of period to be checked,for mode='m'
      mode='s'                           ;single time plot or multiple time plot
;      m=ndoff*long(npd)+td               ;records offset
      m=7

      ncols=40 & zmin=min(flux1[*,*,m])  &  zmax=max(flux1[*,*,m]*2.)  
      ncols=40 & zmin=-4.0e-9  &  zmax=7.0e-9
      step=(zmax-zmin)/ncols
      levs=Indgen(ncols)*step+zmin
;      fac=1.0;1.e8
      fac=1.e9

ctit=syr+' Taka02 Aug. 576x361 CO2 flux plot (*e9,kgC/m2/month)'
;Setting for ps file output(see corresponding color bar below the contour) ******
if(output_ps eq 1) then begin
tvpinit,/col   &  tgcoltab,'spectrum3',ncols
!p.multi=[0,2,2]
endif else begin
;Setting for screen outputs, good for quick check.
device,decomposed=0
window,/free,title=ctit
loadct,33
endelse

tmpflx=fltarr(288,361,12)
tmpflx=flux1(0:287,*,*)
flux1(0:287,*,*)=flux1(288:575,*,*)
flux1(288:575,*,*)=tmpflx

;npp=(resp-npp)*1.e-3/(30.*86400.)
lat0=0  &  lon0=0  &  rot=0
;lat0=45.  &  lon0=-90  &  rot=0

;#for i =0,11 do begin           ;for monthly display
      If(rg_plt eq 1) then begin
;Plot regrided data
        cdat=Wrap(flux1(*,*,m),longitudes,lon1=lon1)     ;for monthly display
;        cdat=Wrap(npp(*,*),longitudes,lon1=lon1)
;        lon1(xdim)=360.000  ;for [0:360]      ; For D_float.yyyy
        lon1(xdim)=180.000  ;for [-180:180]
        Map_set,lat0,lon0,rot,/nob,/adv,tit=ctit,col=10,ymarg=[3,3],xmarg=[3,3],/mer;,scale=35e6
        Contour,cdat,lon1,latitudes,lev=levs,/fill, /over;,$  ; for ocea snf permanent frost
        Map_continents  &  Map_grid,/label,color=150
      endif else begin
;Plot raw data
        cdat0=Wrap(flux0(*,*,m),longitudes0,lon1=lon10)     ;for monthly display
        lon10(xdim0)=180.000  ;for [-180:180]
        Map_set,lat0,lon0,rot,/nob,/adv,tit=ctit,col=10,ymarg=[3,3],xmarg=[3,3],/mer;,scale=35e6    
        Contour,cdat0*fac0,lon10,latitudes0,lev=levs,/fill, /over;,$  ; for ocea snf permanent frost
        Map_continents  &  Map_grid,/label,color=150
      endelse

  if(output_ps eq 1) then begin
    tgcolbar,'ppmv',[levs(0),levs(N_elements(levs)-1)+npc],[220,0,200,15],nt=8,$
    form='(i3)',/labels,font=1
  endif else begin
;    marks=['0','50','100','150','200','250','300']
;    marks=['0','4','8','12','16','20']
;    div=5
    COLORBAR, div=6,POSITION=[0.30, 0.10, 0.70, 0.14],form='(f9.2)',max=zmax*fac, $
    minr=zmin*fac, chars=0.8, col=150
  endelse

      stop,'check'
END
