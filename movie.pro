nx=400
ny=400
lx=20.
ly=20.
dx=lx/nx & dy=ly/ny
x=dx*indgen(400)
y=dy*indgen(400)
t=indgen(40)*0.5
close,1
close,2
close,3
close,4
close,5

for episode =10,200,10 do begin
;episode=1
n=string(episode,format='(i4.4)')
a=episode

filename='omg.gif'

;restore,file='~/ct_uz.sav'
loadct,0
zr=[-1.,1]*4
nl=21
lev=zr[0]+(zr[1]-zr[0])*indgen(nl)/nl & lev=[-100,lev,zr[1]]
c_col=0+254*indgen(nl)/(nl-1) & c_col=[0,c_col,254]
;cbar=[[lev],[lev]]

DEVICE,Decomposed=0
window,1,xs=600,ys=600
openr,1,n+'pos_movie.bin'
;openr,2,'position_y_movie.bin'
openr,3,n+'angle_movie.bin'
;openr,4,'angle_2.bin'
openr,5,n+'omg_movie.bin'
hpx=0.d0 & hpy=0.d0 & ha=0.d0 & hb=3.14d0
w=dblarr(nx,ny)

for n=0,n_elements(t)-1 do begin
readu,1,hpx
;readu,2,hpy
readu,3,ha
;readu,4,hb
readu,5,w

hpy = 10
w=shift(w,(hpx-10)/dx,(hpy-10)/dy)

contour,w,x,y,lev=lev,/fi,/iso,c_col=c_col,col=255,back=128
;plot,indgen(21),indgen(21),/nodata
oplot,[hpx+20*((tanh(hpx*10^10)+1)/2-fix(hpx/20))],[hpy],ps=4
x0=cos(ha) & y0=sin(ha)
x1=cos(hb) & y1=sin(hb)
plots,[hpx+20*((tanh(hpx*10^10)+1)/2-fix(hpx/20))+x1,hpx+20*((tanh(hpx*10^10)+1)/2-fix(hpx/20))],[hpy+y1,hpy],thick=2
plots,[hpx+20*((tanh(hpx*10^10)+1)/2-fix(hpx/20))+x0,hpx+20*((tanh(hpx*10^10)+1)/2-fix(hpx/20))],[hpy+y0,hpy],thick=2
xyouts,1,1,'t='+string(t[n],format='(f6.2)'),chars=2,col=255,chart=2
xyouts,1,3,'episode='+string(a,format='(i3.3)'),chars=2,col=255,chart=2
write_gif,filename,tvrd(),cr,cg,cb,/multi,REPEAT_COUNT=0, DELAY_TIME=10
;print,t[n]
endfor
write_gif,filename,/close
close,1
;close,2
close,3
;close,4
close,5
endfor

end
