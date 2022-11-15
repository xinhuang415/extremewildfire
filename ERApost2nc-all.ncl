load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
begin
;oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooodefination

ny=180
nx=720
yy=2020
lat1=ispan(0,ny-1,1)-90+0.5
lon1=ispan(0,nx-1,1)-180+0.5

do yy0=2000,2000
fils:= systemfunc("ls ERA5/hourly_surface/ERA5_hourly_basic_meteo_para_"+yy0+"-*.nc") 
print(yy0+"	"+(fils))
f:= addfiles(fils+"","r")
T2tmp:=short2flt(f[:]->t2m(::3,:,:)) ;2m temperature
D2tmp:=short2flt(f[:]->d2m(::3,:,:)) ;2m dewpoint temperature
U10tmp:=short2flt(f[:]->u10(::3,:,:)) ;10m u-component of wind
V10tmp:=short2flt(f[:]->v10(::3,:,:)) ;10m v-component of wind

lat0=f[0]->latitude(:) 
lon0tmp=f[0]->longitude(:) 
lon0tmp=where(lon0tmp .ge. 180,lon0tmp-360,lon0tmp) ; longitude from -180 to 180

lon0=lon0tmp
lon0(0:719)=lon0tmp(720:1439)
lon0(720:1439)=lon0tmp(0:719)
T2:=T2tmp
T2(:,:,0:719)=T2tmp(:,:,720:1439)
T2(:,:,720:1439)=T2tmp(:,:,0:719)
D2:=T2tmp
D2(:,:,0:719)=D2tmp(:,:,720:1439)
D2(:,:,720:1439)=D2tmp(:,:,0:719)
U10:=T2tmp
U10(:,:,0:719)=U10tmp(:,:,720:1439)
U10(:,:,720:1439)=U10tmp(:,:,0:719)
V10:=T2tmp
V10(:,:,0:719)=V10tmp(:,:,720:1439)
V10(:,:,720:1439)=V10tmp(:,:,0:719)
print(lon0+"	"+lon0tmp)

dims=dimsizes(T2)
nt=dims(0)
print(dims)

if(isleapyear(yy0)) then
T2_daily:=new((/366,ny,nx/),float)
T2_daily!0="jday"
T2_daily!1="lat"
T2_daily!2="lon"
T2_daily&jday=ispan(1,366,1)
T2_daily&lat=lat1
T2_daily&lon=lon1
D2_daily:=T2_daily
U10_daily:=T2_daily
V10_daily:=T2_daily
T2_max:=T2_daily
D2_max:=T2_daily
U10_max:=T2_daily
V10_max:=T2_daily
T2_min:=T2_daily
D2_min:=T2_daily
U10_min:=T2_daily
V10_min:=T2_daily
else
T2_daily:=new((/365,ny,nx/),float)
T2_daily!0="jday"
T2_daily!1="lat"
T2_daily!2="lon"
T2_daily&jday=ispan(1,365,1)
T2_daily&lat=lat1
T2_daily&lon=lon1
D2_daily:=T2_daily
U10_daily:=T2_daily
V10_daily:=T2_daily
T2_max:=T2_daily
D2_max:=T2_daily
U10_max:=T2_daily
V10_max:=T2_daily
T2_min:=T2_daily
D2_min:=T2_daily
U10_min:=T2_daily
V10_min:=T2_daily
end if

do it=0,nt-1,8
id=toint(it/8.0)
print(id)
T2_daily(id,:,:)=area_hi2lores(lon0,lat0,dim_avg_n(T2(it:it+7,:,:),0),True,1,lon1,lat1,0) ;Interpolates from high resolution rectilinear grids to low resolution recitlinear grids using local area averaging.
D2_daily(id,:,:)=area_hi2lores(lon0,lat0,dim_avg_n(D2(it:it+7,:,:),0),True,1,lon1,lat1,0)  
U10_daily(id,:,:)=area_hi2lores(lon0,lat0,dim_avg_n(U10(it:it+7,:,:),0),True,1,lon1,lat1,0)  
V10_daily(id,:,:)=area_hi2lores(lon0,lat0,dim_avg_n(V10(it:it+7,:,:),0),True,1,lon1,lat1,0)
end do

setfileoption("nc", "Format",  "NetCDF4")
	system("rm -f erapost/erapost"+yy0+".nc")
fout=addfile("erapost/erapost"+yy0+".nc","c")
fout->T2_daily = T2_daily
fout->D2_daily = D2_daily
fout->U10_daily = U10_daily
fout->V10_daily = V10_daily
fout->lat = lat1
fout->lon = lon1
end do

fils:= systemfunc("ls ERA5/hourly_surface/ERA5_hourly_cloud_precip_"+yy0+"-*.nc") 
print(yy0+"	"+(fils))
f:= addfiles(fils+"","r")
tptmp:=short2flt(f[:]->tp(::3,:,:)) ;Total precipitation
cptmp:=short2flt(f[:]->cp(::3,:,:)) ;Convective precipitation
lsptmp:=short2flt(f[:]->lsp(::3,:,:)) ;Large-scale precipitation
tcctmp:=short2flt(f[:]->tcc(::3,:,:)) ;Total cloud cover
fils:= systemfunc("ls ERA5/hourly_surface/ERA5_hourly_ps_"+yy0+"-*.nc") 
print(yy0+"	"+(fils))
f:= addfiles(fils+"","r")
sptmp:=short2flt(f[:]->sp(::3,:,:)) ; Surface pressure

tp:=tptmp
tp(:,:,0:719)=tptmp(:,:,720:1439)
tp(:,:,720:1439)=tptmp(:,:,0:719)
tcc:=tcctmp
tcc(:,:,0:719)=tcctmp(:,:,720:1439)
tcc(:,:,720:1439)=tcctmp(:,:,0:719)
sp:=sptmp
sp(:,:,0:719)=sptmp(:,:,720:1439)
sp(:,:,720:1439)=sptmp(:,:,0:719)
lsp:=lsptmp
lsp(:,:,0:719)=lsptmp(:,:,720:1439)
lsp(:,:,720:1439)=lsptmp(:,:,0:719)
cp:=cptmp
cp(:,:,0:719)=cptmp(:,:,720:1439)
cp(:,:,720:1439)=cptmp(:,:,0:719)

dims=dimsizes(tp)
nt=dims(0)
print(dims)

if(isleapyear(yy0)) then
tp_daily:=new((/366,ny,nx/),float)
tp_daily!0="jday"
tp_daily!1="lat"
tp_daily!2="lon"
tp_daily&jday=ispan(1,366,1)
tp_daily&lat=lat1
tp_daily&lon=lon1
tcc_daily:=tp_daily
sp_daily:=tp_daily
lsp_daily:=tp_daily
cp_daily:=tp_daily
else
tp_daily:=new((/365,ny,nx/),float)
tp_daily!0="jday"
tp_daily!1="lat"
tp_daily!2="lon"
tp_daily&jday=ispan(1,365,1)
tp_daily&lat=lat1
tp_daily&lon=lon1
tcc_daily:=tp_daily
sp_daily:=tp_daily
lsp_daily:=tp_daily
cp_daily:=tp_daily
end if

do it=0,nt-1,8
id=toint(it/8.0)
;print(id)
tp_daily(id,:,:)=area_hi2lores(lon0,lat0,dim_avg_n(tp(it:it+7,:,:),0),True,1,lon1,lat1,0)
sp_daily(id,:,:)=area_hi2lores(lon0,lat0,dim_avg_n(sp(it:it+7,:,:),0),True,1,lon1,lat1,0)
tcc_daily(id,:,:)=area_hi2lores(lon0,lat0,dim_avg_n(tcc(it:it+7,:,:),0),True,1,lon1,lat1,0)  
lsp_daily(id,:,:)=area_hi2lores(lon0,lat0,dim_avg_n(lsp(it:it+7,:,:),0),True,1,lon1,lat1,0)  
cp_daily(id,:,:)=area_hi2lores(lon0,lat0,dim_avg_n(cp(it:it+7,:,:),0),True,1,lon1,lat1,0)  
end do

setfileoption("nc", "Format",  "NetCDF4")
	system("rm -f erapost/erapost-precip"+yy0+".nc")
fout=addfile("erapost/erapost-precip"+yy0+".nc","c")
fout->sp_daily = sp_daily
fout->tp_daily = tp_daily
fout->tcc_daily = tcc_daily
fout->lsp_daily = lsp_daily
fout->cp_daily = cp_daily
fout->lat = lat1
fout->lon = lon1
end
