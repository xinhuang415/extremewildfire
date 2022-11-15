clear
%EMD + FFT <- FINN data
%% --------Region define---------------
lat_S(1)=-16;%South Africa 
lat_N(1)=-4;
lon_W(1)=15;
lon_E(1)=35;
location(1).name='South_Africa';
location_real(1).name='South_Africa';

lat_S(2)=4;%North Africa 
lat_N(2)=11;
lon_W(2)=-10;
lon_E(2)=36;
location(2).name='North_Africa';
location_real(2).name='North_Africa';

lat_S(3)=-20;%Amazon 
lat_N(3)=-5;
lon_W(3)=-70;
lon_E(3)=-40;
location(3).name='Amazon';
location_real(3).name='Amazon';

lat_S(4)=15;%Southeast Asia
lat_N(4)=25;
lon_W(4)=95;
lon_E(4)=105;
location(4).name='Southeast_Asia';
location_real(4).name='Southeast_Asia';

lat_S(5)=35;%USA west coast
lat_N(5)=50;
lon_W(5)=-125;
lon_E(5)=-115;
location(5).name='USA_west_coast';
location_real(5).name='USA_west_coast';

lat_S(6)=12;%South Mexico
lat_N(6)=20;
lon_W(6)=-105;
lon_E(6)=-82;
location(6).name='South_Mexico';
location_real(6).name='South_Mexico';

lat_S(7)=36;%Europe West
lat_N(7)=45;
lon_W(7)=-10;
lon_E(7)=0;
location(7).name='Europe_West';
location_real(7).name='Europe_West';

region_num=size(lat_S,2);
%% ------------------------FINN data input------------------------
year_start=2002;
year_end=2021;
year_num=year_end-year_start+1;

Area_FINN_all(1:720,1:360,year_num)=0;
Area_regions(1:region_num,1:year_num,1:365)=0;

folder_FINN='G:\FINN_burned_area\';
for year=year_start:year_end%2021
    file_in=[folder_FINN,'fireco',num2str(year),'.nc'];
    ncdisp(file_in);
    
    lon_FINN=ncread(file_in,'lon');%720
    lat_FINN=ncread(file_in,'lat');%360
    Area_FINN=ncread(file_in,'area');%720 360 366 lon lat day
    index=find(Area_FINN>1e36);
    Area_FINN(index)=nan;
    
    Area_FINN_all(:,:,year-(year_start-1))=squeeze(mean(Area_FINN,3));
    
    for region_i=1:region_num
        
        lat_S_num=round((lat_S(region_i)-lat_FINN(1))/0.5)+1;
        lat_N_num=round((lat_N(region_i)-lat_FINN(1))/0.5)+1;
        lon_W_num=round((lon_W(region_i)-lon_FINN(1))/0.5)+1;
        lon_E_num=round((lon_E(region_i)-lon_FINN(1))/0.5)+1;
        
        Area_regions(region_i,year-(year_start-1),1:365)=squeeze(mean(Area_FINN(lon_W_num:lon_E_num,lat_S_num:lat_N_num,1:365),[1,2]));
                 
    end
    
end

Area_regions_yearly=squeeze(mean(Area_regions,3));
Area_regions_std=std(Area_regions_yearly,0,2);
Area_regions_mean=mean(Area_regions_yearly,2);
Area_interal_Var=Area_regions_std./Area_regions_mean*100;%


Area_FINN_mean=squeeze(mean(Area_FINN_all,3));

colorname='WhiteBlueGreenYellowRed';
color_path=['D:\study\Matlab_program_2018_09_09\colorbar_NCL\',colorname,'.txt'];
color_path_bias=['D:\study\Matlab_program_2018_09_09\colorbar_NCL\temp_diff_18lev.txt'];

load coast
figure(1)
pcolor(lon_FINN,lat_FINN,Area_FINN_mean');
shading flat
colorbar_want=load(color_path);
colormap(colorbar_want);
colorbar;
axis([-150 170 -40 50]);
caxis([0 2e6]);
hold on
plot(long,lat,'k');

for region_i=1:region_num
    x_line=[lon_W(region_i),lon_E(region_i),lon_E(region_i),lon_W(region_i),lon_W(region_i)];
    y_line=[lat_S(region_i),lat_S(region_i),lat_N(region_i),lat_N(region_i),lat_S(region_i)];
    plot(x_line,y_line,'k');  
end

%% -------------------EMD and FFT---------------------
fs=1;
ll=0;
station_ll_num(1:region_num)=0;
for station_i=1:region_num
    for year_i=1:year_num
        x_tmp=squeeze(Area_regions(station_i,year_i,:));
        x_tmp_nor=squeeze(normalize(x_tmp));
        imf = emd(x_tmp_nor);%EMD
        %plot_hht(data_tmp,1/Fs)
        imf_num=size(imf,2);
        
        station_ll_num(station_i)=station_ll_num(station_i)+imf_num;
        for i=1:imf_num
            imf_tmp=imf{i};
     
            n=length(imf_tmp);
            y_tmp=fft(imf_tmp);%FFT
            f_tmp=(0:n-1)*(fs/n);
            Time_tmp=1./f_tmp(2:n);
            power=abs(y_tmp).^2/n;%power of the DFT    
            
            ll=ll+1;
            data(ll).T=Time_tmp;
            data(ll).Power=power(2:n);    
        end

    end
 
end

ll_num=size(data,2);
for i=1:ll_num 
    T_all(i,1:364)=data(i).T; 
    Power_all(i,1:364)=data(i).Power;
end


%% ------------------NC OUT----------------------------------------
    out_folder='D:\study\NEW_paper\Rain_fire_relation\Data\';
    out_name=[ 'EMD_FFT_est_Period_Power_for_7region.nc'];
    ncid=netcdf.create([out_folder,out_name],'NC_64BIT_OFFSET');
    file_out=[out_folder,out_name];
    %----------  dimension define-----------
    for station_i=1:region_num
        varname_tmp=['dim_imf_',location_real(station_i).name];
        eval([varname_tmp,'=netcdf.defDim(ncid,''',varname_tmp,''',station_ll_num(station_i));']);
    end
    
    dim_Period=netcdf.defDim(ncid,'Period',364);
    %----------   var define    --------------
    for station_i=1:region_num
        varname_out_tmp=['Var_out_',location_real(station_i).name];
        varname_tmp=['dim_imf_',location_real(station_i).name];
        powername_tmp=['Power_',location_real(station_i).name];
        eval([varname_out_tmp,'=netcdf.defVar(ncid,''',powername_tmp,''',''NC_FLOAT'',[',varname_tmp,',dim_Period]);']);
    end
    Period_out=netcdf.defVar(ncid,'Period','NC_FLOAT',[dim_Period]);
    %-------   netcdf end ------------
    netcdf.endDef(ncid)
   %------------ put var -------------
   for station_i=1:region_num
       ll_start=sum(station_ll_num(1:station_i-1))+1;
       ll_end=sum(station_ll_num(1:station_i));
       varname_out_tmp=['Var_out_',location_real(station_i).name];
       eval(['netcdf.putVar(ncid,',varname_out_tmp,',single(Power_all(ll_start:ll_end,:)));']);
   end

   netcdf.putVar(ncid,Period_out,single(squeeze(T_all(1,:))));
   % -------close--------
   netcdf.close(ncid);
   %--------add attribure--------------
   for station_i=1:region_num
       lat_S_str=num2str(lat_S(station_i));
       lat_N_str=num2str(lat_N(station_i));
       lon_W_str=num2str(lon_W(station_i));
       lon_E_str=num2str(lon_E(station_i));
       Attr=['lat= ',lat_S_str,' to ',lat_N_str,' lon=',lon_W_str,' to ',lon_E_str];
       
       powername_tmp=['Power_',location_real(station_i).name];
       eval(['ncwriteatt([out_folder,out_name],''',powername_tmp,''',''Attribute'',Attr);']);
       
   end
   
   


