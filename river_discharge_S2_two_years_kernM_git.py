# -*- coding: utf-8 -*-
"""
@author: p.filippucci
@description: Python code to extract a proxy of river discharge from reflectance 
            data of S2 using GEE platform
@date: 29/08/2024
@version: 1
    
"""


import ee
import numpy as np
from datetime import datetime,timedelta
import time
import pandas
from scipy.io import savemat
from scipy.io import loadmat as loadmat2

                
# %% Input

# USER GEE login https://developers.google.com/earth-engine/guides/service_account
service_account = 'addyourGEEserviceaccount'
credentials = ee.ServiceAccountCredentials(service_account, 'addpathtoyourGEEcredentials')

ee.Initialize(credentials)
folderGEE='addyourprojectpath'; #e.g.: projects/riverdischarge/assets/research

# Select desired features
extr_video=0                                             #put 1 to extract gif of calibration period
pix_W1=1                                                 #put 1 to use std*avg W pixel
pix_W2=1                                                 #put 1 to use avg W pixel
veg=1                                                    #put 1 to extract calibrated M pixels using CV
no_veg=1                                                 #put 1 to extract calibrated M pixels using C
pix_Mcal_k1_CM=1                                         #put 1 to extract the calibrated M pixel kernel xk1 with CM
pix_Mcal_k1_CMW=1                                        #put 1 to extract the calibrated M pixel kernel xk1 with CMW
pix_Mcal_k2_CM=1                                         #put 1 to extract the calibrated M pixel kernel xk2 with CM
pix_Mcal_k2_CMW=1                                        #put 1 to extract the calibrated M pixel kernel xk2 with CMW    

# Parameters (can be changed)
ray=0.020                                                #half the side of the square centered in the coordinates
snowthre=10                                              #threshold on max snow probability
folder_pc='addthepathtoyoursavingfolder'                 #folder where to save the data

# USER insert below station name, central area coordinates and river coordinates
# below and example with data and coordinated for Pontelagoscuro
file='Po_5sta_HQvAB.xlsx'                                 #Exemplary file with observed data 
cols='B:J'
rws=6
sh_nm='LAG' 
Qobs=pandas.read_excel(file, sheet_name=sh_nm, usecols=cols,skiprows=rws)
Dobs0=np.array(Qobs['Date'])                               #observed dates
Qobs=np.array(Qobs['DISCHARGE'])*1.                        #observed discharge
coordstat=[11.606040000915527, 44.88702392578125]          #longitude-latitude of the station
name='Pontelagoscuro'                                      #station name
str2='Po__uro'                                             #station code
river='Po'                                                 #river name
kern1=7                                                    #1st dimension of the kernel
kern2=8                                                    #2nd dimension of the kernel
nar=0                                                      #1 if the river is too narrow for the reduction of the water mask, 0 otherwise
years_cal=[2016,2019]                                      #Calibration period (has to be < 5000 days)
years=[2000,2005,2010,2015,2020,2022]                      #Extraction periods (divided to reduce the computating weight)
self_Wat=0                                                 #equal to 1 if water detecting algorithm do not work. insert river central line in coordW
hydro0='Discharge'                                         #name of the observed variable
hydro='Dis'                                                #code of the observed variable
coordW=[[11.596556646173383, 44.889278173147744],
        [11.625781995599652, 44.888366036430725],
        [11.612436, 44.888439]]                            #water pixel coordinates


# %% Initialization

Qobs[Qobs<0]=np.nan
DD=np.array([datetime.fromordinal(int(d-366)) for d in Dobs0])
prod='S2_HARMONIZED'
prod_name='S2_1C'
ray0=str(int(ray*1000)).zfill(4)

attempt=1
fail=0
while attempt==1:
    try:
        fail=fail+1
        if fail>5:
            attempt=0
        ee.Reset()
        ee.Initialize(credentials)
        print(str2)
        area1=[coordstat[0]-ray,coordstat[1]-ray,coordstat[0]+ray,coordstat[1]+ray]
        mult=0.0001
        cld_thr=50 
        region = ee.Geometry.Rectangle(area1);
        
        #Set observed data in Q if calibration is active (otherwise, -9999)
        try:
            Q=ee.FeatureCollection(folderGEE+'/data_'+str2)
            Q.first().getInfo()
        except:
            Dobs=np.array([datetime.fromordinal(int(d-366)) for d in Dobs0])
            ID=[d.year>=1999 for d in Dobs]
            Dobs=Dobs[ID]
            Qobs[np.isnan(Qobs)]=-9999
            
            polygon = ee.Geometry.Point([coordstat[0], coordstat[1]]);
            feat=[]
            for i in range (Qobs.size):
                feat.append(ee.Feature(polygon).set('Date',(Dobs[i]-datetime(1970,1,1)).days*3600*24*1000).set(hydro0,Qobs[i]))
            
            Q=ee.FeatureCollection(feat)
            task = ee.batch.Export.table.toAsset(collection=Q,
                                             description='Observed_data_'+str2,
                                             assetId=folderGEE+'/data_'+str2)
            task.start()
        
            while task.status()['state']!='COMPLETED' and task.status()['state']!='FAILED':
                print('saving')
                time.sleep(10)
            print(task.status()['state'])
            
            Q=ee.FeatureCollection(folderGEE+'/data_'+str2)

        # %% snow period Calculation 
        
        try:
            file=loadmat2('snowmask_'+str2+'.mat')
            IDclim=file['IDclim']
        except:
            print('snow period calculation')
            area2=[coordstat[0]-0.02,coordstat[1]-0.02,coordstat[0]+0.02,coordstat[1]+0.02]
            region2=ee.Geometry.Rectangle(area2)
            snowMOD=ee.ImageCollection("MODIS/061/MOD10A1").select('NDSI_Snow_Cover').filterDate(ee.Date('2003-01-01'),ee.Date('2020-01-01')).filterBounds(region2)
            snowMYD=ee.ImageCollection("MODIS/061/MYD10A1").select('NDSI_Snow_Cover').filterDate(ee.Date('2003-01-01'),ee.Date('2020-01-01')).filterBounds(region2)
            
            proj=snowMYD.first().projection();
            scale_mod=ee.Number(proj.nominalScale()).getInfo()
            trans_mod=proj.getInfo()['transform'];
            crs_mod=proj.getInfo()['crs'];
            WW=ee.Image('JRC/GSW1_2/GlobalSurfaceWater').select('max_extent')
            
            def masking(image):
                return image.updateMask(WW.Not())
            
            snowMOD=snowMOD.map(masking)
            snowMYD=snowMYD.map(masking)
        
            snowMOD = np.array(snowMOD.getRegion(region2, scale_mod).getInfo())
            snowMYD = np.array(snowMYD.getRegion(region2, scale_mod).getInfo())
            snowMOD[snowMOD==None]=np.nan
            snowMYD[snowMYD==None]=np.nan
            
            Dmod_full=np.array([datetime.strptime(d,'%Y_%m_%d') for d in snowMOD[1:,0]])
            Dmod=np.unique(Dmod_full)
            Dmyd_full=np.array([datetime.strptime(d,'%Y_%m_%d') for d in snowMYD[1:,0]])
            Dmyd=np.unique(Dmyd_full)
            
            snowMOD_2=np.zeros(Dmod.shape)
            snowMYD_2=np.zeros(Dmyd.shape)
            
            for i,day in enumerate(Dmod):
                ID=np.where(Dmod_full==day)[0]+1   
                snowMOD_2[i]=np.nanmean(snowMOD[ID,4])
            for i,day in enumerate(Dmyd):
                ID=np.where(Dmyd_full==day)[0]+1   
                snowMYD_2[i]=np.nanmean(snowMYD[ID,4])
                
            snowMOD_clim=np.zeros([365])
            snowMYD_clim=np.zeros([365])
            snowMOD_clim[:]=np.nan
            snowMYD_clim[:]=np.nan
            
            Dmod=np.array([int(d.strftime('%j')) for d in Dmod])
            Dmyd=np.array([int(d.strftime('%j')) for d in Dmyd])
            
            for i in range(1,366):
                snowMOD_clim[i-1]=np.nanmean(snowMOD_2[Dmod==i])
                snowMYD_clim[i-1]=np.nanmean(snowMYD_2[Dmyd==i])
            
            snowMOD_clim2=np.concatenate((snowMOD_clim,snowMOD_clim[0:7]))
            snowMYD_clim2=np.concatenate((snowMYD_clim,snowMYD_clim[0:7]))
            snowMOD_clim2[np.isnan(snowMOD_clim2)]=100
            snowMYD_clim2[np.isnan(snowMYD_clim2)]=100

            snowMOD_clim_filt=np.zeros([365])
            snowMYD_clim_filt=np.zeros([365])
            snowMOD_clim_filt[:]=np.nan
            snowMYD_clim_filt[:]=np.nan
            
            for i in range(0,365):
                if i-7<0:
                    snowMOD_clim_filt[i]=np.nanmean([np.nanmean(snowMOD_clim2[i-7:-1]),np.nanmean(snowMOD_clim2[0:i+7])])
                    snowMYD_clim_filt[i]=np.nanmean([np.nanmean(snowMYD_clim2[i-7:-1]),np.nanmean(snowMYD_clim2[0:i+7])])
                else:
                    snowMOD_clim_filt[i]=np.nanmean(snowMOD_clim2[i-7:i+7])
                    snowMYD_clim_filt[i]=np.nanmean(snowMYD_clim2[i-7:i+7])
                
            IDclim=np.logical_and(snowMYD_clim_filt[:]>snowthre,snowMOD_clim_filt[:]>snowthre)
            IDclim2=np.zeros(IDclim.size)
            IDclim=np.append(IDclim,IDclim);
            n2=-1
            n1=0
            for j in range(IDclim.size):
                if IDclim[j]:
                    n1=n1+1
                    app=365
                elif n1!=0:   
                    if n2+1>=365:
                        n2=n2-365
                    elif n2+1+n1>=365:
                        IDclim2[n2+1:365]=n1
                        app=n2
                        n2=-1
                    IDclim2[n2+1:n2+1+n1-(364-app)]=n1
                    n2=j
                    n1=0
                else:
                    n2=j
                    
            IDclim2=IDclim2>30
            IDclim=np.where(IDclim2)[0]
            if 364 in IDclim:
                IDclim=np.append(IDclim,365)            
            savemat('snowmask_'+str2+'.mat',{'IDclim':IDclim})
    
        # %% Collection creation and Water Mask Calculation full period
    
        datestartCAL=ee.Date(str(years_cal[0])+'-01-1').millis();
        dateendCAL=ee.Date(str(years_cal[-1])+'-12-31').millis();
        d1=ee.Date(datestartCAL).format('YYYYMMDD').getInfo()
        d2=ee.Date(dateendCAL).format('YYYYMMDD').getInfo()

        print('product: '+prod)
        band=['B8','B4'];

        product1='COPERNICUS/'+prod;
        product2='COPERNICUS/S2_CLOUD_PROBABILITY';

        coll=ee.ImageCollection(product1).filterDate(datestartCAL,dateendCAL).filterBounds(region).select(band);
        proj=coll.first().projection();
        scale_mod=ee.Number(proj.nominalScale()).getInfo()
        trans_mod=proj.getInfo()['transform'];
        crs_mod=proj.getInfo()['crs'];
        im00=ee.Image().reproject(crs_mod,trans_mod).select(['constant'],['Q']).double();
        im01=im00.rename('C');
        im02=im00.rename('W');
        im03=im00.rename('M');
        im04=im00.rename('V');


        coll2=ee.ImageCollection(product2).filterDate(datestartCAL,dateendCAL).filterBounds(region);
        coll=coll.combine(coll2)
        print('original length: '+str(coll.size().getInfo()))

        
        # %% Cloud and snow days filtering
        
        first=ee.List([])
        def calc_Date(image,llist):
            value =image.get('system:time_start')
            value=ee.Number(ee.Date(ee.List([value, -9999]).reduce(ee.Reducer.firstNonNull())).getRange('day').start().millis());
            return ee.List(llist).add(value)
        dlist=ee.List(coll.iterate(calc_Date,first)).distinct().getInfo()  
        dlist=np.sort(dlist).tolist()
        
        band2=band+['probability']    
        
        def selection(nday):
            image=coll.filterDate(ee.Date(nday),ee.Date(ee.Number(nday).add(86400000))).reduce(ee.Reducer.mean()).rename(band2).reproject(crs_mod,trans_mod)
            image=image.select(band).updateMask(image.select('probability').lte(cld_thr)).multiply(mult).set('system:time_start',nday)
            image=image.updateMask(image.gt(0))
            
            value=image.select(band[0]).mask().reduceRegion(**{'reducer':ee.Reducer.mean(),'geometry':region,'crs':crs_mod,'crsTransform':trans_mod,'bestEffort':True})
            value=ee.Number(value.get(value.keys().get(0)))
            return image.set('Valid_perc',value)
        
        collection0=ee.ImageCollection(ee.List(dlist).map(selection))        
                       
        collection0=collection0.filterMetadata('Valid_perc',"not_less_than",0.20)
        print(collection0.size().getInfo())
        
        dlistn=ee.List(collection0.iterate(calc_Date,first)).getInfo()
        

        try:
            unvalid=ee.Image(folderGEE+'/Unvalid_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr)).reproject(crs_mod,trans_mod)
            unvalid.getInfo()
        except:
            coll2bis=coll2.map(lambda image: image.updateMask(image.lte(50)))
            unvalid=(ee.Image(coll2bis.count()).divide(ee.Number(coll2bis.size())).gt(0.05))
            
            task = ee.batch.Export.image.toAsset(image=unvalid,
                                             description='Unvalid_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr),
                                             region=region,
                                             assetId= folderGEE+'/Unvalid_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr),
                                             crs=crs_mod,
                                             crsTransform=trans_mod,
                                             maxPixels=1000000000)
            task.start()
            
            while task.status()['state']!='COMPLETED' and task.status()['state']!='FAILED':
                print('saving')
                time.sleep(10)
            print(task.status()['state'])
        
            unvalid=ee.Image(folderGEE+'/Unvalid_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr)).reproject(crs_mod,trans_mod)
        
                    


        print('after cloud masking length: '+str(len(dlistn)))
        dlist0=np.array(dlistn)/1000/3600/24
        D=np.array ([datetime(1970,1,1)+timedelta(days=d) for d in dlist0])
        dlist0=np.array ([(datetime(1970,1,1)+timedelta(days=d)).timetuple().tm_yday-1 for d in dlist0])
        ID=[np.where([~np.any(d==IDclim) for d in dlist0])[0]][0]
        dlistn=np.array(dlistn)[ID].tolist()
        print('after snow masking length: '+str(len(dlistn)))

        def mosaic(nday):
            image=coll.filterDate(ee.Date(nday),ee.Date(ee.Number(nday).add(86400000))).reduce(ee.Reducer.mean()).rename(band2).reproject(crs_mod,trans_mod)
            image=image.select(band).updateMask(image.select('probability').lte(cld_thr).And(unvalid)).multiply(mult).set('system:time_start',nday)
            image=image.updateMask(image.gt(0))
            
            feat = Q.filter(ee.Filter.eq('Date', nday)).first();
            obs=im00.unmask(ee.Number(feat.get(hydro0)));
            obs=obs.updateMask(obs.neq(-9999));
        
            return image.addBands(obs).set('system:time_start',nday)
        
        collection=ee.ImageCollection(ee.List(dlistn).map(mosaic))
        
        def mosaicNDVI(nday):
            image=coll.filterDate(ee.Date(nday),ee.Date(ee.Number(nday).add(86400000))).reduce(ee.Reducer.mean()).rename(band2).reproject(crs_mod,trans_mod)
            image=image.select(band).updateMask(image.select('probability').lte(cld_thr)).multiply(mult).set('system:time_start',nday)
            image=image.updateMask(image.gt(0))
            
            return image.select(band[0]).subtract(image.select(band[1])).divide(image.select(band[0]).add(image.select(band[1]))).rename('NDVI').set('system:time_start',nday)
        
        collectionN=ee.ImageCollection(ee.List(dlistn).map(mosaicNDVI))
        
        band=[band[0]]
        
        if extr_video==1:
            minM2=collection.min().reduceRegion(**{'reducer': ee.Reducer.min(),'geometry': region,'crs': crs_mod,'crsTransform':trans_mod,'bestEffort': True}).get(band[0]).getInfo();
            maxM2=collection.max().reduceRegion(**{'reducer': ee.Reducer.max(),'geometry': region,'crs': crs_mod,'crsTransform':trans_mod,'bestEffort': True}).get(band[0]).getInfo();
            val_max = maxM2;
            val_min = minM2;
            coll4Video = collection.map(lambda image: image.select(band).visualize(**{'min':minM2,'max':maxM2,'opacity': 1.0,'palette': ["black", "white"]}));
            task = ee.batch.Export.video.toDrive(**{
                'collection': coll4Video,
                'description': ('video_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr)),
                'crs': crs_mod,
                'framesPerSecond': 10,
                'crsTransform':trans_mod,
                'region': region
            });
            task.start()
            
            while task.status()['state']!='COMPLETED' and task.status()['state']!='FAILED':
                print('saving')
                time.sleep(10)
            print(task.status()['state'])


        # %% Water Mask calculation
        print('start water masking')
        cv_reducersC = ee.Reducer.mean().combine(**{
            'reducer2':ee.Reducer.stdDev(),'sharedInputs':True});
    
        try:
            Wat=ee.Image(folderGEE+'/Wat_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr)).reproject(crs_mod,trans_mod)
            Wat.getInfo()

            Wat2=ee.Image(folderGEE+'/Wat2_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr)).reproject(crs_mod,trans_mod)
            Wat2.getInfo()

        except:
            if self_Wat==0:
                Wat=collectionN.reduce(ee.Reducer.percentile([5])).reproject(crs_mod,trans_mod).lte(-0.05)
                Wat2=Wat.reduceNeighborhood(ee.Reducer.mean(), ee.Kernel.square(1-nar,'pixels')).reproject(crs_mod,trans_mod).eq(1).reproject(crs_mod,trans_mod)
                Wat=Wat.reduceNeighborhood(ee.Reducer.mean(), ee.Kernel.square(3,'pixels')).reproject(crs_mod,trans_mod).gt(0).reproject(crs_mod,trans_mod);

                vectors = Wat2.reduceToVectors(**{
                  'geometry': region,
                  'crs': crs_mod,
                  'crsTransform': trans_mod,
                  'geometryType': 'polygon',
                  'eightConnected': False,
                  'labelProperty': 'zone'});
                def river_select(feat):
                    ins=feat.intersects(ee.Geometry.MultiPoint(coordW),1,proj) 
                    return feat.set('inside',ins)

                solv=vectors.map(river_select)
                solv=solv.filter(ee.Filter.eq('inside',True))
                Wat2=Wat2.clip(solv).unmask(0).reduceNeighborhood(ee.Reducer.mean(), ee.Kernel.square(4-nar,'pixels')).reproject(crs_mod,trans_mod).gt(0).reproject(crs_mod,trans_mod);     
                
                task = ee.batch.Export.image.toAsset(image=Wat,
                                                 description='Wat_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr),
                                                 region=region,
                                                 assetId= folderGEE+'/Wat_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr),
                                                 crs=crs_mod,
                                                 crsTransform=trans_mod,
                                                 maxPixels=1000000000)
                task2 = ee.batch.Export.image.toAsset(image=Wat2,
                                                 description='Wat2_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr),
                                                 region=region,
                                                 assetId= folderGEE+'/Wat2_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr),
                                                 crs=crs_mod,
                                                 crsTransform=trans_mod,
                                                 maxPixels=1000000000)
                task.start()
                task2.start()
                
                while (task.status()['state']!='COMPLETED' and task.status()['state']!='FAILED') or (task2.status()['state']!='COMPLETED' and task2.status()['state']!='FAILED'):
                    print('saving')
                    time.sleep(10)
                print(task.status()['state'])
            
                Wat=ee.Image(folderGEE+'/Wat_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr)).reproject(crs_mod,trans_mod)
                Wat2=ee.Image(folderGEE+'/Wat2_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr)).reproject(crs_mod,trans_mod)
                
            else:
                l=ee.Geometry.LineString(coordW)
                l=l.buffer(4*scale_mod)
                im=ee.Image(1).reproject(crs_mod,trans_mod).clip(l)
                im=ee.Image(1).updateMask(im).unmask(0).reproject(crs_mod,trans_mod).gt(0).reproject(crs_mod,trans_mod)

                task = ee.batch.Export.image.toAsset(image=im,
                                                 description='Wat_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr),
                                                 region=region,
                                                 assetId= folderGEE+'/Wat_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr),
                                                 crs=crs_mod,
                                                 crsTransform=trans_mod,
                                                 maxPixels=1000000000)
                task2 = ee.batch.Export.image.toAsset(image=im,
                                                 description='Wat2_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr),
                                                 region=region,
                                                 assetId= folderGEE+'/Wat2_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr),
                                                 crs=crs_mod,
                                                 crsTransform=trans_mod,
                                                 maxPixels=1000000000)
                task.start()
                task2.start()
                
                while (task.status()['state']!='COMPLETED' and task.status()['state']!='FAILED') or (task2.status()['state']!='COMPLETED' and task2.status()['state']!='FAILED'):
                    print('saving')
                    time.sleep(10)
                print(task.status()['state'])
            
                Wat=ee.Image(folderGEE+'/Wat_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr)).reproject(crs_mod,trans_mod)
                Wat2=ee.Image(folderGEE+'/Wat2_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr)).reproject(crs_mod,trans_mod)


        # %% Mask Calculation 
                    
        tasks=[]

        # CALCULATE C2
        print('C extraction')
        try:
            C2mask=ee.Image(folderGEE+'/Cmask_2_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr)).reproject(crs_mod,trans_mod)                    
            C2mask.getInfo()
                           
        except:
            
            im_C = collection.select(band).reduce(cv_reducersC).reproject(crs_mod,trans_mod);
            im_C= im_C.expression('b(1) / b(0)').rename('C').updateMask(Wat.Not());
            mm = im_C.reduceRegion(**{'reducer':ee.Reducer.percentile([2]), 'geometry':region, 'crs':crs_mod, 'crsTransform':trans_mod, 'bestEffort':True})
            C2mask = im_C.lt(ee.Number(mm.get(mm.keys().get(0))))
        
            task1 = ee.batch.Export.image.toAsset(image=C2mask,
                                             description='Cmask_2_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr),
                                             region=region,
                                             assetId= folderGEE+'/Cmask_2_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr),
                                             crs=crs_mod,
                                             crsTransform=trans_mod,
                                             maxPixels=1000000000)
            task1.start()
            tasks.append(task1)
            
        # CALCULATE V
        if veg==1:
            print('V extraction')
            try:
                Vmask=ee.Image(folderGEE+'/Vmask_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr)).reproject(crs_mod,trans_mod)
                Vmask.getInfo()
            except:
                M_NDVI=collectionN.select('NDVI').mean().rename('V').updateMask(Wat.Not());
                mm5 = M_NDVI.reduceRegion(**{'reducer':ee.Reducer.percentile([95]), 'geometry':region, 'crs':crs_mod, 'crsTransform':trans_mod ,'bestEffort':True})
                Vmask= M_NDVI.gt(ee.Number(mm5.get(mm5.keys().get(0))))
                
                task2 = ee.batch.Export.image.toAsset(image=Vmask,
                                             description='Vmask_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr),
                                             region=region,
                                             assetId= folderGEE+'/Vmask_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr),
                                             crs=crs_mod,
                                             crsTransform=trans_mod,
                                             maxPixels=1000000000)
                task2.start()
                tasks.append(task2)

        # CALCULATE W
        if pix_W1==1:        
            print('W1 extraction')
            try:
                Wmask=ee.Image(folderGEE+'/Wmask_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr)).reproject(crs_mod,trans_mod)
                Wmask.getInfo()
            except:
                
                im_W = collection.select(band).reduce(cv_reducersC).reproject(crs_mod,trans_mod);
                im_W= im_W.expression('b(1) * b(0)').updateMask(Wat2).rename('W');
                mm = im_W.reduceRegion(**{'reducer':ee.Reducer.percentile([5]), 'geometry':region, 'crs':crs_mod,'crsTransform':trans_mod,  'bestEffort':True})
                Wmask = im_W.lt(ee.Number(mm.get(mm.keys().get(0))))
                
                task3 = ee.batch.Export.image.toAsset(image=Wmask,
                                                 description='Wmask_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr),
                                                 region=region,
                                                 assetId= folderGEE+'/Wmask_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr),
                                                 crs=crs_mod,
                                                 crsTransform=trans_mod,
                                                 maxPixels=1000000000)
                task3.start()
                tasks.append(task3)

        if pix_W2==1:                        
            print('W2 extraction')
            try:
                Wmask2=ee.Image(folderGEE+'/Wmask2_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr)).reproject(crs_mod,trans_mod)
                Wmask2.getInfo()
            except:
                
                im_W = collection.select(band).reduce(cv_reducersC).reproject(crs_mod,trans_mod);
                im_W= im_W.expression('b(0)').updateMask(Wat2).rename('W');
                        
                mm = im_W.reduceRegion(**{'reducer':ee.Reducer.percentile([5]), 'geometry':region, 'crs':crs_mod,'crsTransform':trans_mod,  'bestEffort':True})
                Wmask2 = im_W.lt(ee.Number(mm.get(mm.keys().get(0))))
                
                task4 = ee.batch.Export.image.toAsset(image=Wmask2,
                                                 description='Wmask2_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr),
                                                 region=region,
                                                 assetId= folderGEE+'/Wmask2_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr),
                                                 crs=crs_mod,
                                                 crsTransform=trans_mod,
                                                 maxPixels=1000000000)
                task4.start()
                tasks.append(task4)
                
                
        if len(tasks)>0:
            while any([t.status()['state']!='COMPLETED' and t.status()['state']!='FAILED' for t in tasks]):
                print('saving')
                time.sleep(10)
            print([t.status()['state'] for t in tasks])

            try:
                C2mask=ee.Image(folderGEE+'/Cmask_2_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr)).reproject(crs_mod,trans_mod)                    
                C2mask.getInfo()
            except:
                C2mask=im01
                            
            if veg==1:
                try:
                    Vmask=ee.Image(folderGEE+'/Vmask_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr)).reproject(crs_mod,trans_mod)
                    Vmask.getInfo()
                except:
                    Vmask=im04

            if pix_W1==1:
                try:
                    Wmask=ee.Image(folderGEE+'/Wmask_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr)).reproject(crs_mod,trans_mod)
                    Wmask.getInfo()
                except:
                    Wmask=im02

            if pix_W2==1:
                try:
                    Wmask2=ee.Image(folderGEE+'/Wmask2_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr)).reproject(crs_mod,trans_mod)
                    Wmask2.getInfo()
                except:
                    Wmask2=im02

                    
        # %% M-cal calculation

        if no_veg==1:                
            if pix_Mcal_k1_CM==1: 
                # calculate calibrated CM M k1 pixel
                print('M CM k1 cal extraction')                     
                try:
                    Mmask_C2M_k1=ee.Image(folderGEE+'/Mmask_C2_calib_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr)+'_var_'+hydro+'_k'+str(kern1)).reproject(crs_mod,trans_mod)
                    Mmask_C2M_k1.getInfo()
               
                except: 
                    ee.Reset()
                    ee.Initialize(credentials)
                    def calc_points(image):
                        value =image.select(band).updateMask(C2mask).reduceRegion(**{'reducer':ee.Reducer.mean(),
                                                                  'geometry': region,
                                                                  'crs': crs_mod,
                                                                  'crsTransform':trans_mod,
                                                                  'bestEffort': True});
                        value=ee.Number(ee.List([value.get(value.keys().get(0)), -9999]).reduce(ee.Reducer.firstNonNull()));
                        im = im01.unmask(value).updateMask(value.neq(-9999));
                        return image.addBands(im);
                    collection0 = collection.map(calc_points);
                
                    coll02=collection0.map(lambda image: image.select(['Q','C']).addBands(image.select(band).reduceNeighborhood(ee.Reducer.mean(), ee.Kernel.square(kern1,'pixels',True)).rename(band).reproject(crs_mod,trans_mod)));            
                    CM=coll02.map(lambda image: image.select('Q').addBands(image.select('C').divide(image.select(band))));
            
                    # calculate correlation between discharge and C/M
                    corr=CM.reduce(ee.Reducer.spearmansCorrelation()).reproject(crs_mod,trans_mod).updateMask(Wat2);
                    
                    # find maximum correlation Calibration and create mask for that pixel
                    maxval = corr.reduceRegion(**{'reducer':ee.Reducer.max(), 'geometry':region, 'crs':crs_mod,'crsTransform':trans_mod, 'bestEffort':True});
                    Mmask_C2M_k1 = corr.select('correlation').eq(ee.Number(maxval.get('correlation')));
                    Mmask_C2M_k1 = Mmask_C2M_k1.reduceNeighborhood(ee.Reducer.mean(), ee.Kernel.square(kern1,'pixels',True)).rename('correlation').reproject(crs_mod,trans_mod).gt(0)
                    
                    task1 = ee.batch.Export.image.toAsset(image=Mmask_C2M_k1,
                                                         description='Mmask_C2_calib_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr)+'_var_'+hydro+'_k'+str(kern1),
                                                         region=region,
                                                         assetId= folderGEE+'/Mmask_C2_calib_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr)+'_var_'+hydro+'_k'+str(kern1),
                                                         crs=crs_mod,
                                                         crsTransform=trans_mod)
                    task1.start()
                    
            if pix_Mcal_k2_CM==1: 
                # calculate calibrated CM M k2 pixel
                print('M CM k2 cal extraction')                     
                try:
                    Mmask_C2M_k2=ee.Image(folderGEE+'/Mmask_C2_calib_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr)+'_var_'+hydro+'_k'+str(kern2)).reproject(crs_mod,trans_mod)
                    Mmask_C2M_k2.getInfo()
               
                except: 
                    ee.Reset()
                    ee.Initialize(credentials)
                    def calc_points(image):
                        value =image.select(band).updateMask(C2mask).reduceRegion(**{'reducer':ee.Reducer.mean(),
                                                                  'geometry': region,
                                                                  'crs': crs_mod,
                                                                  'crsTransform':trans_mod,
                                                                  'bestEffort': True});
                        value=ee.Number(ee.List([value.get(value.keys().get(0)), -9999]).reduce(ee.Reducer.firstNonNull()));
                        im = im01.unmask(value).updateMask(value.neq(-9999));
                        return image.addBands(im);
                    collection0 = collection.map(calc_points);
                
                    coll02=collection0.map(lambda image: image.select(['Q','C']).addBands(image.select(band).reduceNeighborhood(ee.Reducer.mean(), ee.Kernel.square(kern2,'pixels',True)).rename(band).reproject(crs_mod,trans_mod)));            
                    CM=coll02.map(lambda image: image.select('Q').addBands(image.select('C').divide(image.select(band))));
               
                    # calculate correlation between discharge and C/M
                    corr=CM.reduce(ee.Reducer.spearmansCorrelation()).reproject(crs_mod,trans_mod).updateMask(Wat2);
                    
                    # find maximum correlation Calibration and create mask for that pixel
                    maxval = corr.reduceRegion(**{'reducer':ee.Reducer.max(), 'geometry':region, 'crs':crs_mod,'crsTransform':trans_mod, 'bestEffort':True});
                    Mmask_C2M_k2 = corr.select('correlation').eq(ee.Number(maxval.get('correlation')));
                    Mmask_C2M_k2 = Mmask_C2M_k2.reduceNeighborhood(ee.Reducer.mean(), ee.Kernel.square(kern2,'pixels',True)).rename('correlation').reproject(crs_mod,trans_mod).gt(0)

                    task2 = ee.batch.Export.image.toAsset(image=Mmask_C2M_k2,
                                                         description='Mmask_C2_calib_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr)+'_var_'+hydro+'_k'+str(kern2),
                                                         region=region,
                                                         assetId= folderGEE+'/Mmask_C2_calib_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr)+'_var_'+hydro+'_k'+str(kern2),
                                                         crs=crs_mod,
                                                         crsTransform=trans_mod)
                    task2.start()
                    
            if pix_Mcal_k1_CMW==1 and pix_W1==1: 
                # calculate calibrated CMW1 M k1 pixel
                print('M CMW1 k1 cal extraction')                     
                try:
                    Mmask_C2W1M_k1=ee.Image(folderGEE+'/Mmask_C2_calib_W1_VD_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr)+'_var_'+hydro+'_k'+str(kern1)).reproject(crs_mod,trans_mod)
                    Mmask_C2W1M_k1.getInfo()
               
                except: 
                    ee.Reset()
                    ee.Initialize(credentials)
                    def calc_points(image):
                        value =image.select(band).updateMask(C2mask).reduceRegion(**{'reducer':ee.Reducer.mean(),
                                                                  'geometry': region,
                                                                  'crs': crs_mod,
                                                                  'crsTransform':trans_mod,
                                                                  'bestEffort': True});
                        value=ee.Number(ee.List([value.get(value.keys().get(0)), -9999]).reduce(ee.Reducer.firstNonNull()));
                        im = im01.unmask(value).updateMask(value.neq(-9999));
                        value2 =image.select(band).updateMask(Wmask).reduceRegion(**{'reducer':ee.Reducer.mean(),
                                                                  'geometry': region,
                                                                  'crs': crs_mod,
                                                                  'crsTransform':trans_mod,
                                                                  'bestEffort': True});
                        value2=ee.Number(ee.List([value2.get(value2.keys().get(0)), -9999]).reduce(ee.Reducer.firstNonNull()));
                        im2 = im02.unmask(value2).updateMask(value2.neq(-9999));
                        return image.addBands([im,im2]);
                    collection0 = collection.map(calc_points);
                
                    coll02=collection0.map(lambda image: image.select(['W','Q','C']).addBands(image.select(band).reduceNeighborhood(ee.Reducer.mean(), ee.Kernel.square(kern1,'pixels',True)).rename(band).reproject(crs_mod,trans_mod)));            
                    def calc_VD(image):
                        im=image.select('C').subtract(image.select(band)).divide(image.select('C').subtract(image.select('W')))
                        im=im.updateMask(im.gte(3).Or(im.lte(-2)).Not())
                        im=im.where(im.gte(1),1)
                        im=im.where(im.lte(0),0).rename(['VD'])
                        return image.addBands(im)
                        
                    coll02=coll02.map(calc_VD)
            
                    coeff=coll02.map(lambda image: image.select('W').multiply(image.select('VD')).subtract(image.select(band))).max().add(coll02.select(band).min())
                    CMW=coll02.map(lambda image: image.select('Q').addBands(image.select('C').divide(image.select(band).subtract(image.select('W').multiply(image.select('VD'))).add(coeff))));
               
                    
                    # calculate correlation between discharge and C/M
                    corr=CMW.reduce(ee.Reducer.spearmansCorrelation()).reproject(crs_mod,trans_mod).updateMask(Wat2);
                    
                    # find maximum correlation Calibration and create mask for that pixel
                    maxval = corr.reduceRegion(**{'reducer':ee.Reducer.max(), 'geometry':region, 'crs':crs_mod,'crsTransform':trans_mod, 'bestEffort':True});
                    Mmask_C2W1M_k1 = corr.select('correlation').eq(ee.Number(maxval.get('correlation')));
                    Mmask_C2W1M_k1 = Mmask_C2W1M_k1.reduceNeighborhood(ee.Reducer.mean(), ee.Kernel.square(kern1,'pixels',True)).rename('correlation').reproject(crs_mod,trans_mod).gt(0)

                    task3 = ee.batch.Export.image.toAsset(image=Mmask_C2W1M_k1,
                                                         description='Mmask_C2_calib_W1_VD_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr)+'_var_'+hydro+'_k'+str(kern1),
                                                         region=region,
                                                         assetId= folderGEE+'/Mmask_C2_calib_W1_VD_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr)+'_var_'+hydro+'_k'+str(kern1),
                                                         crs=crs_mod,
                                                         crsTransform=trans_mod)
                    task3.start()
                    
            if pix_Mcal_k1_CMW==1 and pix_W2==1: 
                # calculate calibrated CMW2 M k1 pixel
                print('M CMW2 k1 cal extraction')                     
                try:
                    Mmask_C2W2M_k1=ee.Image(folderGEE+'/Mmask_C2_calib_W2_VD_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr)+'_var_'+hydro+'_k'+str(kern1)).reproject(crs_mod,trans_mod)
                    Mmask_C2W2M_k1.getInfo()
               
                except: 
                    ee.Reset()
                    ee.Initialize(credentials)
                    def calc_points(image):
                        value =image.select(band).updateMask(C2mask).reduceRegion(**{'reducer':ee.Reducer.mean(),
                                                                  'geometry': region,
                                                                  'crs': crs_mod,
                                                                  'crsTransform':trans_mod,
                                                                  'bestEffort': True});
                        value=ee.Number(ee.List([value.get(value.keys().get(0)), -9999]).reduce(ee.Reducer.firstNonNull()));
                        im = im01.unmask(value).updateMask(value.neq(-9999));
                        value2 =image.select(band).updateMask(Wmask2).reduceRegion(**{'reducer':ee.Reducer.mean(),
                                                                  'geometry': region,
                                                                  'crs': crs_mod,
                                                                  'crsTransform':trans_mod,
                                                                  'bestEffort': True});
                        value2=ee.Number(ee.List([value2.get(value2.keys().get(0)), -9999]).reduce(ee.Reducer.firstNonNull()));
                        im2 = im02.unmask(value2).updateMask(value2.neq(-9999));
                        return image.addBands([im,im2]);
                    collection0 = collection.map(calc_points);
                
                    coll02=collection0.map(lambda image: image.select(['W','Q','C']).addBands(image.select(band).reduceNeighborhood(ee.Reducer.mean(), ee.Kernel.square(kern1,'pixels',True)).rename(band).reproject(crs_mod,trans_mod)));            
                    def calc_VD(image):
                        im=image.select('C').subtract(image.select(band)).divide(image.select('C').subtract(image.select('W')))
                        im=im.updateMask(im.gte(3).Or(im.lte(-2)).Not())
                        im=im.where(im.gte(1),1)
                        im=im.where(im.lte(0),0).rename(['VD'])
                        return image.addBands(im)
                        
                    coll02=coll02.map(calc_VD)
            
                    coeff=coll02.map(lambda image: image.select('W').multiply(image.select('VD')).subtract(image.select(band))).max().add(coll02.select(band).min())
                    CMW=coll02.map(lambda image: image.select('Q').addBands(image.select('C').divide(image.select(band).subtract(image.select('W').multiply(image.select('VD'))).add(coeff))));
               
                    
                    # calculate correlation between discharge and C/M
                    corr=CMW.reduce(ee.Reducer.spearmansCorrelation()).reproject(crs_mod,trans_mod).updateMask(Wat2);
                    
                    # find maximum correlation Calibration and create mask for that pixel
                    maxval = corr.reduceRegion(**{'reducer':ee.Reducer.max(), 'geometry':region, 'crs':crs_mod,'crsTransform':trans_mod, 'bestEffort':True});
                    Mmask_C2W2M_k1 = corr.select('correlation').eq(ee.Number(maxval.get('correlation')));
                    Mmask_C2W2M_k1 = Mmask_C2W2M_k1.reduceNeighborhood(ee.Reducer.mean(), ee.Kernel.square(kern1,'pixels',True)).rename('correlation').reproject(crs_mod,trans_mod).gt(0)

                    task4 = ee.batch.Export.image.toAsset(image=Mmask_C2W2M_k1,
                                                         description='Mmask_C2_calib_W2_VD_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr)+'_var_'+hydro+'_k'+str(kern1),
                                                         region=region,
                                                         assetId= folderGEE+'/Mmask_C2_calib_W2_VD_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr)+'_var_'+hydro+'_k'+str(kern1),
                                                         crs=crs_mod,
                                                         crsTransform=trans_mod)
                    task4.start()
                                    
            if pix_Mcal_k2_CMW==1 and pix_W1==1: 
                # calculate calibrated CMW1 M k2 pixel
                print('M CMW1 k2 cal extraction')                     
                try:
                    Mmask_C2W1M_k2=ee.Image(folderGEE+'/Mmask_C2_calib_W1_VD_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr)+'_var_'+hydro+'_k'+str(kern2)).reproject(crs_mod,trans_mod)
                    Mmask_C2W1M_k2.getInfo()
               
                except: 
                    ee.Reset()
                    ee.Initialize(credentials)
                    def calc_points(image):
                        value =image.select(band).updateMask(C2mask).reduceRegion(**{'reducer':ee.Reducer.mean(),
                                                                  'geometry': region,
                                                                  'crs': crs_mod,
                                                                  'crsTransform':trans_mod,
                                                                  'bestEffort': True});
                        value=ee.Number(ee.List([value.get(value.keys().get(0)), -9999]).reduce(ee.Reducer.firstNonNull()));
                        im = im01.unmask(value).updateMask(value.neq(-9999));
                        value2 =image.select(band).updateMask(Wmask).reduceRegion(**{'reducer':ee.Reducer.mean(),
                                                                  'geometry': region,
                                                                  'crs': crs_mod,
                                                                  'crsTransform':trans_mod,
                                                                  'bestEffort': True});
                        value2=ee.Number(ee.List([value2.get(value2.keys().get(0)), -9999]).reduce(ee.Reducer.firstNonNull()));
                        im2 = im02.unmask(value2).updateMask(value2.neq(-9999));
                        return image.addBands([im,im2]);
                    collection0 = collection.map(calc_points);
                
                    coll02=collection0.map(lambda image: image.select(['W','Q','C']).addBands(image.select(band).reduceNeighborhood(ee.Reducer.mean(), ee.Kernel.square(kern2,'pixels',True)).rename(band).reproject(crs_mod,trans_mod)));            
                    def calc_VD(image):
                        im=image.select('C').subtract(image.select(band)).divide(image.select('C').subtract(image.select('W')))
                        im=im.updateMask(im.gte(3).Or(im.lte(-2)).Not())
                        im=im.where(im.gte(1),1)
                        im=im.where(im.lte(0),0).rename(['VD'])
                        return image.addBands(im)
                        
                    coll02=coll02.map(calc_VD)
            
                    coeff=coll02.map(lambda image: image.select('W').multiply(image.select('VD')).subtract(image.select(band))).max().add(coll02.select(band).min())
                    CMW=coll02.map(lambda image: image.select('Q').addBands(image.select('C').divide(image.select(band).subtract(image.select('W').multiply(image.select('VD'))).add(coeff))));
               
                    
                    # calculate correlation between discharge and C/M
                    corr=CMW.reduce(ee.Reducer.spearmansCorrelation()).reproject(crs_mod,trans_mod).updateMask(Wat2);
                    
                    # find maximum correlation Calibration and create mask for that pixel
                    maxval = corr.reduceRegion(**{'reducer':ee.Reducer.max(), 'geometry':region, 'crs':crs_mod,'crsTransform':trans_mod, 'bestEffort':True});
                    Mmask_C2W1M_k2 = corr.select('correlation').eq(ee.Number(maxval.get('correlation')));
                    Mmask_C2W1M_k2 = Mmask_C2W1M_k2.reduceNeighborhood(ee.Reducer.mean(), ee.Kernel.square(kern2,'pixels',True)).rename('correlation').reproject(crs_mod,trans_mod).gt(0)

                    task5 = ee.batch.Export.image.toAsset(image=Mmask_C2W1M_k2,
                                                         description='Mmask_C2_calib_W1_VD_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr)+'_var_'+hydro+'_k'+str(kern2),
                                                         region=region,
                                                         assetId= folderGEE+'/Mmask_C2_calib_W1_VD_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr)+'_var_'+hydro+'_k'+str(kern2),
                                                         crs=crs_mod,
                                                         crsTransform=trans_mod)
                    task5.start()
           
            if pix_Mcal_k2_CMW==1 and pix_W2==1: 
                # calculate calibrated CMW2 M k2 pixel
                print('M CMW2 k2 cal extraction')                     
                try:
                    Mmask_C2W2M_k2=ee.Image(folderGEE+'/Mmask_C2_calib_W2_VD_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr)+'_var_'+hydro+'_k'+str(kern2)).reproject(crs_mod,trans_mod)
                    Mmask_C2W2M_k2.getInfo()
               
                except: 
                    ee.Reset()
                    ee.Initialize(credentials)
                    def calc_points(image):
                        value =image.select(band).updateMask(C2mask).reduceRegion(**{'reducer':ee.Reducer.mean(),
                                                                  'geometry': region,
                                                                  'crs': crs_mod,
                                                                  'crsTransform':trans_mod,
                                                                  'bestEffort': True});
                        value=ee.Number(ee.List([value.get(value.keys().get(0)), -9999]).reduce(ee.Reducer.firstNonNull()));
                        im = im01.unmask(value).updateMask(value.neq(-9999));
                        value2 =image.select(band).updateMask(Wmask2).reduceRegion(**{'reducer':ee.Reducer.mean(),
                                                                  'geometry': region,
                                                                  'crs': crs_mod,
                                                                  'crsTransform':trans_mod,
                                                                  'bestEffort': True});
                        value2=ee.Number(ee.List([value2.get(value2.keys().get(0)), -9999]).reduce(ee.Reducer.firstNonNull()));
                        im2 = im02.unmask(value2).updateMask(value2.neq(-9999));
                        return image.addBands([im,im2]);
                    collection0 = collection.map(calc_points);
                
                    coll02=collection0.map(lambda image: image.select(['W','Q','C']).addBands(image.select(band).reduceNeighborhood(ee.Reducer.mean(), ee.Kernel.square(kern2,'pixels',True)).rename(band).reproject(crs_mod,trans_mod)));            
                    def calc_VD(image):
                        im=image.select('C').subtract(image.select(band)).divide(image.select('C').subtract(image.select('W')))
                        im=im.updateMask(im.gte(3).Or(im.lte(-2)).Not())
                        im=im.where(im.gte(1),1)
                        im=im.where(im.lte(0),0).rename(['VD'])
                        return image.addBands(im)
                        
                    coll02=coll02.map(calc_VD)
            
                    coeff=coll02.map(lambda image: image.select('W').multiply(image.select('VD')).subtract(image.select(band))).max().add(coll02.select(band).min())
                    CMW=coll02.map(lambda image: image.select('Q').addBands(image.select('C').divide(image.select(band).subtract(image.select('W').multiply(image.select('VD'))).add(coeff))));
               
                    
                    # calculate correlation between discharge and C/M
                    corr=CMW.reduce(ee.Reducer.spearmansCorrelation()).reproject(crs_mod,trans_mod).updateMask(Wat2);
                        
                    # find maximum correlation Calibration and create mask for that pixel
                    maxval = corr.reduceRegion(**{'reducer':ee.Reducer.max(), 'geometry':region, 'crs':crs_mod,'crsTransform':trans_mod, 'bestEffort':True});
                    Mmask_C2W2M_k2 = corr.select('correlation').eq(ee.Number(maxval.get('correlation')));
                    Mmask_C2W2M_k2 = Mmask_C2W2M_k2.reduceNeighborhood(ee.Reducer.mean(), ee.Kernel.square(kern2,'pixels',True)).rename('correlation').reproject(crs_mod,trans_mod).gt(0)

                    task6 = ee.batch.Export.image.toAsset(image=Mmask_C2W2M_k2,
                                                         description='Mmask_C2_calib_W2_VD_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr)+'_var_'+hydro+'_k'+str(kern2),
                                                         region=region,
                                                         assetId= folderGEE+'/Mmask_C2_calib_W2_VD_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr)+'_var_'+hydro+'_k'+str(kern2),
                                                         crs=crs_mod,
                                                         crsTransform=trans_mod)
                    task6.start()
                        
        if veg==1:                
            if pix_Mcal_k1_CM==1: 
                # calculate calibrated CM M k1 pixel
                print('M CVM k1 cal extraction')                     
                try:
                    Mmask_C2VM_k1=ee.Image(folderGEE+'/Mmask_C2V_calib_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr)+'_var_'+hydro+'_k'+str(kern1)).reproject(crs_mod,trans_mod)
                    Mmask_C2VM_k1.getInfo()
               
                except: 
                    ee.Reset()
                    ee.Initialize(credentials)
                    def calc_points(image):
                        value =image.select(band).updateMask(C2mask).reduceRegion(**{'reducer':ee.Reducer.mean(),
                                                                  'geometry': region,
                                                                  'crs': crs_mod,
                                                                  'crsTransform':trans_mod,
                                                                  'bestEffort': True});
                        value=ee.Number(ee.List([value.get(value.keys().get(0)), -9999]).reduce(ee.Reducer.firstNonNull()));
                        value2 =image.select(band).updateMask(Vmask).reduceRegion(**{'reducer':ee.Reducer.mean(),
                                                                  'geometry': region,
                                                                  'crs': crs_mod,
                                                                  'crsTransform':trans_mod,
                                                                  'bestEffort': True});
                        value2=ee.Number(ee.List([value2.get(value2.keys().get(0)), -9999]).reduce(ee.Reducer.firstNonNull()));
                        im = ee.ImageCollection([im01.unmask(value).updateMask(value.neq(-9999)),im01.unmask(value2).updateMask(value2.neq(-9999))]).mean().rename('C')
                        return image.addBands(im);
                    collection0 = collection.map(calc_points);
                
                    coll02=collection0.map(lambda image: image.select(['Q','C']).addBands(image.select(band).reduceNeighborhood(ee.Reducer.mean(), ee.Kernel.square(kern1,'pixels',True)).rename(band).reproject(crs_mod,trans_mod)));            
                    CVM=coll02.map(lambda image: image.select('Q').addBands(image.select('C').divide(image.select(band))));
            
                    # calculate correlation between discharge and C/M
                    corr=CVM.reduce(ee.Reducer.spearmansCorrelation()).reproject(crs_mod,trans_mod).updateMask(Wat2);
                    
                    # find maximum correlation Calibration and create mask for that pixel
                    maxval = corr.reduceRegion(**{'reducer':ee.Reducer.max(), 'geometry':region, 'crs':crs_mod,'crsTransform':trans_mod, 'bestEffort':True});
                    Mmask_C2VM_k1 = corr.select('correlation').eq(ee.Number(maxval.get('correlation')));
                    Mmask_C2VM_k1 = Mmask_C2VM_k1.reduceNeighborhood(ee.Reducer.mean(), ee.Kernel.square(kern1,'pixels',True)).rename('correlation').reproject(crs_mod,trans_mod).gt(0)

                    task7 = ee.batch.Export.image.toAsset(image=Mmask_C2VM_k1,
                                                         description='Mmask_C2V_calib_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr)+'_var_'+hydro+'_k'+str(kern1),
                                                         region=region,
                                                         assetId= folderGEE+'/Mmask_C2V_calib_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr)+'_var_'+hydro+'_k'+str(kern1),
                                                         crs=crs_mod,
                                                         crsTransform=trans_mod)
                    task7.start()

            if pix_Mcal_k2_CM==1: 
                # calculate calibrated CM M k2 pixel
                print('M CVM k2 cal extraction')                     
                try:
                    Mmask_C2VM_k2=ee.Image(folderGEE+'/Mmask_C2V_calib_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr)+'_var_'+hydro+'_k'+str(kern2)).reproject(crs_mod,trans_mod)
                    Mmask_C2VM_k2.getInfo()
               
                except: 
                    ee.Reset()
                    ee.Initialize(credentials)
                    def calc_points(image):
                        value =image.select(band).updateMask(C2mask).reduceRegion(**{'reducer':ee.Reducer.mean(),
                                                                  'geometry': region,
                                                                  'crs': crs_mod,
                                                                  'crsTransform':trans_mod,
                                                                  'bestEffort': True});
                        value=ee.Number(ee.List([value.get(value.keys().get(0)), -9999]).reduce(ee.Reducer.firstNonNull()));
                        value2 =image.select(band).updateMask(Vmask).reduceRegion(**{'reducer':ee.Reducer.mean(),
                                                                  'geometry': region,
                                                                  'crs': crs_mod,
                                                                  'crsTransform':trans_mod,
                                                                  'bestEffort': True});
                        value2=ee.Number(ee.List([value2.get(value2.keys().get(0)), -9999]).reduce(ee.Reducer.firstNonNull()));
                        im = ee.ImageCollection([im01.unmask(value).updateMask(value.neq(-9999)),im01.unmask(value2).updateMask(value2.neq(-9999))]).mean().rename('C')
                        return image.addBands(im);
                    collection0 = collection.map(calc_points);
                
                    coll02=collection0.map(lambda image: image.select(['Q','C']).addBands(image.select(band).reduceNeighborhood(ee.Reducer.mean(), ee.Kernel.square(kern2,'pixels',True)).rename(band).reproject(crs_mod,trans_mod)));            
                    CVM=coll02.map(lambda image: image.select('Q').addBands(image.select('C').divide(image.select(band))));
               
                    # calculate correlation between discharge and C/M
                    corr=CVM.reduce(ee.Reducer.spearmansCorrelation()).reproject(crs_mod,trans_mod).updateMask(Wat2);
                
                    # find maximum correlation Calibration and create mask for that pixel
                    maxval = corr.reduceRegion(**{'reducer':ee.Reducer.max(), 'geometry':region, 'crs':crs_mod,'crsTransform':trans_mod, 'bestEffort':True});
                    Mmask_C2VM_k2 = corr.select('correlation').eq(ee.Number(maxval.get('correlation')));
                    Mmask_C2VM_k2 = Mmask_C2VM_k2.reduceNeighborhood(ee.Reducer.mean(), ee.Kernel.square(kern2,'pixels',True)).rename('correlation').reproject(crs_mod,trans_mod).gt(0)
                    
                    task8 = ee.batch.Export.image.toAsset(image=Mmask_C2VM_k2,
                                                         description='Mmask_C2V_calib_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr)+'_var_'+hydro+'_k'+str(kern2),
                                                         region=region,
                                                         assetId= folderGEE+'/Mmask_C2V_calib_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr)+'_var_'+hydro+'_k'+str(kern2),
                                                         crs=crs_mod,
                                                         crsTransform=trans_mod)
                    task8.start()
                                           
            if pix_Mcal_k1_CMW==1 and pix_W1==1: 
                # calculate calibrated CMW1 M k1 pixel
                print('M CVMW1 k1 cal extraction')                     
                try:
                    Mmask_C2VW1M_k1=ee.Image(folderGEE+'/Mmask_C2V_calib_W1_VD_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr)+'_var_'+hydro+'_k'+str(kern1)).reproject(crs_mod,trans_mod)
                    Mmask_C2VW1M_k1.getInfo()
               
                except: 
                    ee.Reset()
                    ee.Initialize(credentials)
                    def calc_points(image):
                        value =image.select(band).updateMask(C2mask).reduceRegion(**{'reducer':ee.Reducer.mean(),
                                                                  'geometry': region,
                                                                  'crs': crs_mod,
                                                                  'crsTransform':trans_mod,
                                                                  'bestEffort': True});
                        value=ee.Number(ee.List([value.get(value.keys().get(0)), -9999]).reduce(ee.Reducer.firstNonNull()));
                        value2 =image.select(band).updateMask(Vmask).reduceRegion(**{'reducer':ee.Reducer.mean(),
                                                                  'geometry': region,
                                                                  'crs': crs_mod,
                                                                  'crsTransform':trans_mod,
                                                                  'bestEffort': True});
                        value2=ee.Number(ee.List([value2.get(value2.keys().get(0)), -9999]).reduce(ee.Reducer.firstNonNull()));
                        im = ee.ImageCollection([im01.unmask(value).updateMask(value.neq(-9999)),im01.unmask(value2).updateMask(value2.neq(-9999))]).mean().rename('C')
                        
                        value2 =image.select(band).updateMask(Wmask).reduceRegion(**{'reducer':ee.Reducer.mean(),
                                                                  'geometry': region,
                                                                  'crs': crs_mod,
                                                                  'crsTransform':trans_mod,
                                                                  'bestEffort': True});
                        value2=ee.Number(ee.List([value2.get(value2.keys().get(0)), -9999]).reduce(ee.Reducer.firstNonNull()));
                        im2 = im02.unmask(value2).updateMask(value2.neq(-9999));
                        return image.addBands([im,im2]);
                    collection0 = collection.map(calc_points);
                
                    coll02=collection0.map(lambda image: image.select(['W','Q','C']).addBands(image.select(band).reduceNeighborhood(ee.Reducer.mean(), ee.Kernel.square(kern1,'pixels',True)).rename(band).reproject(crs_mod,trans_mod)));            
                    def calc_VD(image):
                        im=image.select('C').subtract(image.select(band)).divide(image.select('C').subtract(image.select('W')))
                        im=im.updateMask(im.gte(3).Or(im.lte(-2)).Not())
                        im=im.where(im.gte(1),1)
                        im=im.where(im.lte(0),0).rename(['VD'])
                        return image.addBands(im)
                        
                    coll02=coll02.map(calc_VD)
            
                    coeff=coll02.map(lambda image: image.select('W').multiply(image.select('VD')).subtract(image.select(band))).max().add(coll02.select(band).min())
                    CVMW=coll02.map(lambda image: image.select('Q').addBands(image.select('C').divide(image.select(band).subtract(image.select('W').multiply(image.select('VD'))).add(coeff))));
               
                    
                    # calculate correlation between discharge and C/M
                    corr=CVMW.reduce(ee.Reducer.spearmansCorrelation()).reproject(crs_mod,trans_mod).updateMask(Wat2);
                    
                    # find maximum correlation Calibration and create mask for that pixel
                    maxval = corr.reduceRegion(**{'reducer':ee.Reducer.max(), 'geometry':region, 'crs':crs_mod,'crsTransform':trans_mod, 'bestEffort':True});
                    Mmask_C2VW1M_k1 = corr.select('correlation').eq(ee.Number(maxval.get('correlation')));
                    Mmask_C2VW1M_k1 = Mmask_C2VW1M_k1.reduceNeighborhood(ee.Reducer.mean(), ee.Kernel.square(kern1,'pixels',True)).rename('correlation').reproject(crs_mod,trans_mod).gt(0)

                    task9 = ee.batch.Export.image.toAsset(image=Mmask_C2VW1M_k1,
                                                         description='Mmask_C2V_calib_W1_VD_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr)+'_var_'+hydro+'_k'+str(kern1),
                                                         region=region,
                                                         assetId= folderGEE+'/Mmask_C2V_calib_W1_VD_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr)+'_var_'+hydro+'_k'+str(kern1),
                                                         crs=crs_mod,
                                                         crsTransform=trans_mod)
                    task9.start()

            if pix_Mcal_k1_CMW==1 and pix_W2==1: 
                # calculate calibrated CMW2 M k1 pixel
                print('M CMW2 k1 cal extraction')                     
                try:
                    Mmask_C2VW2M_k1=ee.Image(folderGEE+'/Mmask_C2V_calib_W2_VD_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr)+'_var_'+hydro+'_k'+str(kern1)).reproject(crs_mod,trans_mod)
                    Mmask_C2VW2M_k1.getInfo()
               
                except: 
                    ee.Reset()
                    ee.Initialize(credentials)
                    def calc_points(image):
                        value =image.select(band).updateMask(C2mask).reduceRegion(**{'reducer':ee.Reducer.mean(),
                                                                  'geometry': region,
                                                                  'crs': crs_mod,
                                                                  'crsTransform':trans_mod,
                                                                  'bestEffort': True});
                        value=ee.Number(ee.List([value.get(value.keys().get(0)), -9999]).reduce(ee.Reducer.firstNonNull()));
                        value2 =image.select(band).updateMask(Vmask).reduceRegion(**{'reducer':ee.Reducer.mean(),
                                                                  'geometry': region,
                                                                  'crs': crs_mod,
                                                                  'crsTransform':trans_mod,
                                                                  'bestEffort': True});
                        value2=ee.Number(ee.List([value2.get(value2.keys().get(0)), -9999]).reduce(ee.Reducer.firstNonNull()));
                        im = ee.ImageCollection([im01.unmask(value).updateMask(value.neq(-9999)),im01.unmask(value2).updateMask(value2.neq(-9999))]).mean().rename('C')

                        value2 =image.select(band).updateMask(Wmask2).reduceRegion(**{'reducer':ee.Reducer.mean(),
                                                                  'geometry': region,
                                                                  'crs': crs_mod,
                                                                  'crsTransform':trans_mod,
                                                                  'bestEffort': True});
                        value2=ee.Number(ee.List([value2.get(value2.keys().get(0)), -9999]).reduce(ee.Reducer.firstNonNull()));
                        im2 = im02.unmask(value2).updateMask(value2.neq(-9999));
                        return image.addBands([im,im2]);
                    collection0 = collection.map(calc_points);
                
                    coll02=collection0.map(lambda image: image.select(['W','Q','C']).addBands(image.select(band).reduceNeighborhood(ee.Reducer.mean(), ee.Kernel.square(kern1,'pixels',True)).rename(band).reproject(crs_mod,trans_mod)));            
                    def calc_VD(image):
                        im=image.select('C').subtract(image.select(band)).divide(image.select('C').subtract(image.select('W')))
                        im=im.updateMask(im.gte(3).Or(im.lte(-2)).Not())
                        im=im.where(im.gte(1),1)
                        im=im.where(im.lte(0),0).rename(['VD'])
                        return image.addBands(im)
                        
                    coll02=coll02.map(calc_VD)
            
                    coeff=coll02.map(lambda image: image.select('W').multiply(image.select('VD')).subtract(image.select(band))).max().add(coll02.select(band).min())
                    CVMW=coll02.map(lambda image: image.select('Q').addBands(image.select('C').divide(image.select(band).subtract(image.select('W').multiply(image.select('VD'))).add(coeff))));
               
                    
                    # calculate correlation between discharge and C/M
                    corr=CVMW.reduce(ee.Reducer.spearmansCorrelation()).reproject(crs_mod,trans_mod).updateMask(Wat2);
                    
                    # find maximum correlation Calibration and create mask for that pixel
                    maxval = corr.reduceRegion(**{'reducer':ee.Reducer.max(), 'geometry':region, 'crs':crs_mod,'crsTransform':trans_mod, 'bestEffort':True});
                    Mmask_C2VW2M_k1 = corr.select('correlation').eq(ee.Number(maxval.get('correlation')));
                    Mmask_C2VW2M_k1 = Mmask_C2VW2M_k1.reduceNeighborhood(ee.Reducer.mean(), ee.Kernel.square(kern1,'pixels',True)).rename('correlation').reproject(crs_mod,trans_mod).gt(0)

                    task10 = ee.batch.Export.image.toAsset(image=Mmask_C2VW2M_k1,
                                                         description='Mmask_C2V_calib_W2_VD_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr)+'_var_'+hydro+'_k'+str(kern1),
                                                         region=region,
                                                         assetId= folderGEE+'/Mmask_C2V_calib_W2_VD_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr)+'_var_'+hydro+'_k'+str(kern1),
                                                         crs=crs_mod,
                                                         crsTransform=trans_mod)
                    task10.start()
                                        
            if pix_Mcal_k2_CMW==1 and pix_W1==1: 
                # calculate calibrated CMW1 M k2 pixel
                print('M CVMW1 k2 cal extraction')                     
                try:
                    Mmask_C2VW1M_k2=ee.Image(folderGEE+'/Mmask_C2V_calib_W1_VD_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr)+'_var_'+hydro+'_k'+str(kern2)).reproject(crs_mod,trans_mod)
                    Mmask_C2VW1M_k2.getInfo()
               
                except: 
                    ee.Reset()
                    ee.Initialize(credentials)
                    def calc_points(image):
                        value =image.select(band).updateMask(C2mask).reduceRegion(**{'reducer':ee.Reducer.mean(),
                                                                  'geometry': region,
                                                                  'crs': crs_mod,
                                                                  'crsTransform':trans_mod,
                                                                  'bestEffort': True});
                        value=ee.Number(ee.List([value.get(value.keys().get(0)), -9999]).reduce(ee.Reducer.firstNonNull()));
                        value2 =image.select(band).updateMask(Vmask).reduceRegion(**{'reducer':ee.Reducer.mean(),
                                                                  'geometry': region,
                                                                  'crs': crs_mod,
                                                                  'crsTransform':trans_mod,
                                                                  'bestEffort': True});
                        value2=ee.Number(ee.List([value2.get(value2.keys().get(0)), -9999]).reduce(ee.Reducer.firstNonNull()));
                        im = ee.ImageCollection([im01.unmask(value).updateMask(value.neq(-9999)),im01.unmask(value2).updateMask(value2.neq(-9999))]).mean().rename('C')

                        value2 =image.select(band).updateMask(Wmask).reduceRegion(**{'reducer':ee.Reducer.mean(),
                                                                  'geometry': region,
                                                                  'crs': crs_mod,
                                                                  'crsTransform':trans_mod,
                                                                  'bestEffort': True});
                        value2=ee.Number(ee.List([value2.get(value2.keys().get(0)), -9999]).reduce(ee.Reducer.firstNonNull()));
                        im2 = im02.unmask(value2).updateMask(value2.neq(-9999));
                        return image.addBands([im,im2]);
                    collection0 = collection.map(calc_points);
                
                    coll02=collection0.map(lambda image: image.select(['W','Q','C']).addBands(image.select(band).reduceNeighborhood(ee.Reducer.mean(), ee.Kernel.square(kern2,'pixels',True)).rename(band).reproject(crs_mod,trans_mod)));            
                    def calc_VD(image):
                        im=image.select('C').subtract(image.select(band)).divide(image.select('C').subtract(image.select('W')))
                        im=im.updateMask(im.gte(3).Or(im.lte(-2)).Not())
                        im=im.where(im.gte(1),1)
                        im=im.where(im.lte(0),0).rename(['VD'])
                        return image.addBands(im)
                        
                    coll02=coll02.map(calc_VD)
            
                    coeff=coll02.map(lambda image: image.select('W').multiply(image.select('VD')).subtract(image.select(band))).max().add(coll02.select(band).min())
                    CVMW=coll02.map(lambda image: image.select('Q').addBands(image.select('C').divide(image.select(band).subtract(image.select('W').multiply(image.select('VD'))).add(coeff))));
               
                    
                    # calculate correlation between discharge and C/M
                    corr=CVMW.reduce(ee.Reducer.spearmansCorrelation()).reproject(crs_mod,trans_mod).updateMask(Wat2);
                    
                    # find maximum correlation Calibration and create mask for that pixel
                    maxval = corr.reduceRegion(**{'reducer':ee.Reducer.max(), 'geometry':region, 'crs':crs_mod,'crsTransform':trans_mod, 'bestEffort':True});
                    Mmask_C2VW1M_k2 = corr.select('correlation').eq(ee.Number(maxval.get('correlation')));
                    Mmask_C2VW1M_k2 = Mmask_C2VW1M_k2.reduceNeighborhood(ee.Reducer.mean(), ee.Kernel.square(kern2,'pixels',True)).rename('correlation').reproject(crs_mod,trans_mod).gt(0)
                    
                    task11 = ee.batch.Export.image.toAsset(image=Mmask_C2VW1M_k2,
                                                         description='Mmask_C2V_calib_W1_VD_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr)+'_var_'+hydro+'_k'+str(kern2),
                                                         region=region,
                                                         assetId= folderGEE+'/Mmask_C2V_calib_W1_VD_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr)+'_var_'+hydro+'_k'+str(kern2),
                                                         crs=crs_mod,
                                                         crsTransform=trans_mod)
                    task11.start()
           
            if pix_Mcal_k2_CMW==1 and pix_W2==1: 
                # calculate calibrated CMW2 M k2 pixel
                print('M CVMW2 k2 cal extraction')                     
                try:
                    Mmask_C2VW2M_k2=ee.Image(folderGEE+'/Mmask_C2V_calib_W2_VD_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr)+'_var_'+hydro+'_k'+str(kern2)).reproject(crs_mod,trans_mod)
                    Mmask_C2VW2M_k2.getInfo()
               
                except: 
                    ee.Reset()
                    ee.Initialize(credentials)
                    def calc_points(image):
                        value =image.select(band).updateMask(C2mask).reduceRegion(**{'reducer':ee.Reducer.mean(),
                                                                  'geometry': region,
                                                                  'crs': crs_mod,
                                                                  'crsTransform':trans_mod,
                                                                  'bestEffort': True});
                        value=ee.Number(ee.List([value.get(value.keys().get(0)), -9999]).reduce(ee.Reducer.firstNonNull()));
                        value2 =image.select(band).updateMask(Vmask).reduceRegion(**{'reducer':ee.Reducer.mean(),
                                                                  'geometry': region,
                                                                  'crs': crs_mod,
                                                                  'crsTransform':trans_mod,
                                                                  'bestEffort': True});
                        value2=ee.Number(ee.List([value2.get(value2.keys().get(0)), -9999]).reduce(ee.Reducer.firstNonNull()));
                        im = ee.ImageCollection([im01.unmask(value).updateMask(value.neq(-9999)),im01.unmask(value2).updateMask(value2.neq(-9999))]).mean().rename('C')

                        value2 =image.select(band).updateMask(Wmask2).reduceRegion(**{'reducer':ee.Reducer.mean(),
                                                                  'geometry': region,
                                                                  'crs': crs_mod,
                                                                  'crsTransform':trans_mod,
                                                                  'bestEffort': True});
                        value2=ee.Number(ee.List([value2.get(value2.keys().get(0)), -9999]).reduce(ee.Reducer.firstNonNull()));
                        im2 = im02.unmask(value2).updateMask(value2.neq(-9999));
                        return image.addBands([im,im2]);
                    collection0 = collection.map(calc_points);
                
                    coll02=collection0.map(lambda image: image.select(['W','Q','C']).addBands(image.select(band).reduceNeighborhood(ee.Reducer.mean(), ee.Kernel.square(kern2,'pixels',True)).rename(band).reproject(crs_mod,trans_mod)));            
                    def calc_VD(image):
                        im=image.select('C').subtract(image.select(band)).divide(image.select('C').subtract(image.select('W')))
                        im=im.updateMask(im.gte(3).Or(im.lte(-2)).Not())
                        im=im.where(im.gte(1),1)
                        im=im.where(im.lte(0),0).rename(['VD'])
                        return image.addBands(im)
                        
                    coll02=coll02.map(calc_VD)
            
                    coeff=coll02.map(lambda image: image.select('W').multiply(image.select('VD')).subtract(image.select(band))).max().add(coll02.select(band).min())
                    CVMW=coll02.map(lambda image: image.select('Q').addBands(image.select('C').divide(image.select(band).subtract(image.select('W').multiply(image.select('VD'))).add(coeff))));
               
                    
                    # calculate correlation between discharge and C/M
                    corr=CVMW.reduce(ee.Reducer.spearmansCorrelation()).reproject(crs_mod,trans_mod).updateMask(Wat2);
                    
                    # find maximum correlation Calibration and create mask for that pixel
                    maxval = corr.reduceRegion(**{'reducer':ee.Reducer.max(), 'geometry':region, 'crs':crs_mod,'crsTransform':trans_mod, 'bestEffort':True});
                    Mmask_C2VW2M_k2 = corr.select('correlation').eq(ee.Number(maxval.get('correlation')));
                    Mmask_C2VW2M_k2 = Mmask_C2VW2M_k2.reduceNeighborhood(ee.Reducer.mean(), ee.Kernel.square(kern2,'pixels',True)).rename('correlation').reproject(crs_mod,trans_mod).gt(0)
                    
                    task12 = ee.batch.Export.image.toAsset(image=Mmask_C2VW2M_k2,
                                                         description='Mmask_C2V_calib_W2_VD_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr)+'_var_'+hydro+'_k'+str(kern2),
                                                         region=region,
                                                         assetId= folderGEE+'/Mmask_C2V_calib_W2_VD_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr)+'_var_'+hydro+'_k'+str(kern2),
                                                         crs=crs_mod,
                                                         crsTransform=trans_mod)
                    task12.start()
      
        # %%
        asd=1
        while asd==1:
            tasks=[]
            if no_veg==1:
                if pix_Mcal_k1_CM==1: 
                    try:
                        Mmask_C2M_k1=ee.Image(folderGEE+'/Mmask_C2_calib_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr)+'_var_'+hydro+'_k'+str(kern1)).reproject(crs_mod,trans_mod)
                        Mmask_C2M_k1.getInfo()
                    except:
                        tasks.append(task1.status()['state'])
                if pix_Mcal_k2_CM==1:
                    try:
                        Mmask_C2M_k2=ee.Image(folderGEE+'/Mmask_C2_calib_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr)+'_var_'+hydro+'_k'+str(kern2)).reproject(crs_mod,trans_mod)
                        Mmask_C2M_k2.getInfo()
                    except:
                        tasks.append(task2.status()['state'])
                if pix_Mcal_k1_CMW==1 and pix_W1==1: 
                    try:
                        Mmask_C2W1M_k1=ee.Image(folderGEE+'/Mmask_C2_calib_W1_VD_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr)+'_var_'+hydro+'_k'+str(kern1)).reproject(crs_mod,trans_mod)
                        Mmask_C2W1M_k1.getInfo()
                    except:
                        tasks.append(task3.status()['state'])
                if pix_Mcal_k1_CMW==1 and pix_W2==1: 
                    try:
                        Mmask_C2W2M_k1=ee.Image(folderGEE+'/Mmask_C2_calib_W2_VD_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr)+'_var_'+hydro+'_k'+str(kern1)).reproject(crs_mod,trans_mod)
                        Mmask_C2W2M_k1.getInfo()
                    except:
                        tasks.append(task4.status()['state'])
                if pix_Mcal_k2_CMW==1 and pix_W1==1: 
                    try:
                        Mmask_C2W1M_k2=ee.Image(folderGEE+'/Mmask_C2_calib_W1_VD_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr)+'_var_'+hydro+'_k'+str(kern2)).reproject(crs_mod,trans_mod)
                        Mmask_C2W1M_k2.getInfo()
                    except:
                        tasks.append(task5.status()['state'])
                if pix_Mcal_k2_CMW==1 and pix_W2==1: 
                    try:
                        Mmask_C2W2M_k2=ee.Image(folderGEE+'/Mmask_C2_calib_W2_VD_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr)+'_var_'+hydro+'_k'+str(kern2)).reproject(crs_mod,trans_mod)
                        Mmask_C2W2M_k2.getInfo()
                    except:
                        tasks.append(task6.status()['state'])
            if veg==1:
                if pix_Mcal_k1_CM==1: 
                    try:
                        Mmask_C2VM_k1=ee.Image(folderGEE+'/Mmask_C2V_calib_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr)+'_var_'+hydro+'_k'+str(kern1)).reproject(crs_mod,trans_mod)
                        Mmask_C2VM_k1.getInfo()
                    except:
                        tasks.append(task7.status()['state'])
                if pix_Mcal_k2_CM==1:
                    try:
                        Mmask_C2VM_k2=ee.Image(folderGEE+'/Mmask_C2V_calib_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr)+'_var_'+hydro+'_k'+str(kern2)).reproject(crs_mod,trans_mod)
                        Mmask_C2VM_k2.getInfo()
                    except:
                        tasks.append(task8.status()['state'])
                if pix_Mcal_k1_CMW==1 and pix_W1==1: 
                    try:
                        Mmask_C2VW1M_k1=ee.Image(folderGEE+'/Mmask_C2V_calib_W1_VD_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr)+'_var_'+hydro+'_k'+str(kern1)).reproject(crs_mod,trans_mod)
                        Mmask_C2VW1M_k1.getInfo()
                    except:
                        tasks.append(task9.status()['state'])
                if pix_Mcal_k1_CMW==1 and pix_W2==1: 
                    try:
                        Mmask_C2VW2M_k1=ee.Image(folderGEE+'/Mmask_C2V_calib_W2_VD_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr)+'_var_'+hydro+'_k'+str(kern1)).reproject(crs_mod,trans_mod)
                        Mmask_C2VW2M_k1.getInfo()
                    except:
                        tasks.append(task10.status()['state'])
                if pix_Mcal_k2_CMW==1 and pix_W1==1: 
                    try:
                        Mmask_C2VW1M_k2=ee.Image(folderGEE+'/Mmask_C2V_calib_W1_VD_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr)+'_var_'+hydro+'_k'+str(kern2)).reproject(crs_mod,trans_mod)
                        Mmask_C2VW1M_k2.getInfo()
                    except:
                        tasks.append(task11.status()['state'])
                if pix_Mcal_k2_CMW==1 and pix_W2==1: 
                    try:
                        Mmask_C2VW2M_k2=ee.Image(folderGEE+'/Mmask_C2V_calib_W2_VD_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr)+'_var_'+hydro+'_k'+str(kern2)).reproject(crs_mod,trans_mod)
                        Mmask_C2VW2M_k2.getInfo()
                    except:
                        tasks.append(task12.status()['state'])
                        
            if len(tasks)==0 or np.all(np.logical_or(np.unique(tasks)=='COMPLETED',np.unique(tasks)=='FAILED')):
                asd=0
                print('COMPLETE')
                if no_veg==1:
                    if pix_Mcal_k1_CM==1: 
                        try:
                            Mmask_C2M_k1=ee.Image(folderGEE+'/Mmask_C2_calib_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr)+'_var_'+hydro+'_k'+str(kern1)).reproject(crs_mod,trans_mod)
                            Mmask_C2M_k1.getInfo()
                        except:
                            Mmask_C2M_k1=im03
                    if pix_Mcal_k2_CM==1:
                        try:
                            Mmask_C2M_k2=ee.Image(folderGEE+'/Mmask_C2_calib_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr)+'_var_'+hydro+'_k'+str(kern2)).reproject(crs_mod,trans_mod)
                            Mmask_C2M_k2.getInfo()
                        except:
                            Mmask_C2M_k2=im03
                    if pix_Mcal_k1_CMW==1 and pix_W1==1: 
                        try:
                            Mmask_C2W1M_k1=ee.Image(folderGEE+'/Mmask_C2_calib_W1_VD_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr)+'_var_'+hydro+'_k'+str(kern1)).reproject(crs_mod,trans_mod)
                            Mmask_C2W1M_k1.getInfo()
                        except:
                            Mmask_C2W1M_k1=im03
                    if pix_Mcal_k1_CMW==1 and pix_W2==1: 
                        try:
                            Mmask_C2W2M_k1=ee.Image(folderGEE+'/Mmask_C2_calib_W2_VD_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr)+'_var_'+hydro+'_k'+str(kern1)).reproject(crs_mod,trans_mod)
                            Mmask_C2W2M_k1.getInfo()
                        except:
                            Mmask_C2W2M_k1=im03
                    if pix_Mcal_k2_CMW==1 and pix_W1==1: 
                        try:
                            Mmask_C2W1M_k2=ee.Image(folderGEE+'/Mmask_C2_calib_W1_VD_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr)+'_var_'+hydro+'_k'+str(kern2)).reproject(crs_mod,trans_mod)
                            Mmask_C2W1M_k2.getInfo()
                        except:
                            Mmask_C2W1M_k2=im03
                    if pix_Mcal_k2_CMW==1 and pix_W2==1: 
                        try:
                            Mmask_C2W2M_k2=ee.Image(folderGEE+'/Mmask_C2_calib_W2_VD_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr)+'_var_'+hydro+'_k'+str(kern2)).reproject(crs_mod,trans_mod)
                            Mmask_C2W2M_k2.getInfo()
                        except:
                            Mmask_C2W2M_k2=im03
                if veg==1:
                    if pix_Mcal_k1_CM==1: 
                        try:
                            Mmask_C2VM_k1=ee.Image(folderGEE+'/Mmask_C2V_calib_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr)+'_var_'+hydro+'_k'+str(kern1)).reproject(crs_mod,trans_mod)
                            Mmask_C2VM_k1.getInfo()
                        except:
                            Mmask_C2VM_k1=im03
                    if pix_Mcal_k2_CM==1:
                        try:
                            Mmask_C2VM_k2=ee.Image(folderGEE+'/Mmask_C2V_calib_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr)+'_var_'+hydro+'_k'+str(kern2)).reproject(crs_mod,trans_mod)
                            Mmask_C2VM_k2.getInfo()
                        except:
                            Mmask_C2VM_k2=im03
                    if pix_Mcal_k1_CMW==1 and pix_W1==1: 
                        try:
                            Mmask_C2VW1M_k1=ee.Image(folderGEE+'/Mmask_C2V_calib_W1_VD_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr)+'_var_'+hydro+'_k'+str(kern1)).reproject(crs_mod,trans_mod)
                            Mmask_C2VW1M_k1.getInfo()
                        except:
                            Mmask_C2VW1M_k1=im03
                    if pix_Mcal_k1_CMW==1 and pix_W2==1: 
                        try:
                            Mmask_C2VW2M_k1=ee.Image(folderGEE+'/Mmask_C2V_calib_W2_VD_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr)+'_var_'+hydro+'_k'+str(kern1)).reproject(crs_mod,trans_mod)
                            Mmask_C2VW2M_k1.getInfo()
                        except:
                            Mmask_C2VW2M_k1=im03
                    if pix_Mcal_k2_CMW==1 and pix_W1==1: 
                        try:
                            Mmask_C2VW1M_k2=ee.Image(folderGEE+'/Mmask_C2V_calib_W1_VD_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr)+'_var_'+hydro+'_k'+str(kern2)).reproject(crs_mod,trans_mod)
                            Mmask_C2VW1M_k2.getInfo()
                        except:
                            Mmask_C2VW1M_k2=im03
                    if pix_Mcal_k2_CMW==1 and pix_W2==1: 
                        try:
                            Mmask_C2VW2M_k2=ee.Image(folderGEE+'/Mmask_C2V_calib_W2_VD_'+prod_name+'_'+str2+'_'+d1+'_'+d2+'_ray_'+ray0+'_cloudthr_'+str(cld_thr)+'_var_'+hydro+'_k'+str(kern2)).reproject(crs_mod,trans_mod)
                            Mmask_C2VW2M_k2.getInfo()
                        except:
                            Mmask_C2VW2M_k2=im03
            else:
                print('saving')
                time.sleep(10)
                
        #%% Extract Data

        #preallocate selected variables
        D=[]
        C2=[]
        if veg==1:
            V=[]
        if pix_W1==1:
            W=[]
        if pix_W2==1:    
            W2=[]
        if no_veg==1:
            if pix_Mcal_k1_CM==1:
                M_cal_k1_CM=[]
            if pix_Mcal_k2_CM==1:
                M_cal_k2_CM=[]
            if pix_Mcal_k1_CMW==1 and pix_W1==1:
                M_cal_k1_CMW1=[]
            if pix_Mcal_k2_CMW==1 and pix_W1==1:              
                M_cal_k2_CMW1=[]
            if pix_Mcal_k1_CMW==1 and pix_W2==1:
                M_cal_k1_CMW2=[]
            if pix_Mcal_k2_CMW==1 and pix_W2==1:              
                M_cal_k2_CMW2=[]
        if veg==1:
            if pix_Mcal_k1_CM==1:
                M_cal_k1_CVM=[]
            if pix_Mcal_k2_CM==1:
                M_cal_k2_CVM=[]
            if pix_Mcal_k1_CMW==1 and pix_W1==1:
                M_cal_k1_CVMW1=[]
            if pix_Mcal_k2_CMW==1 and pix_W1==1:              
                M_cal_k2_CVMW1=[]
            if pix_Mcal_k1_CMW==1 and pix_W2==1:
                M_cal_k1_CVMW2=[]
            if pix_Mcal_k2_CMW==1 and pix_W2==1:              
                M_cal_k2_CVMW2=[]

        #extraction of the periods
        for k in range(len(years)-1):
            print(years[k])
            print(k)
 
            ee.Reset()
            ee.Initialize(credentials)
            band=['B8','B4'];
            print('extraction from '+str(years[k])+' to '+str(years[k+1]))
            d01=ee.Date(str(years[k])+'-01-01').millis()
            if k==len(years)-2:
                d02=ee.Date(str(years[k+1]+1)+'-01-01').millis()
            else:
                d02=ee.Date(str(years[k+1])+'-01-01').millis()
                
            coll=ee.ImageCollection(product1).filterDate(d01,d02).filterBounds(region).select(band).filter(ee.Filter.lte('CLOUDY_PIXEL_PERCENTAGE', 70));
            coll2=ee.ImageCollection(product2).filterDate(d01,d02).filterBounds(region);
            coll=coll.combine(coll2)     

            dlist=ee.List(coll.iterate(calc_Date,first)).distinct().getInfo()  
            dlist=np.sort(dlist).tolist()
            collection0=ee.ImageCollection(ee.List(dlist).map(selection))    
            collection0=collection0.filterMetadata('Valid_perc',"not_less_than",0.20)
            dlistn=ee.List(collection0.iterate(calc_Date,first)).getInfo()
            
            print('after cloud masking length: '+str(len(dlistn)))

            def mosaic2(nday):
                image=coll.filterDate(ee.Date(nday),ee.Date(ee.Number(nday).add(86400000))).reduce(ee.Reducer.mean()).rename(band2)
                image=image.select(band).updateMask(image.select('probability').lte(cld_thr)).multiply(mult).set('system:time_start',nday)
                image=image.updateMask(image.gt(0))
                
                return image.set('system:time_start',nday)
                
            # EXTRACTING COLLECTION WITH VALID IMAGES
            collection=ee.ImageCollection(ee.List(dlistn).map(mosaic2))
            band=[band[0]]

            #extraction satellite date
            def calc_Date_oss(image,llist):
                value =image.get('system:time_start')
                value=ee.Date(ee.List([value, -9999]).reduce(ee.Reducer.firstNonNull())).millis();
                return ee.List(llist).add(value) 
            asd=1
            while asd==1:
                print('calc_D')
                try:
                    D0=collection.iterate(calc_Date_oss,first).getInfo()
                    asd=0
                except:
                    print('failure')
            D=np.append(D,np.array(D0).astype('float'))

            #extraction C2
            def calc_C2(image,llist):
                value =image.select(band).updateMask(C2mask).reduceRegion(**{'reducer':ee.Reducer.mean(),
                                                          'geometry':region,
                                                          'crs':crs_mod,
                                                          'crsTransform':trans_mod,
                                                          'bestEffort':True}).get(band[0]);
                value=ee.Number(ee.List([value, -9999]).reduce(ee.Reducer.firstNonNull()));
                return ee.List(llist).add(value)
            asd=1
            while asd==1:
                print('calc_C2')
                try:
                    C20=collection.iterate(calc_C2,first).getInfo()
                    asd=0
                except:
                    print('failure')
            C2=np.append(C2,np.array(C20).astype('float'))
        
            #extraction V
            if veg==1:
                def calc_V(image,llist):
                    value =image.select(band).updateMask(Vmask).reduceRegion(**{'reducer':ee.Reducer.mean(),
                                                              'geometry':region,
                                                              'crs':crs_mod,
                                                              'crsTransform':trans_mod,
                                                              'bestEffort':True}).get(band[0]);
                    value=ee.Number(ee.List([value, -9999]).reduce(ee.Reducer.firstNonNull()));
                    return ee.List(llist).add(value)
                asd=1
                while asd==1:
                    print('calc_V')
                    try:
                        V0=collection.iterate(calc_V,first).getInfo()
                        asd=0
                    except:
                        print('failure')                
                V=np.append(V,np.array(V0).astype('float'))
            
            #extraction W1 
            if pix_W1==1:
                def calc_W(image,llist):
                    value =image.select(band).updateMask(Wmask).reduceRegion(**{'reducer':ee.Reducer.mean(),
                                                              'geometry':region,
                                                              'crs':crs_mod,
                                                              'crsTransform':trans_mod,
                                                              'bestEffort':True}).get(band[0]);
                    value=ee.Number(ee.List([value, -9999]).reduce(ee.Reducer.firstNonNull()));
                    return ee.List(llist).add(value)
                asd=1
                while asd==1:
                     print('calc_W')
                     try:
                         W0=collection.iterate(calc_W,first).getInfo()
                         asd=0
                     except:
                         print('failure')
                W=np.append(W,np.array(W0).astype('float'))

            #extraction W2 
            if pix_W2==1:
                def calc_W2(image,llist):
                    value =image.select(band).updateMask(Wmask2).reduceRegion(**{'reducer':ee.Reducer.mean(),
                                                              'geometry':region,
                                                              'crs':crs_mod,
                                                              'crsTransform':trans_mod,
                                                              'bestEffort':True}).get(band[0]);
                    value=ee.Number(ee.List([value, -9999]).reduce(ee.Reducer.firstNonNull()));
                    return ee.List(llist).add(value)
                asd=1
                while asd==1:
                    print('calc_W2')
                    try:
                        W02=collection.iterate(calc_W2,first).getInfo()
                        asd=0
                    except:
                        print('failure')
                W2=np.append(W2,np.array(W02).astype('float'))

            if no_veg==1:                
                #extraction M cal CM kernel 1
                if pix_Mcal_k1_CM==1:
                    def calc_M_cal_k1_CM(image,llist):
                        value =image.select(band).updateMask(Mmask_C2M_k1).reduceRegion(**{'reducer':ee.Reducer.mean(),
                                                                  'geometry':region,
                                                                  'crs':crs_mod,
                                                                  'crsTransform':trans_mod,
                                                                  'bestEffort':True}).get(band[0]);
                        value=ee.Number(ee.List([value, -9999]).reduce(ee.Reducer.firstNonNull()));
                        return ee.List(llist).add(value)
                    asd=1
                    while asd==1:
                        print('calc_M_cal_k1_CM')
                        try:
                            M_cal_k1_CM0=collection.iterate(calc_M_cal_k1_CM,first).getInfo()
                            asd=0
                        except:
                            print('failure')
                    M_cal_k1_CM=np.append(M_cal_k1_CM,np.array(M_cal_k1_CM0).astype('float'))        

                #extraction M cal CM kernel 2
                if pix_Mcal_k2_CM==1:
                    def calc_M_cal_k2_CM(image,llist):
                        value =image.select(band).updateMask(Mmask_C2M_k2).reduceRegion(**{'reducer':ee.Reducer.mean(),
                                                                  'geometry':region,
                                                                  'crs':crs_mod,
                                                                  'crsTransform':trans_mod,
                                                                  'bestEffort':True}).get(band[0]);
                        value=ee.Number(ee.List([value, -9999]).reduce(ee.Reducer.firstNonNull()));
                        return ee.List(llist).add(value)
                    asd=1
                    while asd==1:
                        print('calc_M_cal_k2_CM')
                        try:
                            M_cal_k2_CM0=collection.iterate(calc_M_cal_k2_CM,first).getInfo()
                            asd=0
                        except:
                            print('failure')
                    M_cal_k2_CM=np.append(M_cal_k2_CM,np.array(M_cal_k2_CM0).astype('float'))        

                #extraction M cal CMW1 kernel 1
                if pix_Mcal_k1_CMW==1 and pix_W1==1:
                    def calc_M_cal_k1_CMW1(image,llist):
                        value =image.select(band).updateMask(Mmask_C2W1M_k1).reduceRegion(**{'reducer':ee.Reducer.mean(),
                                                                  'geometry':region,
                                                                  'crs':crs_mod,
                                                                  'crsTransform':trans_mod,
                                                                  'bestEffort':True}).get(band[0]);
                        value=ee.Number(ee.List([value, -9999]).reduce(ee.Reducer.firstNonNull()));
                        return ee.List(llist).add(value)
                    asd=1
                    while asd==1:
                        print('calc_M_cal_k1_CMW1')
                        try:
                            M_cal_k1_CMW10=collection.iterate(calc_M_cal_k1_CMW1,first).getInfo()
                            asd=0
                        except:
                            print('failure')
                    M_cal_k1_CMW1=np.append(M_cal_k1_CMW1,np.array(M_cal_k1_CMW10).astype('float'))        

                #extraction M cal CMW1 kernel 2
                if pix_Mcal_k2_CMW==1 and pix_W1==1:              
                    def calc_M_cal_k2_CMW1(image,llist):
                        value =image.select(band).updateMask(Mmask_C2W1M_k2).reduceRegion(**{'reducer':ee.Reducer.mean(),
                                                                  'geometry':region,
                                                                  'crs':crs_mod,
                                                                  'crsTransform':trans_mod,
                                                                  'bestEffort':True}).get(band[0]);
                        value=ee.Number(ee.List([value, -9999]).reduce(ee.Reducer.firstNonNull()));
                        return ee.List(llist).add(value)
                    asd=1
                    while asd==1:
                        print('calc_M_cal_k2_CMW1')
                        try:
                            M_cal_k2_CMW10=collection.iterate(calc_M_cal_k2_CMW1,first).getInfo()
                            asd=0
                        except:
                            print('failure')
                    M_cal_k2_CMW1=np.append(M_cal_k2_CMW1,np.array(M_cal_k2_CMW10).astype('float'))
                    
                #extraction M cal CMW2 kernel 1
                if pix_Mcal_k1_CMW==1 and pix_W2==1:
                    def calc_M_cal_k1_CMW2(image,llist):
                        value =image.select(band).updateMask(Mmask_C2W2M_k1).reduceRegion(**{'reducer':ee.Reducer.mean(),
                                                                  'geometry':region,
                                                                  'crs':crs_mod,
                                                                  'crsTransform':trans_mod,
                                                                  'bestEffort':True}).get(band[0]);
                        value=ee.Number(ee.List([value, -9999]).reduce(ee.Reducer.firstNonNull()));
                        return ee.List(llist).add(value)
                    asd=1
                    while asd==1:
                        print('calc_M_cal_k1_CMW2')
                        try:
                            M_cal_k1_CMW20=collection.iterate(calc_M_cal_k1_CMW2,first).getInfo()
                            asd=0
                        except:
                            print('failure')
                    M_cal_k1_CMW2=np.append(M_cal_k1_CMW2,np.array(M_cal_k1_CMW20).astype('float'))        

                #extraction M cal CMW2 kernel 2
                if pix_Mcal_k2_CMW==1 and pix_W2==1:              
                    def calc_M_cal_k2_CMW2(image,llist):
                        value =image.select(band).updateMask(Mmask_C2W2M_k2).reduceRegion(**{'reducer':ee.Reducer.mean(),
                                                                  'geometry':region,
                                                                  'crs':crs_mod,
                                                                  'crsTransform':trans_mod,
                                                                  'bestEffort':True}).get(band[0]);
                        value=ee.Number(ee.List([value, -9999]).reduce(ee.Reducer.firstNonNull()));
                        return ee.List(llist).add(value)
                    asd=1
                    while asd==1:
                        print('calc_M_cal_k2_CMW2')
                        try:
                            M_cal_k2_CMW20=collection.iterate(calc_M_cal_k2_CMW2,first).getInfo()
                            asd=0
                        except:
                            print('failure')
                    M_cal_k2_CMW2=np.append(M_cal_k2_CMW2,np.array(M_cal_k2_CMW20).astype('float'))


            if veg==1:                
                #extraction M cal CM kernel 1
                if pix_Mcal_k1_CM==1:
                    def calc_M_cal_k1_CVM(image,llist):
                        value =image.select(band).updateMask(Mmask_C2VM_k1).reduceRegion(**{'reducer':ee.Reducer.mean(),
                                                                  'geometry':region,
                                                                  'crs':crs_mod,
                                                                  'crsTransform':trans_mod,
                                                                  'bestEffort':True}).get(band[0]);
                        value=ee.Number(ee.List([value, -9999]).reduce(ee.Reducer.firstNonNull()));
                        return ee.List(llist).add(value)
                    asd=1
                    while asd==1:
                        print('calc_M_cal_k1_CVM')
                        try:
                            M_cal_k1_CVM0=collection.iterate(calc_M_cal_k1_CVM,first).getInfo()
                            asd=0
                        except:
                            print('failure')
                    M_cal_k1_CVM=np.append(M_cal_k1_CVM,np.array(M_cal_k1_CVM0).astype('float'))        

                #extraction M cal CM kernel 2
                if pix_Mcal_k2_CM==1:
                    def calc_M_cal_k2_CVM(image,llist):
                        value =image.select(band).updateMask(Mmask_C2VM_k2).reduceRegion(**{'reducer':ee.Reducer.mean(),
                                                                  'geometry':region,
                                                                  'crs':crs_mod,
                                                                  'crsTransform':trans_mod,
                                                                  'bestEffort':True}).get(band[0]);
                        value=ee.Number(ee.List([value, -9999]).reduce(ee.Reducer.firstNonNull()));
                        return ee.List(llist).add(value)
                    asd=1
                    while asd==1:
                        print('calc_M_cal_k2_CVM')
                        try:
                            M_cal_k2_CVM0=collection.iterate(calc_M_cal_k2_CVM,first).getInfo()
                            asd=0
                        except:
                            print('failure')
                    M_cal_k2_CVM=np.append(M_cal_k2_CVM,np.array(M_cal_k2_CVM0).astype('float'))        

                #extraction M cal CMW1 kernel 1
                if pix_Mcal_k1_CMW==1 and pix_W1==1:
                    def calc_M_cal_k1_CVMW1(image,llist):
                        value =image.select(band).updateMask(Mmask_C2VW1M_k1).reduceRegion(**{'reducer':ee.Reducer.mean(),
                                                                  'geometry':region,
                                                                  'crs':crs_mod,
                                                                  'crsTransform':trans_mod,
                                                                  'bestEffort':True}).get(band[0]);
                        value=ee.Number(ee.List([value, -9999]).reduce(ee.Reducer.firstNonNull()));
                        return ee.List(llist).add(value)
                    asd=1
                    while asd==1:
                        print('calc_M_cal_k1_CVMW1')
                        try:
                            M_cal_k1_CVMW10=collection.iterate(calc_M_cal_k1_CVMW1,first).getInfo()
                            asd=0
                        except:
                            print('failure')
                    M_cal_k1_CVMW1=np.append(M_cal_k1_CVMW1,np.array(M_cal_k1_CVMW10).astype('float'))        

                #extraction M cal CMW1 kernel 2
                if pix_Mcal_k2_CMW==1 and pix_W1==1:              
                    def calc_M_cal_k2_CVMW1(image,llist):
                        value =image.select(band).updateMask(Mmask_C2VW1M_k2).reduceRegion(**{'reducer':ee.Reducer.mean(),
                                                                  'geometry':region,
                                                                  'crs':crs_mod,
                                                                  'crsTransform':trans_mod,
                                                                  'bestEffort':True}).get(band[0]);
                        value=ee.Number(ee.List([value, -9999]).reduce(ee.Reducer.firstNonNull()));
                        return ee.List(llist).add(value)
                    asd=1
                    while asd==1:
                        print('calc_M_cal_k2_CVMW1')
                        try:
                            M_cal_k2_CVMW10=collection.iterate(calc_M_cal_k2_CVMW1,first).getInfo()
                            asd=0
                        except:
                            print('failure')
                    M_cal_k2_CVMW1=np.append(M_cal_k2_CVMW1,np.array(M_cal_k2_CVMW10).astype('float'))
                    
                #extraction M cal CMW2 kernel 1
                if pix_Mcal_k1_CMW==1 and pix_W2==1:
                    def calc_M_cal_k1_CVMW2(image,llist):
                        value =image.select(band).updateMask(Mmask_C2VW2M_k1).reduceRegion(**{'reducer':ee.Reducer.mean(),
                                                                  'geometry':region,
                                                                  'crs':crs_mod,
                                                                  'crsTransform':trans_mod,
                                                                  'bestEffort':True}).get(band[0]);
                        value=ee.Number(ee.List([value, -9999]).reduce(ee.Reducer.firstNonNull()));
                        return ee.List(llist).add(value)
                    asd=1
                    while asd==1:
                        print('calc_M_cal_k1_CVMW2')
                        try:
                            M_cal_k1_CVMW20=collection.iterate(calc_M_cal_k1_CVMW2,first).getInfo()
                            asd=0
                        except:
                            print('failure')
                    M_cal_k1_CVMW2=np.append(M_cal_k1_CVMW2,np.array(M_cal_k1_CVMW20).astype('float'))        

                #extraction M cal CMW2 kernel 2
                if pix_Mcal_k2_CMW==1 and pix_W2==1:              
                    def calc_M_cal_k2_CVMW2(image,llist):
                        value =image.select(band).updateMask(Mmask_C2VW2M_k2).reduceRegion(**{'reducer':ee.Reducer.mean(),
                                                                  'geometry':region,
                                                                  'crs':crs_mod,
                                                                  'crsTransform':trans_mod,
                                                                  'bestEffort':True}).get(band[0]);
                        value=ee.Number(ee.List([value, -9999]).reduce(ee.Reducer.firstNonNull()));
                        return ee.List(llist).add(value)
                    asd=1
                    while asd==1:
                        print('calc_M_cal_k2_CVMW2')
                        try:
                            M_cal_k2_CVMW20=collection.iterate(calc_M_cal_k2_CVMW2,first).getInfo()
                            asd=0
                        except:
                            print('failure')
                    M_cal_k2_CVMW2=np.append(M_cal_k2_CVMW2,np.array(M_cal_k2_CVMW20).astype('float'))



        D[D==-9999]=np.nan
        C2[C2==-9999]=np.nan
        if veg==1:
            V[V==-9999]=np.nan
        if pix_W1==1:
            W[W==-9999]=np.nan
        if pix_W2==1:
            W2[W2==-9999]=np.nan
        if no_veg==1:               
            if pix_Mcal_k1_CM==1:
                M_cal_k1_CM[M_cal_k1_CM==-9999]=np.nan
            if pix_Mcal_k2_CM==1:
                M_cal_k2_CM[M_cal_k2_CM==-9999]=np.nan
            if pix_Mcal_k1_CMW==1 and pix_W1==1:
                M_cal_k1_CMW1[M_cal_k1_CMW1==-9999]=np.nan
            if pix_Mcal_k2_CMW==1 and pix_W1==1:              
                M_cal_k2_CMW1[M_cal_k2_CMW1==-9999]=np.nan
            if pix_Mcal_k1_CMW==1 and pix_W2==1:
                M_cal_k1_CMW2[M_cal_k1_CMW2==-9999]=np.nan
            if pix_Mcal_k2_CMW==1 and pix_W2==1:              
                M_cal_k2_CMW2[M_cal_k2_CMW2==-9999]=np.nan
        if veg==1:               
            if pix_Mcal_k1_CM==1:
                M_cal_k1_CVM[M_cal_k1_CVM==-9999]=np.nan
            if pix_Mcal_k2_CM==1:
                M_cal_k2_CVM[M_cal_k2_CVM==-9999]=np.nan
            if pix_Mcal_k1_CMW==1 and pix_W1==1:
                M_cal_k1_CVMW1[M_cal_k1_CVMW1==-9999]=np.nan
            if pix_Mcal_k2_CMW==1 and pix_W1==1:              
                M_cal_k2_CVMW1[M_cal_k2_CVMW1==-9999]=np.nan
            if pix_Mcal_k1_CMW==1 and pix_W2==1:
                M_cal_k1_CVMW2[M_cal_k1_CVMW2==-9999]=np.nan
            if pix_Mcal_k2_CMW==1 and pix_W2==1:              
                M_cal_k2_CVMW2[M_cal_k2_CVMW2==-9999]=np.nan

        D=D/1000/3600/24
        D=np.array ([datetime(1970,1,1)+timedelta(days=d) for d in D])
        Dmat=np.array([d.toordinal()+366 for d in D])
        dic={'D':Dmat,'C2':C2,'IDclim':IDclim}
        if veg==1:
            dic['V']=V
        if pix_W1==1:
            dic['W']=W
        if pix_W2==1:
            dic['W2']=W2
        if no_veg==1:
            if pix_Mcal_k1_CM==1:
                dic['M_cal_k1_CM']=M_cal_k1_CM
            if pix_Mcal_k2_CM==1:
                dic['M_cal_k2_CM']=M_cal_k1_CMW1
            if pix_Mcal_k1_CMW==1 and pix_W1==1:
                dic['M_cal_k1_CMW1']=M_cal_k1_CMW1
            if pix_Mcal_k2_CMW==1 and pix_W1==1:              
                dic['M_cal_k2_CMW1']=M_cal_k2_CMW1
            if pix_Mcal_k1_CMW==1 and pix_W2==1:
                dic['M_cal_k1_CMW2']=M_cal_k1_CMW2
            if pix_Mcal_k2_CMW==1 and pix_W2==1:              
                dic['M_cal_k2_CMW2']=M_cal_k2_CMW2
        if veg==1:
            if pix_Mcal_k1_CM==1:
                dic['M_cal_k1_CVM']=M_cal_k1_CVM
            if pix_Mcal_k2_CM==1:
                dic['M_cal_k2_CVM']=M_cal_k2_CVM
            if pix_Mcal_k1_CMW==1 and pix_W1==1:
                dic['M_cal_k1_CVMW1']=M_cal_k1_CVMW1
            if pix_Mcal_k2_CMW==1 and pix_W1==1:              
                dic['M_cal_k2_CVMW1']=M_cal_k2_CVMW1
            if pix_Mcal_k1_CMW==1 and pix_W2==1:
                dic['M_cal_k1_CVMW2']=M_cal_k1_CVMW2
            if pix_Mcal_k2_CMW==1 and pix_W2==1:              
                dic['M_cal_k2_CVMW2']=M_cal_k2_CVMW2

        savemat(folder_pc+'/data_full_'+prod_name+'_'+str2+'_calval2.mat',dic)
        
        attempt=0
    except:
        print('failure analysis')
    