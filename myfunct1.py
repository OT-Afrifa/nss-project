import numpy as np
import xarray as xr
from scipy import stats as sc
import pandas as pd
import os

#calculate Rainfall Anomaly Index 
def rai(ds, dimension, method='ordinary'):
    
    ds_ma = ds.mean(dimension)
    
    ds_anom = ds - ds_ma
    
    if method.lower() == 'percentile':
        l_thresh = ds.reduce(np.nanpercentile,q=10,dim=dimension)
        u_thresh = ds.reduce(np.nanpercentile,q=90,dim=dimension)
        ds_low_10 = ds.where(ds<l_thresh).mean(dimension)
        ds_high_10 = ds.where(ds>u_thresh).mean(dimension)
    
    elif method.lower() == 'ordinary':
        thresh = ds.reduce(np.sort,dim=dimension)
        ds_low_10 = thresh[:10].mean(dimension)
        ds_high_10 = thresh[:-10:-1].mean(dimension)
    
    else:
        print('Wrong/No method selected.')
    
    negatives = -3*( (ds_anom.where(ds_anom<0)) / (ds_low_10-ds_ma) )
    positives = 3*( (ds_anom.where(ds_anom>0)) / (ds_high_10-ds_ma) )
    RAI = ds_anom.where(ds_anom>=0, negatives).where(ds_anom<0, positives)
    
    return RAI



#calculate Standardized Precipitation Index (Function for 3D)
def spi(ds, thresh, dimension):
    #ds - data ; thresh - time interval / scale; dimension - dimension as a string
    
    #Rolling Mean / Moving Averages
    ds_ma = ds.rolling(time = thresh, center=False).mean(dim=dimension)
    
    #Natural log of moving averages
    ds_In = np.log(ds_ma)
    ds_In = ds_In.where(np.isinf(ds_In) == False) #= np.nan  #Change infinity to NaN
    
    #Overall Mean of Moving Averages
    ds_mu = ds_ma.mean(dimension)
    
    #Summation of Natural log of moving averages
    ds_sum = ds_In.sum(dimension)
    
    #Computing essentials for gamma distribution
    n = ds_In[thresh-1:, :, :].count(dimension)                  #size of data
    A = np.log(ds_mu) - (ds_sum/n)             #Computing A
    alpha = (1/(4*A))*(1+(1+((4*A)/3))**0.5)   #Computing alpha  (a)
    beta = ds_mu/alpha                         #Computing beta (scale)
    
    #Gamma Distribution (CDF) 
    gamma_func = lambda data, a, scale: sc.gamma.cdf(data, a=a, scale=scale)
    gamma = xr.apply_ufunc(gamma_func, ds_ma, alpha, beta)
    
    #Standardized Precipitation Index   (Inverse of CDF)
    norminv = lambda data: sc.norm.ppf(data, loc=0, scale=1)
    norm_spi = xr.apply_ufunc(norminv, gamma)  #loc is mean and scale is standard dev.
    
    return ds_ma, ds_In , ds_mu, ds_sum,n, A, alpha, beta, gamma, norm_spi


#Calculate the Standardized Anomaly Index
def sai(ds, dimension):
    return (ds-ds.mean(dimension))/ds.std(dimension)
    




#unzip function
def unzip(z_file, model, print_info=True):
    with zipfile.ZipFile(z_file, 'r') as zip_ref:
        zip_ref.extractall(model)
        if print_info==True:
            print('Done with '+model+' model')
            



#filter values between x and y and count values
def filtercounter(ds, x, y):
    return ds.where((ds > x) & (ds <= y)).count().values.tolist()


#filter values between x and y and find mean of values
def filtermean(ds, x, y):
    return ds.where((ds > x) & (ds <= y)).reduce(np.nanmean).values.tolist()

                
#create background for polar plot
def background(ax,x):
    dx = 20
    theta = np.linspace(0,np.pi,100)#*100
    r = np.sin(theta)*(100-dx)

    #ax = plt.subplot(2,2,(o+1), polar=True)
    #ax = x.plot(polar=True)
    c = ax.scatter(theta, r, c=r, s=20, cmap='Blues', alpha=0.0)

    ax.set_thetamin(0)
    ax.set_thetamax(180)

    ax.set_yticks(np.arange(0,60.1,10),)#np.arange(-50,100.1,10))
    ax.text( -0.35, 40, 'Events (%)', color='k', fontsize = 15)
    ax.set_xticks(np.arange(0,np.pi+0.1,np.pi/8))

    perc = ['< -3', '-3', '-2', '-1', '', '1', '2', '3', '> 3']
    ax.set_xticklabels(perc[::-1])
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')

    #Descriptors
    ax.text( 1.65, 116-dx, 'Normal', color='r', fontsize = 15)
    ax.text( 2.2, 117-dx, 'Dry', color='r', fontsize = 14, rotation=35)
    ax.text( 2.65, 117-dx, 'Very Dry', color='r', fontsize = 14, rotation=60)
    ax.text( 3.075, 120-dx, "Extremely\n Dry", color='r', fontsize = 14, rotation=82.5)

    ax.text( 1, 114-dx, 'Wet', color='r', fontsize = 14, rotation=-35)
    ax.text( 0.57, 108-dx, 'Very Wet', color='r', fontsize = 14, rotation=-60)
    ax.text( 0.03, 108-dx, "Extremely\n Wet", color='r', fontsize = 14, rotation=-82.5)


    ax.tick_params(axis='both', labelsize=12)
    ax.spines['inner'].set_color('r')
    ax.spines['polar'].set_linewidth(2)
    #ax.spines['end'].set_linewidth(0)

    ax.grid(True, which='major', axis='both', linestyle='--', color = 'b', linewidth=1)
    ax.get_xgridlines()[4].set_linestyle('None')
    ax.get_xgridlines()[-1].set_linestyle('None')
    ax.set_title(x, fontsize = 15)
    colors=['r']

    return ax



#Plot pointers on polar plot background
def pointers(ax, num, percentage, linestyle = None, color = None, marker=None, marksize=None, label=None, show_legend=None):
    b=num
    pos=np.int_(b)+(np.ones(len(b))*4)
    #pos = int(b)+4  #position of the whole number component of the number within the perc array.
    perc = ['< -3', '-3', '-2', '-1', '', '1', '2', '3', '> 3']
    whole = (((pos)/(len(perc)-1))*np.pi)       #change whole part to angle
    bot = np.ones(len(b))*(len(perc)-1)
    bot1 = np.ones(len(b))*np.pi
    frac = ((b-np.int_(b))/bot)*bot1    #change fractional part to angle

    if color == None:
        ax.plot(np.pi-whole-frac,percentage,linestyle=linestyle, linewidth=2, marker = marker, markersize=marksize, label=label)
    else:
        ax.plot(np.pi-whole-frac,percentage,linestyle=linestyle, color=color, linewidth=5, marker = marker, markersize=marksize, label=label)
   
    #ax.legend(loc=0, ncol=3, bbox_to_anchor=(1.25, 1))

    if show_legend==True:
        ax.legend(loc=0, ncol=2, bbox_to_anchor=(2.6, 0.8), fontsize=14)
    #else:
    #    ax.legend(loc=0, ncol=10, bbox_to_anchor=(1.25, 1), alpha=0.0)
    

    
    
def cumulative_mean(x, dim):
    output = xr.DataArray()
    output = xr.concat([x.cumsum(dim=dim)[i-1] for i in np.arange(1,x[dim].size+1)], dim)
    return output
    

    
def trends_eval(data, case='drought', threshold=1, ttype='spatial', dim='time', lons=None, lats=None, start_year=None, stop_year=None, set_historical=False, historical_year=None):
    '''
    data = Data
    case = drought / flood
    threshold = thrshold for drought (-1) or flood (1)
    ttype = can be either spatial (XD) or point (1D)
    dim = dimension along which to compute the assessment metrics.
    lons = longitude values
    lats = latitude values
    start_year = the starting year for metric computation
    stop_year = the end year for metric computation
    '''
    
    
    import numpy as np
    
    error = 'Longitude or latitude not a point. Check and revert'
    if ttype=='spatial':
        if len(np.shape(lons))>1:
            success_type = 1
            left = lons[0]; right = lons[1]
            bottom = lats[0]; top = lats[1]
            try:
                data = data.sel(longitude=slice(left,right), latitude=slice(bottom,top))
            except KeyError:
                data = data.sel(lon=slice(left,right), lat=slice(bottom,top))
        else:
            success_type = 1
            left = lons; bottom =lats
            try:
                data = data
                #data = data.sel(longitude=left,latitude=bottom, method='nearest')
            except:
                #KeyError:
                data = data.sel(lon=left,lat=bottom, method='nearest')
            
    elif ttype=='point':
        if len(np.shape(lons))>1:
            error_type = 1
        else:
            success_type = 1
            left = lons; bottom =lats
            try:
                data = data
                #data = data.sel(lon=left,lat=bottom, method='nearest')
                #data = data.sel(longitude=left,latitude=bottom, method='nearest')
            except KeyError:
                #data = data.sel(lon=left,lat=bottom, method='nearest')
                data = data.sel(longitude=left,latitude=bottom, method='nearest')        
        
    ####Check Time Limits
    if start_year == None:
        try:
            start_year = np.nanmin(data['time.year'])
        except:
            start_year = np.nanmin(data['year'])
    
    if stop_year == None:
        try:
            stop_year = np.nanmax(data['time.year'])
        except:
            stop_year = np.nanmax(data['year'])
         
        
    ### Check Historical Criteria  
    if set_historical==True:
        if historical_year != None:
            try:
                historical_year = int(historical_year)
            except:
                output='Check the data format of historical year. Should be an integer and not None.'
                pass
            break_year = historical_year    #### Split data into historical period and projection period.
            historical_period = [start_year,break_year]; projection_period = [break_year, stop_year]
            
            if start_year>break_year:
                try:
                    int(None)
                except:
                    output = 'Start_year can not be later than the historical year. Check and revert.'
                    pass

            if stop_year<break_year:
                try:
                    int(None)
                except:
                    output = 'Stop_year can not be earlier than the historical year. Check and revert.'
                    pass
  
        else:
            print('Historical year can not be None when set_historical is True.')
            try:
                int(historical_year)
            except:
                output = 'Historical year can not be None when set_historical is True.'
                pass
            
        
      
    
  

    
    if success_type == 1:
        if set_historical == True:
            tmp = data.sel(time=slice(str(historical_period[0]),str(historical_period[1])))
            if case == 'drought':
                output = tmp.where(tmp<=threshold).count(dim) / tmp.count(dim)
            elif case == 'flood':
                output = tmp.where(tmp>=threshold).count(dim) / tmp.count(dim)       
            start_year = break_year
            
        else:
            output = xr.DataArray()
            
        for i,year in enumerate(range(start_year, stop_year+1)):
            #if dim == 'time':
            tmp = data.sel(time=str(year))
            if case == 'drought':
                #output.append(tmp.where(tmp<=threshold).count('time') / tmp.count('time'))
                if i ==0:
                    output = tmp.where(tmp<=threshold).count('time') / tmp.count('time')
                else:
                    output = xr.concat([output, tmp.where(tmp<=threshold).count('time') / tmp.count('time')], 'time')
            elif case == 'flood':
                #output.append(tmp.where(tmp>=threshold).count('time') / tmp.count('time'))
                if i ==0:
                    output = tmp.where(tmp>=threshold).count('time') / tmp.count('time')
                else:
                    output = xr.concat([output, tmp.where(tmp>=threshold).count('time') / tmp.count('time')], 'time')

                
    else:
        output = error
        
    
    #### Generating the Cumulative Means
    try:
        output = (output.cumsum(dim=dim)/range(1,output[dim].size+1))
    except:
        #output = cumulative_mean(output, dim)
        pass
    
    
    return output
   


def evaluate_dght_rai(c1):
    ######## Drought Assessment
    eval_dght_c1 = trends_eval(c1, case='drought', threshold=-1, ttype='point', dim='time', 
                      lons=None, lats=None, start_year=None, stop_year=None, 
                      set_historical=False, historical_year=None)

    return eval_dght_c1


def evaluate_dght_sai(c2):
    ######## Drought Assessment sai
    eval_dght_c2 = trends_eval(c2, case='drought', threshold=-1, ttype='point', dim='time', 
                      lons=None, lats=None, start_year=None, stop_year=None, 
                      set_historical=False, historical_year=None)
    
    return eval_dght_c2

    
    

def evaluate_fld_sai(c2):
    eval_fld_c2 = trends_eval(c2, case='flood', threshold=1, ttype='point', dim='time', 
                  lons=None, lats=None, start_year=None, stop_year=None, 
                  set_historical=False, historical_year=None)
    return eval_fld_c2


def evaluate_fld_rai(c1):
    eval_fld_c1 = trends_eval(c1, case='flood', threshold=1, ttype='point', dim='time', 
                  lons=None, lats=None, start_year=None, stop_year=None, 
                  set_historical=False, historical_year=None)
    return eval_fld_c1


def dif(ds, p1):
    return (ds.sel(year=p1).mean()-ds.sel(year=slice('1950','2020')).mean()).round(2)


def save_plot(path, filename):
    import matplotlib.pyplot as plt
    if os.path.isdir(path) == False:
        os.mkdir(path)
        
    plt.savefig(path+'/'+filename+'.jpg')