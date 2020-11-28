import numpy as np
import pandas as pd
import sys
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib
matplotlib.rc('text',usetex=True)

def readin_table(infile,header_list,header_passed=False):
    #This will read the table into a pandas dataframe
    #We need to read the table in from exoplanets archive with some parameters
    if header_passed==False:
        dt = pd.read_csv(infile,sep=',',comment='#',header=0,index_col='loc_rowid')
    else:
        dt = pd.read_csv(infile,sep=',',comment='#',header=None,names=header_list,index_col=False)
    return dt

def scat_plotter(table_in,xcol,ycol,xlab,ylab,logscaleX=True,logscaleY=True,svfig=False,cdense=False):
    matplotlib.rc('font',size=15)
    fig,ax = plt.subplots(figsize=(10,8))
    ax.scatter(table_in[xcol].values[:],table_in[ycol].values[:],marker='.',color='black')
    if cdense==True:
        x_limits = ax.get_xlim()
        y_limits = ax.get_ylim()
        wgrad = np.where(np.isfinite(table_in[xcol].values[:]))[0]
        min_rad, max_rad = np.amin(table_in[xcol].values[:][wgrad]),np.amax(table_in[xcol].values[:][wgrad])
        mj = 317.8*5.9736e27 #in gm
        rj = 11.209*6.378136e8 #in cm
        vj = (4/3)*np.pi*rj**(3.0)
        rho_factor = vj / mj #converts densities into jupiter mass/radius units
        rcd = np.linspace(min_rad,max_rad,1000)
        vcd = (4/3)*np.pi*rcd**(3.0)
        rho_arr = np.array([0.5,1,2,3,5,7])
        for i in range(6):
            mcd = rho_arr[i] * vcd
            rho_str = r'$\rho=$%.1f g cm$^{-3}$'%rho_arr[i]
            ax.plot(rcd,mcd,label=rho_str)
        ax.set_xlim(x_limits)
        ax.set_ylim(y_limits)
    ax.set_xlabel(xlab,fontsize=15)
    ax.set_ylabel(ylab,fontsize=15)
    ax.tick_params(axis='both',direction='in')
    ax.tick_params(axis='both',which='minor',direction='in')
    ax.tick_params(top=True,right=True)
    ax.tick_params(which='minor',top=True,right=True)
    if logscaleX==True:
        ax.set_xscale('log',nonposx='clip')
    if logscaleY==True:
        ax.set_yscale('log',nonposy='clip')
    if cdense==True:
        ax.legend()
    if svfig==True:
        fgname = '{}_v_{}.png'.format(ycol,xcol)
        fig.savefig(fgname,format='png')
    else:
        plt.show()

def histo_plotter(table_in,cname,xlab,num_bins,logscaleX=True,logscaleY=True,svfig=False):
    matplotlib.rc('font',size=15)
    figh,axh = plt.subplots(figsize=(10,8))
    temp_arr = np.array(table_in[cname].values[:])
    #wnan = np.where(np.isnan(temp_arr))[0]
    #temp_arr[wnan] = 0
    wfin = np.where(np.isfinite(temp_arr))[0]
    temp_arr2 = np.sort(temp_arr[wfin])
    min_val = np.amin(temp_arr2)
    max_val = np.amax(temp_arr2)

    if logscaleX==True:
        axh.hist(temp_arr2,bins=10**np.linspace(np.log10(min_val),np.log10(max_val),num_bins),histtype='step',color='black')
    else:
        axh.hist(temp_arr2,bins=num_bins,histtype='step',color='black')
    axh.set_xlabel(xlab,fontsize=15)
    axh.set_ylabel('Counts',fontsize=15)
    axh.tick_params(axis='both',direction='in')
    axh.tick_params(axis='both',which='minor',direction='in')
    axh.tick_params(top=True,right=True)
    axh.tick_params(which='minor',top=True,right=True)
    if logscaleX==True:
        axh.set_xscale('log',nonposx='clip')
    if logscaleY==True:
        axh.set_yscale('log',nonposy='clip')
    if svfig==True:
        fgname = '{}_histo.png'.format(cname)
        figh.savefig(fgname,format='png')
    else:
        plt.show()

if __name__=='__main__':
    data_csv = 'planetsFull.csv'
    opendb_csv = 'open_exoplanet_catalogue.csv'
    exop_header = np.array([])
    dtable = readin_table(data_csv,exop_header)
    #print(dtable.iloc[0]) #To call all data in the first row
    #print(dtable['fpl_hostname'].values[:]) #to call all names in col WITHOUT index column
    #Note that hostname is NOT unique, fpl_name is the unique names of planets
    opencat_colnames = np.array(['opl_name','opl_binflag','opl_pmass','opl_prad',
                            'opl_period','opl_sma','opl_ecc','opl_periastron',
                            'opl_longitude','opl_asnode','opl_inc','opl_eqtemp',
                            'opl_sysage','opl_discmethod','opl_discyear','opl_lastupdate',
                            'opl_ra','opl_dec','opl_dist','opl_starmass','opl_starrad',
                            'opl_starmet','opl_startemp','opl_starage','opl_pstatus'])
    dtable2 = readin_table(opendb_csv,opencat_colnames,header_passed=True)
    #print(dtable2.iloc[0])
    #print(dtable2['opl_name'].values[:])
    #print('----------------\n')
    #num_obj_exop = len(np.unique(dtable['fpl_name'].values[:]))
    #num_obj_open = len(np.unique(dtable2['opl_name'].values[:]))
    #print('Num EXP: {} | Num OPEN: {}'.format(num_obj_exop,num_obj_open))
    #test_planet = dtable2['opl_name'].values[1]
    #print(test_planet)
    #dloc = np.where(dtable['fpl_name'].values[:]==test_planet)[0]
    #print(dloc)


    #For the first plot in HW3, part 2:
    x_str1 = r'Semi-major Axis (au)'
    y_str1 = r'Planet Mass (M$_{\textrm{J}}$)'
    scat_plotter(dtable,'pl_orbsmax','pl_bmassj',x_str1,y_str1,logscaleY=False)
    #For the second plot in HW3, part 2:
    x_str2 = r'Period (days)'
    y_str2 = r'Planet Radius (R$_{\textrm{J}}$)'
    scat_plotter(dtable,'pl_orbper','pl_radj',x_str2,y_str2,logscaleY=False)
    #for the third plot:
    x_str3 = r'Period (days)'
    y_str3 = r'Eccentricity'
    scat_plotter(dtable,'pl_orbper','pl_orbeccen',x_str3,y_str3,logscaleY=False)

    #for the fourth plot:
    x_str4 = r'Planet Radius (R$_{\textrm{J}}$)'
    y_str4 = r'Planet Mass (M$_{\textrm{J}}$)'
    scat_plotter(dtable,'pl_radj','pl_bmassj',x_str4,y_str4,logscaleX=False,logscaleY=False,cdense=True)

    #And now for the histograms:
    x_str5 = r'Period (days)'
    histo_plotter(dtable,'pl_orbper',x_str5,50)

    x_str6 = r'Planet Mass (M$_{\textrm{J}}$)'
    histo_plotter(dtable,'pl_bmassj',x_str6,50)
