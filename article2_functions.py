root = "/media/carlos/6E34D2CD34D29783/2015-02_SerrationPIV/TR_Data/"

cases = [
    "Slit20R21_a0_p0_U20_z00_tr.h5"  ,
    "Sr20R21_a0_p0_U20_z05_tr.h5"    ,
    "Slit20R21_a0_p0_U20_z05_tr.h5"  ,
    "Sr20R21_a0_p0_U20_z10_tr.h5"    ,
    "Slit20R21_a0_p0_U20_z10_tr.h5"  ,
    "Sr20R21_a12_p0_U20_z00_tr.h5"   ,
    "Slit20R21_a12_p0_U20_z00_tr.h5" ,
    "Sr20R21_a12_p0_U20_z05_tr.h5"   ,
    "Slit20R21_a12_p0_U20_z05_tr.h5" ,
    "Sr20R21_a12_p0_U20_z10_tr.h5"   ,
    "Slit20R21_a12_p0_U20_z10_tr.h5" ,
    "STE_a0_p0_U20_z00_tr.h5"        ,
    "Sr20R21_a0_p0_U20_z00_tr.h5"    ,
    "STE_a12_p0_U20_z00_tr.h5"       ,
]

def load_time_series(device="Sr20R21",alpha='0',phi='0',z='10', p=(-0.1,0.1)):
    from time_data_functions import read_hdf5_time_series
    import os

    #case = "{0}_a{1}_p{2}_U20_z{3}_tr_planar".format(
    #        device,
    #        alpha,
    #        phi,
    #        z
    #)
    case = "{0}_a{1}_p{2}_U20_z{3}_tr_NewProcessing".format(
            device,
            alpha,
            phi,
            z
    )

    hdf5_file = os.path.join(
        root,
        'TimeData_NewProcessing.hdf5'
        #'planar.hdf5'
    )
    
    return read_hdf5_time_series(
            hdf5_file,
            case,
            loc = p
        ).interpolate()

def check_PSD(device="Sr20R21",alpha='0',phi='0',z='10',
              p=(-0.1,0.1),component='vx',plot_name='test_psd.png',
             sampling=10000,NFFT=256):
    from time_series_functions import get_spectral_power_density#,\
            #butter_lowpass_filter
    from matplotlib import pyplot as plt
    import seaborn as sns
    sns.__version__

    time_series = load_time_series(device=device,alpha=alpha,phi=phi,z=z,p=p)

    #time_series.vx = butter_lowpass_filter(time_series.vx,cutoff=2000,fs=10000)
    #time_series.vy = butter_lowpass_filter(time_series.vy,cutoff=2000,fs=10000)
    #time_series.vz = butter_lowpass_filter(time_series.vz,cutoff=2000,fs=10000)
    
    Pxx,freqs = get_spectral_power_density(time_series.vx,
                                           NFFT=NFFT,sampling=sampling)
    Pyy,freqs = get_spectral_power_density(time_series.vy,
                                           NFFT=NFFT,sampling=sampling)
    #Pzz,freqs = get_spectral_power_density(time_series.vz,
    #                                       NFFT=NFFT,sampling=sampling)

    if plot_name:
        fig = plt.figure()

        plt.plot(freqs, Pxx, label='$u$', alpha=0.6)
        plt.plot(freqs, Pyy, label='$v$', alpha=0.6)
        #plt.plot(freqs, Pzz, label='$z$', alpha=0.6)

        plt.ylabel("Power spectral density [p$^2/$Hz]")
        plt.xlabel("Frequency [Hz]")
        plt.yscale('log')
        plt.xscale('log')
        plt.xlim(xmin=50,xmax=5000)
        #plt.ylim(ymin=10e-4)
        plt.legend(loc='upper right')
        plt.savefig(plot_name)
        fig.clear()
    return [Pxx,Pyy,freqs]


def check_autocorrelation(device="Sr20R21",alpha='0',phi='0',z='10',
                          p=(-0.1,0.1),component='vx',plot_name='test.png'):
    from time_series_functions import get_autocorrelation,butter_lowpass_filter
    from matplotlib import pyplot as plt
    import seaborn as sns
    from numpy import arange
    sns.__version__

    time_series = load_time_series(device=device,alpha=alpha,phi=phi,z=z,p=p)

    time_series.vx = butter_lowpass_filter(time_series.vx,cutoff=2000,fs=10000)
    time_series.vy = butter_lowpass_filter(time_series.vy,cutoff=2000,fs=10000)
    time_series.vz = butter_lowpass_filter(time_series.vz,cutoff=2000,fs=10000)

    autocorr_u = get_autocorrelation(
        time_series.vx
    )
    autocorr_v = get_autocorrelation(
        time_series.vy
    )
    autocorr_w = get_autocorrelation(
        time_series.vz
    )

    if plot_name:
        fig = plt.figure()
        plt.plot(arange(len(autocorr_u)),
                 autocorr_u/autocorr_u.max(),label='$u$',alpha=0.6)
        plt.plot(arange(len(autocorr_u)),
                 autocorr_v/autocorr_v.max(),label='$v$',alpha=0.6)
        plt.plot(arange(len(autocorr_u)),
                 autocorr_w/autocorr_w.max(),label='$w$',alpha=0.6)
        plt.xlabel("Lag [time steps (1:1/10000 s)]")
        plt.ylabel("Autocorrelation")
        plt.xlim(0,20) 
        plt.legend(loc='upper right')
        plt.savefig(plot_name)
        fig.clear()

    return autocorr_u,autocorr_v

def plot_time_series(device="Sr20R21",alpha='0',phi='0',z='10',
                     p=(-0.1,0.1),component='vx',plot_name='test.png'):
    from matplotlib import pyplot as plt
    import seaborn as sns
    from time_series_functions import butter_lowpass_filter
    sns.__version__

    time_series = load_time_series(device=device,alpha=alpha,
                                   phi=phi,z=z,p=p)

    
    time_series_low_passed = butter_lowpass_filter(time_series.vx,
                                                   cutoff=2000,fs=10000)

    fig = plt.figure()
    plt.plot(time_series.t,time_series.vx,label='$u$',alpha=0.6)
    plt.plot(time_series.t,time_series_low_passed,
             label='Low passed $u$',alpha=1.0,lw=3,color='k')
    #plt.plot(time_series.t,time_series.vy,label='$v$',alpha=0.6)
    #plt.plot(time_series.t,time_series.vz,label='$w$',alpha=0.6)
    plt.xlabel("t [s]")
    plt.ylabel("Velocity [m/s]")
    plt.xlim(0,500/10000.)
    plt.ylim(-5,22)
    plt.legend(loc='lower right')
    plt.savefig(plot_name)
    fig.clear()

def plot_surface_at_t(device="Sr20R21",alpha='0',phi='0',z='10',
                      p=[(-0.1,0.1)],component='Vx',plot_name='test.png',
                      t=0):
    from time_data_functions import read_time_series_range
    from matplotlib import pyplot as plt
    from numpy import meshgrid
    import os

    points = p

    case = "{0}_a{1}_p{2}_U20_z{3}_tr_NewProcessing".format(
            device,
            alpha,
            phi,
            z
    )

    hdf5_file = os.path.join(
        root,
        'TimeData_NewProcessing.hdf5'
    )
    df = read_time_series_range(
        hdf5_file = hdf5_file,
        case      = case,
        variable  = "Vx",
        ti        = t,
        tf        = t+1
    ).interpolate()

    return df

    X,Y = meshgrid( df.x.unique(), df.y.unique() )
    print X.shape
    print df[df.t==t].vx.shape,df[df.t==t].x.shape
    #U   = df[df.t==t].vx.reshape(X.shape)
    #V   = df[df.t==t].vy.reshape(X.shape)
    #W   = df[df.t==t].vz.reshape(X.shape)

    fig = plt.figure()
    ax = plt.subplot(111,aspect=1)
    #levels = list(linspace(0,20,21))+[25]
    #cnt = ax.contourf(Y,-X,U,levels=levels)
    #ax.quiver( Y[::6,::3], -X[::6,::3], U[::6,::3], V[::6,::3],
    #          linewidths=(1,), edgecolors=('k'), scale=700 )
    for p in points:
        ax.scatter(p[1],-p[0],marker='x',s=300,color='k')
        ax.scatter(p[1],-p[0],marker='o',s=200,color='k')
    ax.fill_between([df.y.min(),df.y.max()],-df.x.max(),0,facecolor='k')
    plt.savefig(plot_name,bbox_inches='tight')
    fig.clear()

    return df

def popular_points():
    from numpy import linspace
    px = linspace(0,1,5)
    # 2h = 4.0cm; the BL is about 1.5 cm, so 0.375*2h, plus ignore
    # the first 0.4 cm, so it starts at 0.1
    py = linspace(0.1,0.375,6)
    return px,py

def build_popular_point_matrix(root,case,save_folder=0):
    """ Get a 5x5 matrix of the most important (used) points 
    for processing, and save them as pickles of the time series

    Input
        case:           the HDF5 case
        save_folder:    folder where to save the pickle
    """

    if not save_folder:
        save_folder = root

    px,py = popular_points()

    for x in px:
        for y in py:
            make_hdf5_time_series(root=root,
                                  case=case,
                                  x=x,
                                  y=y,
                                  save_folder=save_folder
                                 )

def make_hdf5_time_series(root,case,x,y,save_folder):
    import time_data_functions as tdf
    from os.path import join
    df = tdf.read_hdf5_time_series(join(root,case),
                              case.replace('.hdf5',''),
                              loc=(-y,x)
                             )
    if df is not None:
        df['x'] = x
        df['y'] = y
        df.to_pickle(join(save_folder,"{0}_px{1}_py{2}.p".format(
            case.replace('.hdf5',''),
            x,
            y
        )))

def check_quick_point_existence(requested_data_points,root,data_target):
    """ This function checks if the points requested exist amongst
    the quick points popular points matrix, otherwise go and make them

    Input: a file array of the available points for this device and
        conditions

    Returns: the files that have been approved for processing

    """
    from numpy import argmin
    from re import findall

    px,py = popular_points()
    data_points  = []
    non_existant = []
    for all_data in requested_data_points:
        x = float(findall('px[0-9]+.[0-9]+',all_data)[0].replace('px',''))
        y = float(findall('py[0-9]+.[0-9]+',all_data)[0].replace('py',''))
        dy = abs(y-py)
        dx = abs(x-px)
        if dx[argmin(dx)]<1e-5 and dy[argmin(dy)]<1e-5:
            data_points.append(all_data)
        else:
            non_existant.append(all_data)

    if len(non_existant):
        for n_existant in non_existant:
            case = findall("^[A-Za-z0-9.-]_tr",n_existant)[0]+".hdf5"
            make_hdf5_time_series(root,case,x,y,data_target)
            data_points.append(n_existant)

    return data_points


#p = (-0.05,0.3)
#check_PSD(p=p,NFFT=512)
#acorr = check_autocorrelation(plot_name='test_autocorrelation_new.png')
#plot_time_series(plot_name='test_timeseries_new.png')
#df = plot_surface_at_t(p=p,plot_name='test_surface.png')
#for i in range(100):
#    df = plot_surface_at_t(t=i,plot_name='test_{0:03d}.png'.format(i))

