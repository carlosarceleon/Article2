import article2_functions as art2fnc
reload(art2fnc)
delta = 0.0159
Uinf  = 20.

root = '/media/carlos/6E34D2CD34D29783/2015-02_SerrationPIV/TR_Data_NewProcessing/'
data_root = './PointData/'

#case='Sr20R21_a0_p0_U20_z10_tr.hdf5') 

# Take the Sr20R21 a0 p0 case and try to understand how turbulence
# changes as it convects over the surface (z00) and the root (z10)

def compare_SPD(device="Sr20R21",alpha=0,phi=0,z=10):
    """ Gets the spectral power density for the 'popular' points
    in the root folder (defined in the function)
    """
    from os import listdir
    from os.path import join
    import pandas as pd
    from numpy import argmin,linspace,log10
    import time_series_functions as tdf
    from re import findall
    import matplotlib.pyplot as plt
    import seaborn as sns
    sns.__version__

    all_data_points = sorted(
        [f for f in listdir(data_root) \
         if f.endswith('.p')\
         and "a{0}".format(alpha) in f\
         and "p{0}".format(phi) in f\
         and "z{0:02d}".format(z) in f\
         and device in f
        ]
    )

    px,py = art2fnc.popular_points()

    data_points = art2fnc.check_quick_point_existence(
        all_data_points,
        root,
        data_root
    )

    fig,ax = plt.subplots(len(py),1,sharey=True,figsize=(5,10)) 
    for data in data_points:
        df = pd.read_pickle(join(data_root,data))
        spd_x,freq = tdf.get_spectral_power_density(df.vx)
        spd_y,freq = tdf.get_spectral_power_density(df.vy)
        x = float(findall('px[0-9]+.[0-9]+',data)[0].replace('px',''))
        y = float(findall('py[0-9]+.[0-9]+',data)[0].replace('py',''))
        y_ix = len(py)-1-argmin(map(abs,py-[y]*len(py)))

        St = freq*delta/Uinf

        slope_x = linspace(0.2,1,100)
        slope_y = 10*log10(slope_x**(-5/3.))-40

        spd = 10*log10(spd_x+spd_y)

        if x == px[0]:
            ax[y_ix].plot(slope_x,slope_y,color='k',lw=1)
            ax[y_ix].text(0.6,-35,"$\\mathrm{St}^{-5/3}$")
            ax[y_ix].text(-0.25,0.5,
                          "$y/\\delta = {0:.2f}$".format(y*4/1.5),
                          transform=ax[y_ix].transAxes,
                          ha='center',
                          va='center',
                          rotation=90
                         )
        ax[y_ix].plot(St,spd,label="$x/2h=${0}".format(x))

        ax[y_ix].set_xscale('log')
        ax[y_ix].set_xlim( 5e-2  , 2 )
        ax[y_ix].set_ylim( -55,-20  )
        if not y_ix == len(ax)-1:
            plt.setp( ax[y_ix].get_xticklabels(), visible=False )
        else:
            ax[y_ix].set_xlabel('$\\mathrm{St} = f\\delta/U_\\infty$')
        ax[y_ix].set_ylabel(
            '$10\\,\\mathrm{log}_{10}(\\Phi_x+\\Phi_y)$ \n[dB]'
        )

    ax[0].legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    fig.tight_layout()
    fig.savefig(
        "images/SPD_{0}_a{1}_p{2}_z{3}.png".format(
            device,alpha,phi,z),
        bbox_inches='tight'
    )
    fig.clear()



