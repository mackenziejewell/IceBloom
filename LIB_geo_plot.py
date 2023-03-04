#////////////////////////////
#  fix_cartopy_vectors   ///
#//////////////////////////
#---------------------------------------------------------------------
# function to output vector components for plotting in cartopy. 
#---------------------------------------------------------------------
#////////////////
#  add_land  ///
#//////////////
#---------------------------------------------------------------------
# Add land feature to cartopy figure
#---------------------------------------------------------------------
#/////////////////
#  add_coast  ///
#///////////////
#---------------------------------------------------------------------
# Add coast feature to cartopy figure
#---------------------------------------------------------------------
#////////////////
#  add_grid  ///
#//////////////
#---------------------------------------------------------------------
# Add specified gridlines to cartopy figure
#---------------------------------------------------------------------
#////////////////
#  add_date  ///
#//////////////
#---------------------------------------------------------------------
# Add date label to cartopy plot.
#---------------------------------------------------------------------


#////////////////////////////
#  fix_cartopy_vectors   ///
#//////////////////////////
#---------------------------------------------------------------------
# function to output vector components for plotting in cartopy. 
#---------------------------------------------------------------------
# DEPENDENCIES:
import numpy as np
import numpy.ma as ma
#---------------------------------------------------------------------
def fix_cartopy_vectors(u, v, uvlats):
    
    """Function to output vector components for plotting in cartopy. 
    
    Reads in vectors and associated latitudes, return fixed vectors. 
    
    Cartopy doesn't know meters per degree increase toward the pole 
    in zonal direction as cosine(lat) when reprojecting for vectors 
    given as m/s (where  a meter covers different amount of a degree 
    depending on latitude), we need to rescale the u (east-west) speeds. 
    otherwise at high latitudes, u will be drawn much larger than 
    they should which will give incorrect angles

INPUT: 
- u: (N x M) array of eastward velocity component (m/s)
- v: (N x M) array of northward velocity component (m/s)
- uvlats: (N x M) array latitudes associated with u,v vectors

OUTPUT:
- u_fixed: (N x M) array of u with correct angle and magnitude for plotting
- v_fixed: (N x M) array of v with correct angle and magnitude for plotting


DEPENDENCIES:
import numpy as np
import numpy.ma as ma

Latest recorded update:
04-20-2022
    """
    
    # FIX ANGLE
    # fix u scale to be correct relative to v scale
    #----------------------------------------------
    # for u (m/s), m/degree decreases as cos(lat)
    u_fixed = u/np.cos(uvlats/180*np.pi) 
    v_fixed = v  # v does not change with latitude

    # FIX MAGNITUDE
    # scale fixed u,v to have correct magnitude 
    #-----------------------
    # original magnitude 
    orig_mag = ma.sqrt(u**2+v**2)
    # new magnitude
    fixed_mag = ma.sqrt(u_fixed**2+v_fixed**2)
    u_fixed = u_fixed*(orig_mag/fixed_mag)
    v_fixed = v_fixed*(orig_mag/fixed_mag)

    return u_fixed, v_fixed


#////////////////
#  add_land  ///
#//////////////
#---------------------------------------------------------------------
# Add land features to cartopy figure
#---------------------------------------------------------------------
# DEPENDENCIES:
#-------------
import numpy as np
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeat
import matplotlib as mpl
from matplotlib import pyplot as plt
import matplotlib.colors
from shapely import wkt
#---------------------------------------------------------------------
def add_land(ax, scale = '50m', color='gray', alpha=1, fill_dateline_gap = True, zorder=2):
    
    """Add land feature to cartopy figure
    
INPUT:
- ax: cartopy figure axis
- scale = NaturalEarthFeature land feature scale (e.g. '10m', '50m', '110m')
        (default: '50m')
- color = land color (e.g. 'k' or [0.9,0.6,0.5]) (default: 'gray')
- alpha = land opacity (default: 1)
- zorder: drawing order of land layer (default: 2)
- fill_dateline_gap: specify whether to fill gap in cartopy land feature along 
   dateline that crosses Russia and Wrangel Island (default: True)

OUTPUT:
- input plot with added land layer

DEPENDENCIES:
import numpy as np
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeat
import matplotlib as mpl
from matplotlib import pyplot as plt
from shapely import wkt

Latest recorded update:
05-11-2022

    """
    
        
    # grab land from cfeat.NaturalEarthFeature
    #-----------------------------------------
    ax.add_feature(cfeat.NaturalEarthFeature(category='physical', name='land', 
                                             scale=scale, facecolor=color),
                                             alpha = alpha, zorder = zorder)

    # if specified, fill dateline gap in land feature with shapely polygons
    if fill_dateline_gap == True:
        # generate polygon to fill line across Wrangel Island and line across Russia
        WKT_fill_Wrangel = 'POLYGON ((-180.1 71.51,-180.1 71.01,-179.9 71.01,-179.9 71.51,-180.1 71.51))'
        poly1 = wkt.loads(WKT_fill_Wrangel)
        ax.add_geometries([poly1], crs=ccrs.PlateCarree(), 
              facecolor=color, edgecolor=color, alpha = alpha, zorder=zorder)
        WKT_fill_Russia = 'POLYGON ((-180.1 65.1,-180.1 68.96,-179.9 68.96,-179.9 65.1,-180.1 65.1))'
        poly2 = wkt.loads(WKT_fill_Russia)
        ax.add_geometries([poly2], crs=ccrs.PlateCarree(), 
              facecolor=color, edgecolor=color, alpha = alpha, zorder=zorder)


        
#/////////////////
#  add_coast  ///
#///////////////
#---------------------------------------------------------------------
# Add coast feature to cartopy figure
#---------------------------------------------------------------------
# DEPENDENCIES:
import numpy as np
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeat
import matplotlib as mpl
from matplotlib import pyplot as plt
import matplotlib.colors
#---------------------------------------------------------------------
def add_coast(ax, scale = '50m', color='gray', linewidth = 1, alpha=1, zorder=3):

    """Add land feature to cartopy figure
    
INPUT:
- ax: cartopy figure axis
- scale = NaturalEarthFeature coast feature scale (e.g. '10m', '50m', '110m')
        (default: '50m')
- color = coastline color (e.g. 'k' or [0.9,0.6,0.5]) (default: 'gray')
- linewidth = coastline linewidth (default: 1)
- alpha = coastline opacity (default: 1)
- zorder: drawing order of coast layer (default: 3)

OUTPUT:
- input plot with added coast layer

DEPENDENCIES:
import numpy as np
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeat
import matplotlib as mpl
from matplotlib import pyplot as plt
from shapely import wkt

Latest recorded update:
05-06-2022

    """

    # coastline
    #----------
    ax.coastlines(scale, color=color, linewidth=linewidth, alpha = alpha, zorder = zorder)



#////////////////
#  add_grid  ///
#//////////////
#---------------------------------------------------------------------
# Add specified gridlines to cartopy figure
#---------------------------------------------------------------------
# DEPENDENCIES
import numpy as np
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeat
import matplotlib.ticker as mticker
import matplotlib as mpl
from matplotlib import pyplot as plt
#---------------------------------------------------------------------

def add_grid(ax, lats = None, lons = None, linewidth = 1, color = 'gray', alpha=0.5, zorder = 4): 
    
    """Add specified gridlines to cartopy figure.
    
INPUT:
- ax: cartopy figure axis
- lats: None or array of latitudes to plot lines (default: None)
- lons: None or array of latitudes to plot lines (default: None)
- linewdith: grid line linewidths (default: 1)
- color: grid line color (default: 'gray')
- alpha: line transparency (default: 0.5)
- zorder: drawing order of gridlines layer (default: 4)

OUTPUT:
- input plot with added grid

DEPENDENCIES:
import numpy as np
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeat
import matplotlib.ticker as mticker
import matplotlib as mpl
from matplotlib import pyplot as plt

Latest recorded update:
04-22-2022

    """
        
    # give gridline specifications
    #-----------------------------
    gl = ax.gridlines(crs=ccrs.PlateCarree(), linewidth=linewidth, color=color, alpha=alpha, zorder = zorder)

    # add the longitude gridlines
    #----------------------------
    if lons is None:
        gl.xlocator = mticker.FixedLocator([])
    else:
        # shift all longitudes from [0,360] to [180,-180]
        lons = np.concatenate((lons[(lons>180)]-360,lons[(lons<=180)]))
        gl.xlocator = mticker.FixedLocator(lons)

        
    # add the latitude gridlines
    #----------------------------
    if lats is None:
        gl.ylocator = mticker.FixedLocator([])
    else:
        gl.ylocator = mticker.FixedLocator(lats)
        




        

#////////////////
#  add_date  ///
#//////////////
#---------------------------------------------------------------------
# Add date label to cartopy plot.
#---------------------------------------------------------------------
# DEPENDENCIES
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
import matplotlib.colors
from datetime import datetime
from matplotlib.offsetbox import AnchoredText
#---------------------------------------------------------------------

def add_date(fig, ax, dt_obj, date_format = '%b %d, %Y (%H:%M UTC)', method = 'anchor',
             boxstyle="round,pad=0.,rounding_size=0.2", facecolor = 'black', edgecolor = 'black',
             zorder = 10,
             
             anchor_loc = 4, anchor_prop = {'size': 20, 'color':'white'},
             

             x = 0.02, y= 0.05, textcolor = 'white',fontsize=15): 
    
    """Add date label to cartopy plot.
    
INPUT:
- fig: cartopy figure
- ax: cartopy figure axis
- dt_obj: datetime object of date for plotted data 
            OR
          string with text to show (date format already provided (e.g. 'Dec 20, 2018 (6:00 UTC)')
          
IF dt_obj IS DATETIME OBJECT:
- date_format: str, format to display date (default: '%b %d, %Y (%H:%M UTC)')
    - example 1: '%b %d, %Y (%H:%M UTC)' could give 'Dec 20, 2018 (6:00 UTC)'
    - example 2: '%m-%d-%Y' could give '12-20-2018'
    
- method: method to place the date label (either 'anchor' for AnchoredText or 'manual' to place manually).
        (default: 'anchor')
- boxstyle: anchor box shape style (default: "round,pad=0.,rounding_size=0.2")
- facecolor: color of bounding box (default: 'black')
- edgecolor: color of bounding box edge (default: 'black')
- zorder: drawing order of date layer (default: 10)

IF METHOD = 'anchor':
- anchor_loc: anchor text location (default: 4)
- anchor_prop: anchor properties dictionary (default: {'size': 20, 'color':'white'})

IF METHOD = 'manual':
- x: x-location of figure extent to place date
- y: y-location of figure extent to place date
- textcolor: color oftext (default: 'white')
- fontsize: fontsize of text (defult: 15)

OUTPUT:
- input plot with added date label

DEPENDENCIES:
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
import matplotlib.colors
from datetime import datetime
from matplotlib.offsetbox import AnchoredText

Latest recorded update:
06-03-2022
    """

    
    assert method in ['anchor', 'manual'], f">>> method should be 'manual' or 'anchor', given: '{method}'"
    
    assert str(type(dt_obj)) in ["<class 'datetime.datetime'>", "<class 'str'>"], f">>> dt_obj should be datetime object or string, given: {str(type(dt_obj))}"
    
    
    # if given as datetime object, convert to specified date format
    if str(type(dt_obj)) == "<class 'datetime.datetime'>":
        date_text = dt_obj.strftime(date_format)
    
    # else, set date directly to given string object
    else:
        date_text = dt_obj
    
    
    # add text
    #---------
    if str(method) == 'anchor':
        at = AnchoredText(date_text, loc=anchor_loc, prop=anchor_prop)
        at.patch.set_boxstyle(boxstyle)
        at.patch.set_facecolor(facecolor)
        at.patch.set_edgecolor(edgecolor)
        at.zorder = zorder
        ax.add_artist(at)
    
    elif str(method) == 'manual':
        ax.text(x, y, date_text, 
                bbox=dict(boxstyle = boxstyle, facecolor=facecolor, edgecolor = edgecolor), 
                transform=ax.transAxes, fontsize=fontsize, 
                c=textcolor, verticalalignment='top', zorder = zorder);

        
        
