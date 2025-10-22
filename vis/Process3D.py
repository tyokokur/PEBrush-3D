import numpy as np, pandas as pd, matplotlib.pyplot as plt, matplotlib as mpl, plotly.graph_objects as go 
import os.path
from Pha3D import Pha3D

class Process3D:
    ''' 
    kwargs : { fname, dims, fprefix, ( default: isomin, n_coarse, zmax, reflect_over ) }
    '''
    def __init__(self, **kwargs):
        pha3d = Pha3D(**kwargs)
        self.data = pha3d.get_vol(**kwargs)
        self.dims = kwargs['dims']
        self.zmax = kwargs['zmax']
        self.reflect=kwargs['reflect_over']
        print('Data processed into self.data')
            
    def Set_camera(self, opt='d', d = 3.0, a = 25, zcenter = 0.0, old=None, eye_kws=None):
        '''
        'd'istance, 'x'z, 'y'z, 'custom'
        Distance from center
        Angle, if applicable
        Camera focal point
        '''
        
        match opt:
            case 'd': # distance-based
                if not old: 
                    raise ValueError("Need to provide old = {x,y,z} of reference")
                from math import sqrt, atan, acos, cos, sin
                dd = sqrt(old['x']**2 + old['y']**2 + old['z']**2)
                print("Old distance was {:.3f}".format(dd))
                th, ph = atan(old['y']/old['x']), acos(old['z']/dd)
                
                x = d * cos(th)*sin(ph)
                y = d * sin(th)*sin(ph)
                z = d * cos(ph)
                
                print("New: ({:.3f}, {:.3f}, {:.3f})".format(x, y, z))
                eye=dict(x=x, y=y, z=z)
            case 'x': eye=dict(x=0, y=-d*np.cos(np.pi/180*a), z=d*np.sin(np.pi/180*a)) # xz at given angle
            case 'y': eye=dict(x=d*np.cos(np.pi/180*a), y=0, z=d*np.sin(np.pi/180*a)) # yz at given angle
            case 'custom': 
                eye = eye_kws
            case  _ : eye=dict(x=-d/np.sqrt(3), y=-d/np.sqrt(3), z=d/np.sqrt(3)) #default
        
        self.camera = dict(
            up=dict(x=0, y=0, z=1), center=dict(x=0, y=0, z= zcenter),
            eye = eye
        )
        print('Camera set to self.camera')

    def Figure(self, 
               xrange=[], xticks=[], 
               yrange=[], yticks=[],
               zrange=[], zticks=[],
               write_html=True, fprefix = None,
               width = None
               ):
        if not width: width = 1200
        if not xrange: 
            if 'w' in self.reflect: xrange = [-self.dims[0]-0.01, self.dims[0]+0.01]
            else:                   xrange = [-0.01, 2*self.dims[0]+0.01]
        if not yrange: 
            if 'n' in self.reflect: yrange = [-0.01, 2*self.dims[1]+0.01]
            else:                   yrange = [-self.dims[1]-0.01, self.dims[1]+0.01]
        if not zrange: zrange = [0, self.zmax+0.2]
        
        xn, yn, zn = len(xticks), len(yticks), len(zticks)

        try: camera = self.camera
        except NameError: 
            print('Camera not set, doing default')
            self.Set_camera()

        scene = dict(aspectmode='data')

        if xn > 0:
            scene.update( 
                xaxis = dict(nticks=xn, range=xrange, title='',
                             ticks='outside',tickangle=0,tickwidth=1,ticklen=5,
                             ticktext= ['']*xn, tickvals= xticks, 
                             spikethickness=0,
                             zeroline=False, showgrid=False, backgroundcolor='white',
                             )
            )
        else: 
            scene.update(
                xaxis = dict(showticklabels=False, title='', zeroline=False, showgrid=False, backgroundcolor='white',), 
            )

            
        if yn > 0:
            scene.update( 
                yaxis = dict(nticks=yn, range=yrange, title='',
                             ticks='outside',tickangle=0,tickwidth=1,ticklen=5,
                             ticktext= ['']*yn, tickvals= yticks, 
                             spikethickness=0,
                             zeroline=False, showgrid=False, backgroundcolor='white',
                             )
            )
        else: 
             scene.update(
                 yaxis = dict(showticklabels=False, title='', zeroline=False, showgrid=False, backgroundcolor='white',), 
             )

        if zn > 0:
            scene.update(
                zaxis = dict(nticks=zn, range=zrange, title='',
                             ticks='outside',tickangle=0,tickwidth=1,ticklen=5,
                             ticktext= ['']*zn, tickvals= zticks, 
                             spikethickness=0,
                             zeroline=False, showgrid=False,backgroundcolor='white',
                             ) 
            )
        else: 
            scene.update(
                zaxis = dict(showticklabels=False, title='', zeroline=False, showgrid=False, backgroundcolor='white',), 
            )

        
        layout = go.Layout(
            autosize=False,
            font= dict(size=24, family='serif', color='black'),
            hovermode=False,
            scene = scene,
            scene_camera=self.camera,
            width=width,
            height=width/1.5,
            margin=dict(l=20, r=20, t=20, b=20)
        )
        fig = go.Figure(self.data, layout)
        fig.update_traces(showscale=False)

        if write_html:
            fname = fprefix+'pha_vol.html'
            f = open(fname, "w")
            f.close()
            with open(fname, 'a') as f:
                f.write(fig.to_html(full_html=False, include_plotlyjs='cdn'))
            f.close()
            print('Wrote to {}'.format(fname))
                
            import webbrowser
            webbrowser.open(fname)
            print('Opening...')
        
        return fig

class Local:
    def __init__(self):
        self.files = []
    
    def load(self, fname, dims, params=(0,0), figs=(), fprefix=None, quiet=True):
        self.files.append(self.file(fname, dims, params=params, figs=figs,fprefix=fprefix,quiet=quiet))
        
    def search(self, search_kws, quiet=True):
        ret = []
        kfigs = search_kws.pop('figs', None)
        for i in self.files:
            if all( [getattr(i, k)==search_kws[k] for k in search_kws.keys()] ):
                if kfigs: 
                    if kfigs in getattr(i, 'figs'):
                        ret.append(i)
                else: 
                    ret.append(i)
        if not quiet: 
            print("Search results: ", [i.fname for i in ret])
        return ret
    
    class file:
        def __init__(self, fname, dims, params=(0,0), figs=(), fprefix=None, quiet=True):
            if type(fname)==str : 
                self.fname = fname
            else: 
                raise TypeError("fname should be filename as string")
            if not os.path.isfile(fprefix+self.fname): 
                raise TypeError("{} not found".format(fprefix+self.fname))
            if type(dims)==tuple and len(dims)==3 :  
                self.dims = dims
            else: 
                raise TypeError("dims should be tuple of len 3")
            if type(params)==tuple and len(params)==2 : 
                self.params= params
                self.alpha = params[0]
                self.sigma = params[1]
            else: 
                raise TypeError("params should be tupe of len 2")
            if all(isinstance(n, int) for n in figs) : 
                self.figs = figs
            else: 
                raise TypeError("figs should be tuple/list of ints")
            if type(quiet)==bool: 
                if not quiet : print("{} recorded.".format(fname))
            else: 
                raise TypeError("quiet should be bool")