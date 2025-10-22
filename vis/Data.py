import numpy as np, matplotlib.pyplot as plt
import tjy_udf as tjy

class Data: 
    def __init__(self, alp, sig, morph, multi=False, done=False):
        self.alpha = alp
        self.sigma = sig
        self.morph = morph
        self.multi = multi
        self.done  = done
        self.data  = []
        
    class freeE:
        ''' Class of data associated with a converged run '''
        def __init__(self, dims, F, icm=None, err=None): 
            self.F = F
            self.icm = icm 
            self.err = err 
            self.dims = dims
            if np.size(dims) == 3: self.lx, self.ly, self.lz = dims[0], dims[1], dims[2] 
            else: print('\n\nDIMS NOT SIZE 3: (lx, ly, lz)\n\n')
            
    def show(self):
        ''' Print summary of contained runs '''
        for i in self.data: print('({:5.2f}, {:5.2f}, {:5.2f}): {:.5f}'.format(*i.dims, i.F))
        
    def new(self, dims, F, icm=None, err=None): 
        ''' Add new converged run '''
        if dims not in [i.dims for i in self.data]:
            self.data.append(self.freeE(dims, F, icm=icm, err=err))
        else:
            self.data.pop([i.dims for i in self.data].index(dims))
            self.data.append(self.freeE(dims, F, icm=icm, err=err))
    
    def pop(self, dims):
        self.data.pop([i.dims for i in self.data].index(dims))
    
    def minF(self, dims=False):
        ''' Get minimum F of all converged runs '''
        fs = [i.F for i in self.data]
        if len(fs) == 0: 
            return np.NaN
        if dims: return (min(fs), self.data[np.argmin(fs)].dims)
        else: 
            return min(fs)
    
    def plotF(self, opt, guess=None, fig=None):
        if not fig: fig,ax = plt.subplots() # Make new figure if one is not passed
        
        if opt == 'A' or opt == 'LxLy' or opt == 'lxly':
            self.data.sort(key=lambda data: data.lx*data.ly)
            plt.plot([i.lx*i.ly for i in self.data], [i.F for i in self.data],'ko-', markerfacecolor='w', markersize=8)
            plt.xlabel(r'$L_xL_y$')
            plt.ylabel(r'$\langle F \rangle$')
            
        elif opt == 'lx':
            self.data.sort(key=lambda data: data.lx)
            plt.plot([i.lx for i in self.data], [i.F for i in self.data],'ko-', markerfacecolor='w', markersize=8)
            plt.xlabel(r'$L_x$')
            plt.ylabel(r'$\langle F \rangle$')
            
        elif opt == 'ly': 
            self.data.sort(key=lambda data: data.ly)
            plt.plot([i.ly for i in self.data], [i.F for i in self.data],'ko-', markerfacecolor='w', markersize=8)
            plt.xlabel(r'$L_y$')
            plt.ylabel(r'$\langle F \rangle$')
            
        elif opt == '2d' or '2D': 
            import warnings
            import scipy.optimize as optimize
            from scipy.optimize import OptimizeWarning
            warnings.simplefilter("ignore", OptimizeWarning)
            
            lxs = np.array([i.lx for i in self.data])
            lys = np.array([i.ly for i in self.data])
            fs  = np.array([i.F  for i in self.data])
            
            def func5(_in, a, b, c, d, e): ## Assume quadratic in x and y, 5 params
                return a* _in[0]*_in[0] + b* _in[0] + c* _in[1]*_in[1] + d* _in[1] + e
            def func3(_in, a, b, c): ## Assume linear in x and y, 3 params
                return a* _in[0] + b* _in[1] + c
            
            if len(fs) >= 5: func = func5
            elif len(fs)< 3: 
                print('\n\n NEED 3 OR MORE POINTS TO FIT 2D FUNCTION \n\n')
                return
            else:            func = func3
            params, pcov = optimize.curve_fit(func, [lxs,lys], fs)
            
            step = 0.01
            xrange = max(lxs) - min(lxs)
            yrange = max(lys) - min(lys)
            X = np.arange(np.floor(min(lxs)-0.25*xrange), np.ceil(max(lxs)+0.25*xrange)+1e-06, step)
            Y = np.arange(np.floor(min(lys)-0.25*yrange), np.ceil(max(lys)+0.25*yrange)+1e-06, step)
            XX, YY = np.meshgrid(X, Y)
            cs = plt.contourf(XX, YY, func([XX, YY], *params))
            plt.plot(lxs, lys, 'ko', markerfacecolor='w', markersize=8)
            plt.plot(lxs[np.argmin(fs)], lys[np.argmin(fs)], 'ko', markerfacecolor='r', markersize=8)
            cbar = fig.colorbar(cs)
            cbar.ax.set_ylabel(r'$\langle F \rangle$')
            tjy.ticks([cbar.ax])
            plt.sca(ax)
            plt.xlabel(r'$L_x$')
            plt.ylabel(r'$L_y$')

        else: print('\n\nINVALID OPTION\n\n')
                
        plt.text(0.5, 0.9, r'$minF = {:.5f}$'.format(self.minF()), transform=plt.gca().transAxes, fontsize=18, horizontalalignment='center')
        plt.title(r'$\alpha = {:.2f},\ \sigma = {:.3f}$'.format(self.alpha, self.sigma))
        tjy.ticks()
        return