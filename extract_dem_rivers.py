#!/usr/bin/env python3
#this is the python library used to extracted river network from DEM
from pylib import *
from pysheds.grid import Grid
from affine import Affine

def read_deminfo(path,name='dem',size_subdomain=5e5, **kwargs):
    if path.endswith('tif'):
        gd=SGrid.from_raster(path, data_name, **kwargs);
    elif path.endswith('asc'):
        gd=dem_data(); gd.add_deminfo(path,name,size_subdomain=size_subdomain);
        
        # gd=SGrid.from_ascii(path, data_name, header_only=True,**kwargs)
        #gd=SGrid.from_ascii(path, data_name, header_only=header_only,window=window,**kwargs);           
    return gd 

class dem_data(object):
    def __init__(self):
        self.name=[]
        pass
    
    def read_dem_data(self, name=None,index_domain=None,dem_data=None):
        #read dem data for each subdomain        
        if name==None:
           sname=self.name
        elif type(name)==str:
           sname=[name,]
        elif type(name)==list:
           sname=sname
        else:
           sys.exit('unknow name type')
                   
                
        for snamei in sname:
            S=npz_data();exec('S.dinfo=self.{}_info'.format(snamei)); dinfo=S.dinfo;
            
            #read sample              
            if dem_data is not None:
                data0=SGrid.from_ascii(dinfo.path,'data',skiprows=6,usecols=arange(5),max_rows=5)
            
            #get indices of subdomain
            if index_domain==None:
                sind=arange(dinfo.nsubdomain)
            elif type(index_domain)==list:
                sind=index_domain
            else:
                sind=[index_domain]
             
            #read data for each subdomain
            data=[]
            for sindi in sind:
                dm=dinfo.domains[sindi]
                
                if dem_data is None:
                    datai=SGrid.from_ascii(dm.path,'data',skiprows=dm.skiprows,usecols=dm.usecols,max_rows=dm.max_rows, affine_new=dm.affine);
                else:           
                    datai=data0;
                    dem=dem_data[dm.ixy[0]:dm.ixy[1],dm.ixy[2]:dm.ixy[3]]                 
                    datai.add_gridded_data(data=dem, data_name='data', affine=dm.affine, shape=dm.shape,
                               crs=data0.crs, nodata=data0.nodata, metadata={}) 
                    datai.shape=dm.shape; datai.mask=ones(dm.shape)==1
                  
                data.append(datai)
            exec('self.{}_data=data'.format(snamei))
                            
    def add_deminfo(self,path,name,size_subdomain=5e6):
        # 1G~=6e7 pts; 100M~=5e6 pts
                
        #read *asc info
        gd=SGrid.from_ascii(path,'name',header_only=True)
        
        #dem information
        nrows=gd.shape[0]; ncols=gd.shape[1]; 
        xll=gd.affine[2]; yll=gd.affine[5]; cellsize=gd.affine[0]
        
        #save dem information 
        deminfo=npz_data(); deminfo.type='ascii'; deminfo.path=path
        deminfo.ncols=ncols; deminfo.nrows=nrows; deminfo.xll=xll; deminfo.yll=yll; deminfo.cellsize=cellsize; 
        deminfo.affine=gd.affine; deminfo.bbox=gd.bbox; deminfo.extent=gd.extent
        
        #divide the domain into subdomains
        ndem0=max(around(gd.size/size_subdomain),1); nx=int(round(sqrt(ndem0))); ny=int(ndem0/nx); ndem=nx*ny
        deminfo.nsubdomain=ndem; deminfo.nx=nx; deminfo.ny=ny; deminfo.size=gd.size; deminfo.domains=[]; 
    
        dy0=int(floor(nrows/ny)); dx0=int(floor(ncols/nx))
        offset=int(min(dy0/2,dx0/2,500))
        for i in arange(ny):
            #subdomain index
            if ny==1:
                dy=dy0; iy=0; ylli=yll; indy=[0,dy];                                
            else:
                if i==0:
                    dy=dy0+offset; iy=0; ylli=yll; indy=[0,dy-offset]
                elif i==(ny-1):
                    dy=nrows-(ny-1)*dy0+offset; iy=i*dy0-offset; ylli=yll-iy*cellsize; indy=[offset,dy]
                else:
                    dy=dy0+2*offset; iy=i*dy0-offset; ylli=yll-iy*cellsize; indy=[offset,-offset]                    
            
            for j in arange(nx):
                #subdomain index
                if nx==1:
                    dx=dx0; ix=0; xlli=xll; indx=[0,dx]                               
                else:
                    if j==0:
                        dx=dx0+offset; ix=0; xlli=xll; indx=[0,dx-offset]
                    elif j==(nx-1):
                        dx=ncols-(nx-1)*dx0+offset; ix=j*dx0-offset; xlli=xll+ix*cellsize; indx=[offset,dx]
                    else:
                        dx=dx0+2*offset; ix=j*dx0-offset; xlli=xll+ix*cellsize; indx=[offset,-offset]    
                
                #save subdomain info
                datai=npz_data()
                
                datai.path=path
                datai.affine=Affine(cellsize,0,xlli,0,-cellsize,ylli)        
                datai.bbox=(xlli,ylli-(dy-1)*cellsize,xlli+(dx-1)*cellsize,ylli)
                datai.extent=(xlli,xlli+(dx-1)*cellsize,ylli-(dy-1)*cellsize,ylli)
                datai.size=dx*dy
                datai.shape=(dy,dx)
                datai.skiprows=6+iy
                datai.usecols=range(ix,ix+dx)
                datai.max_rows=dy
                datai.rind=[*indy,*indx]    #extract data from subdomain
                datai.ixy=[iy,iy+dy,ix,ix+dx] #extract data for subdomain
                
                deminfo.domains.append(datai)
          
        exec('self.{}_info=deminfo'.format(name))
        self.name.append(name)
    
class SGrid(Grid): 
        
    @classmethod
    def from_ascii(cls, path, data_name,header_only=False,affine_new=None,**kwargs):          
          newinstance = cls()
          newinstance.read_ascii(path, data_name,header_only=header_only,affine_new=affine_new,**kwargs)          
          return newinstance
                  
    #define new object for adding new method to original Grid Obect
    def read_ascii(self, data, data_name, skiprows=6, crs=Proj('epsg:4326'),
                    xll='lower', yll='lower', metadata={},header_only=False, affine_new=None, **kwargs):
        #copy from pyshed's function
         import ast                               
         with open(data) as header:
             ncols = int(header.readline().split()[1])
             nrows = int(header.readline().split()[1])
             xll = ast.literal_eval(header.readline().split()[1])
             yll = ast.literal_eval(header.readline().split()[1])
             cellsize = ast.literal_eval(header.readline().split()[1])
             nval=header.readline().split()[1];
             if nval.lower() in ('nan','-nan'):
                 nodata=nan
             else:
                 nodata = ast.literal_eval(nval)
             shape = (nrows, ncols)
         if header_only:
             data=array([]);
         else:
             data = np.loadtxt(data, skiprows=skiprows, **kwargs)
             nodata = data.dtype.type(nodata)         
         if affine_new==None:
             affine = Affine(cellsize, 0, xll, 0, -cellsize, yll + (nrows-1)*cellsize)
         else:
             affine=affine_new; shape=data.shape
         self.add_gridded_data(data=data, data_name=data_name, affine=affine, shape=shape,
                               crs=crs, nodata=nodata, metadata=metadata)
         
class dem_dir(object):
    def __init__(self):
        pass
    
    def read_data(self,path):
        #read saved dir data
        data=loadz(path)
        
        #get variable name
        svars=list(data.__dict__.keys())
        if 'VINFO' in svars: svars.remove('VINFO')
        
        #put variables in self
        for svar in svars:
            exec('self.{}=data.{}'.format(svar,svar))
            
    def save_data(self,fname,svars):
        #function to save attributes of svars 
        S=npz_data(); 
        for svar in svars:
            exec('S.{}=self.{}'.format(svar,svar))
        save_npz(fname,S)
        
    def compute_river(self,seg=None,sind=None,acc_limit=1e4,nodata=None,apply_mask=False):   
        #compute river network for watersheds 
        #seg is segment number, sind is catchment indices. 
        
        #compute river mouths        
        if sind is not None:
            sind0=sind.copy(); seg0=self.seg.ravel()[sind0]

        if seg is not None: #may rewrite sind
            if hasattr(seg,'__len__'):
                seg=array(seg)
            else:
                seg=array([seg,])
                
            sind0=nonzero(self.dir.ravel()==0)[0]; seg0=self.seg.ravel()[sind0]
            seg_in,ia,ib=intersect1d(seg,seg0,return_indices=True)
            seg0=seg0[ib]; sind0=sind0[ib];
        
        if (sind is None) and (seg is None):
            sind0=nonzero(self.dir.ravel()==0)[0]; seg0=self.seg.ravel()[sind0]
                    
        #pre-define varibles, will update in the loop
        if nodata is None: nodata=self.nodata
        sind=sind0; slen=len(sind); 
        num=arange(slen).astype('int'); pind0=None; 
        
        print('computing river network')
        
        #each pts cooresponding to a river network
        self.rivers=[[] for i in arange(len(sind0))]; 
                
        while len(sind)!=0:
            #search the largest stream
            sind_list=self.search_upstream(sind,ireturn=7,acc_limit=acc_limit)
                                    
            #add rivers         
            for i in arange(slen):            
                if pind0 is not None: self.rivers[num[i]].append(pind0[i])
                self.rivers[num[i]].extend(sind_list[i])
                self.rivers[num[i]].append(self.nodata)
            
            #create new pind for search 
            pind=[]; pnum=[]; 
            for i in arange(slen):
                pind.extend(sind_list[i])
                pnum.extend(ones(len(sind_list[i]))*num[i])
            pind=array(pind); pnum=array(pnum).astype('int')
        
            #create new sind
            sind,pnum_ind=self.search_upstream(pind,ireturn=5,acc_limit=acc_limit)
            slen=len(sind); num=pnum[pnum_ind];  pind0=pind[pnum_ind]
                    
        #format
        inum=[]
        for i in arange(len(sind0)):
            river=array(self.rivers[i]).astype('int')            
            #apply mask
            if apply_mask:
                fpm=self.mask.ravel()[river]
                river[~fpm]=self.nodata
                
                #only one nodata in between
                nind=nonzero(river==self.nodata)[0]; dind=diff(nind); 
                fpd=nonzero(dind==1)[0]+1; nind=nind[fpd]
                
                fpn=ones(len(river)).astype('bool'); fpn[nind]=False;  
                river=river[fpn]                
                
            self.rivers[i]=river
            if len(self.rivers[i])>3: inum.append(i)
        self.rivers=array(self.rivers); inum=array(inum)
            
        #exclude rivers with len<3
        self.rivers=self.rivers[inum]
        if not hasattr(self,'info'): self.info=npz_data()
        self.info.rseg=seg0[inum]
        self.info.rind=sind0[inum]
        
    
    def compute_boundary(self,seg=None,sind=None,acc_limit=None,level_max=500):
                        
        #pre-define variables
        ds=self.ds; ym,xm=ds
        
        print('computing watershed boundary')
         
        #compute river mouths        
        if sind is not None:
            sind0=sind.copy(); seg0=self.seg.ravel()[sind0]

        if seg is not None: #may rewrite sind
            if hasattr(seg,'__len__'):
                seg=array(seg)
            else:
                seg=array([seg,])
                
            sind0=nonzero(self.dir.ravel()==0)[0]; seg0=self.seg.ravel()[sind0]
            seg_in,ia,ib=intersect1d(seg,seg0,return_indices=True)
            seg0=seg0[ib]; sind0=sind0[ib];
        
        if (sind is None) and (seg is None):
            sind0=nonzero(self.dir.ravel()==0)[0]; seg0=self.seg.ravel()[sind0]
        
        #compute watershed bnd
        sind=sind0.copy()
        
        #exclude pts with acc<acc_limit
        if acc_limit is not None:
            fpa=self.acc.ravel()[sind]>=acc_limit
            sind=sind[fpa]; seg=seg0[fpa]
        
        #already on the boundary
        iy0,ix0=unravel_index(sind,ds); 
        fpn=~((iy0==(ym-1))|(ix0==(xm-1))|(iy0==0)|(ix0==0))        
        
        #search for bounary elements, upward
        sindn=sind[fpn]; segn=self.seg.ravel()[sindn]
        inum=arange(len(sindn)).astype('int'); offset=ds[1]
        while len(inum)!=0:
            sind_next=sindn[inum]-offset
            seg_next=self.seg.ravel()[sind_next]
            iy,ix=unravel_index(sind_next,ds)
                        
            #reach the domain boundary
            fp1=(seg_next==segn[inum])*((iy==(ym-1))|(ix==(xm-1))|(iy==0)|(ix==0))
            
            #already the boudary
            fp2=seg_next!=segn[inum]
            sind_next[fp2]=sindn[inum][fp2]
            
            #exclude fp1 and fp2
            sindn[inum]=sind_next
            inum=inum[~(fp1|fp2)]        
        
        #westward    
        inum=arange(len(sindn)).astype('int'); offset=1
        while len(inum)!=0:
            sind_next=sindn[inum]-offset
            seg_next=self.seg.ravel()[sind_next]
            iy,ix=unravel_index(sind_next,ds)
                        
            #reach the domain boundary
            fp1=(seg_next==segn[inum])*((iy==(ym-1))|(ix==(xm-1))|(iy==0)|(ix==0))
            
            #already the boudary
            fp2=seg_next!=segn[inum]
            sind_next[fp2]=sindn[inum][fp2]
            
            #exclude fp1 and fp2
            sindn[inum]=sind_next
            inum=inum[~(fp1|fp2)]      
            
        sind[fpn]=sindn
                   
        #search boundary
        self.boundary=self.search_boundary(sind,level_max=level_max)  
        
        #add boundary info        
        if not hasattr(self,'info'): self.info=npz_data()
        self.info.bseg=seg
    
    def compute_watershed(self):        
        #acc_max: compute acc and segments
        if not hasattr(self,'dir'): sys.exit('dir not exist')

        #all the catchments        
        sind0=nonzero(self.dir.ravel()==0)[0];         
        
        #initialize acc        
        self.search_upstream(sind0,ireturn=3,level_max=100)
        
        #reorder segment number
        acc0=self.acc.ravel()[sind0]; ind=flipud(argsort(acc0)); sind=sind0[ind]; 
        seg=arange(len(sind)).astype('int')+1
        self.search_upstream(sind,ireturn=3,seg=seg,level_max=100) 
        
        return
                                
    def compute_dir(self,data='dem',outname='dir',subdomain_size=1e5,zlimit=0,method=0):
        #when diff(dem)<zlimit, regard them the same elevation
        #method=0: don't assign dir to flat cells; method=1: assign dir to flat cells
        
        #dem data
        if type(data)==str: #assume data is dem
            self.dem=self.dem.ravel()
            dem0=self.dem
        else:
            dem0=data.ravel()
            
        #pre-calculation                
        ds=self.ds; ym,xm=ds; nsize=prod(ds); nlim=subdomain_size
        dir=zeros(nsize).astype('int32'); nodata=self.nodata
        nsubdomain=int(max(round(nsize/nlim),1))
        
        #convert nan to nodata
        if isnan(nodata):
            nodata=-99999
            fpn=isnan(dem0)
            dem0[fpn]=nodata
            
        offsets_all=array([1,-xm,-1,xm, -xm+1,-xm-1,xm-1,xm+1])
        ndirs_all=array([1, 64, 16, 4, 128, 32, 8,  2])
        pdirs_all=array([16, 4,  1, 64, 8,   2, 128,32])
           
        #calculate dir for each subdomain      
        for it in arange(3):
            if it==0: 
                nsub=nsubdomain #for subdomain
            else:
                nsub=4 #for 4 sides and 4 corners   

            #for each section   
            for i in arange(nsub):
                if it==0:
                    #get index for subdomain
                    if i==(nsub-1):
                        sind0=arange(i*nlim,nsize).astype('int')
                    else:
                        sind0=arange(i*nlim,(i+1)*nlim).astype('int')
                    flag=arange(8)
                elif it==1:
                    #get index for each side
                    if i==0:
                        sind0=arange(1,xm-1);                       
                        flag=array([0,2,3,6,7 ])
                    elif i==1:
                        sind0=arange(1,xm-1)+(ym-1)*xm; 
                        flag=array([0,1,2,4,5])
                    elif i==2:
                        sind0=arange(1,ym-1)*xm
                        flag=array([0,1,3,4,7])
                    elif i==3:
                        sind0=arange(1,ym-1)*xm+(xm-1)                        
                        flag=array([1,2,3,5,6])
                    dir[sind0]=0
                elif it==2:
                    #get index for each corner
                    if i==0:
                        sind0=array([0])
                        flag=array([0,3,7])
                    elif i==1:
                        sind0=array([xm-1])                        
                        flag=array([2,3,6])
                    elif i==2:
                        sind0=array([(ym-1)*xm])
                        flag=array([0,1,4])
                    elif i==3:
                        sind0=array([nsize-1])                        
                        flag=array([1,2,5])
                    dir[sind0]=0
                                                
                #begin compute dir-----------------------------------------------
                #define offsets for each sid-e
                offsets=offsets_all[flag]
                ndirs=ndirs_all[flag]
                #pdirs=pdirs_all[flag]
                               
                #exclude nodata
                sdem=dem0[sind0]; 
                fp=sdem!=nodata; dir[sind0[~fp]]=-1                
                sind=sind0[fp]; dem=sdem[fp];  slen=sum(fp)
                if len(sind)==0: continue
                    
                #find the neighbor with min. depth
                ndem=ones(slen)*nodata; nind=zeros(slen).astype('int');  ndir=nind.copy() # pdir=nind.copy()
                for m in arange(len(offsets)):
                    sindi=sind+offsets[m]; 
                    fpo=sindi>=nsize; sindi[fpo]=nsize-1
                    demi=dem0[sindi]
                    
                    #check values
                    fpm1=(demi!=nodata)*(ndem==nodata)
                    fpm2=(demi!=nodata)*(ndem!=nodata)*(demi<ndem)
                    fpm=fpm1|fpm2                    
                        
                    #assign new ndem if ndem==nodata and demi!=nodata, and if demi<ndem            
                    nind[fpm]=sindi[fpm]
                    ndem[fpm]=demi[fpm]
                    ndir[fpm]=ndirs[m]
                    #pdir[fpm]=pdirs[m]
                                     
                #when dem>ndem
                fpm=(ndem!=nodata)*(dem>(ndem+zlimit))                
                dir[sind[fpm]]=ndir[fpm]
                
                #when dem=ndem, assign dir first
                if method==1:
                    fpm=(ndem!=nodata)*(abs(dem-ndem)<=zlimit)*(sind<nind)
                    dir[sind[fpm]]=ndir[fpm]
        
        #convert nodata to nan
        if isnan(self.nodata): dem0[fpn]=self.nodata
        if type(data)==str: self.dem=self.dem.reshape(ds) #assume data is dem
                                                        
        #reshape
        dir=dir.reshape(ds)
        exec('self.{}=dir'.format(outname))        
        
    def fill_depression(self):
        
        #change nan nodata to -99999
        S.remove_nan_nodata();
        
        #resolve cells with dir==0, but with no upstream cells
        self.compute_dir(method=1);
        sind0=nonzero(self.dir.ravel()==0)[0]
        isum=self.search_upstream(sind0,ireturn=2)
        self.search_flat(sind0[isum==0],ireturn=6)
           
    def resolve_flat(self,sind0=None,zlimit=0):
        #resolve flat based on following paper
        #An Efficient Assignment of Drainage Direction Over Flat Surfaces In 
        #Raster Digital Elevation Models". Computers & Geosciences. doi:10.1016/j.cageo.2013.01.009"
        
        #pre-define variables
        ds=self.ds; ym,xm=ds; nsize=ym*xm
        
        #find indices with (dir==0)
        if sind0 is None:
            sind0=nonzero(self.dir.ravel()==0)[0]       
        
        #return how many neighboring indices with same elevation
        isum=self.search_flat(sind0,ireturn=1,zlimit=zlimit) #indices of flat
        fp=isum>0; sind0=sind0[fp]
        
        #obtain flat cell indices
        sind_flat=self.search_flat(sind0,ireturn=3,zlimit=zlimit); #sind_flat=unique(sind_flat)
        dir_flat0=self.dir.ravel()[sind_flat];
        
        #find the low edges
        fpl=dir_flat0!=0; sind_low=sind_flat[fpl]
        
        #find the high edges
        fp=dir_flat0==0; sind_high=self.search_flat(sind_flat[fp],ireturn=2,zlimit=zlimit); 
                
        #initialize new dem, and assing dem begining from high edges        
        dem_high=self.search_flat(sind_high,ireturn=4,zlimit=zlimit)[sind_flat]        
                         
        #initialize new dem, and assing dem begining from low edges        
        dem_low=self.search_flat(sind_low,ireturn=5,zlimit=zlimit)[sind_flat]
                        
        #combine dem_high and dem_low 
        #dem_new=zeros(nsize).astype('int')*1e5
        dem_new=2*dem_low+dem_high
        
        #calcuate dir for sind_flat
        dir_flat=self.search_flat(sind_flat,ireturn=9,zlimit=zlimit,dem_flat=dem_new)
        dir_flat[fpl]=dir_flat0[fpl]; self.dir.ravel()[sind_flat]=dir_flat
        
        #combine dem_high and dem_low 
        #dem_new=ones(nsize)*self.nodata; 
        #dem_new[sind_flat]=2*dem_low+dem_high
        #replace original dem, this method fails
        #self.compute_dir(data=dem_new.reshape(ds),outname='dir_tmp',zlimit=zlimit)             
        #sind_new=setdiff1d(sind_flat,sind_low)        
        #self.dir.ravel()[sind_new]=self.dir_tmp.ravel()[sind_new]                
        #delattr(self,'dir_tmp')
        
    def search_flat(self,sind0,ireturn=0,zlimit=0,wlevel=0,level=0,level_max=500,dem_flat=None):
        #ireturn=0: return one-level nearby indices with same elevation
        #ireturn=1: return how many neighboring indices with same elevation        
        #ireturn=2: return high edges of flat
        #ireturn=3: return all the nearby indice with same elevation (all the cells in flat)   
        #ireturn=4: calculate elevation from high edges
        #ireturn=5: calculate elevation from low edges
        #ireturn=6: assign dir to the dir=0 neighbor that is without upstream cells
        #ireturn=7: return neighboring indices in seg==0 & with minimum depth.
        #ireturn=8: return neighboring indices in neighboring seg & with minimum depth &dir.
        #ireturn=9: assign dir for sind_flat based on dem_flat
        
        #pre-define variables
        slen=len(sind0); ds=self.ds; ym,xm=ds; nsize=ym*xm
        dem0=self.dem.ravel()[sind0]
        
        #convert index
        iy0,ix0=unravel_index(sind0,ds)
        
        #construct maxtrix for neighbor
        #yind=r_[iy0-1,iy0-1,iy0-1,iy0,iy0+1,iy0+1,iy0+1,iy0]
        #xind=r_[ix0+1,ix0,ix0-1,ix0-1,ix0-1,ix0,ix0+1,ix0+1]        
        yind=r_[iy0, iy0-1,iy0,iy0+1, iy0-1,iy0-1, iy0+1, iy0+1]
        xind=r_[ix0+1,ix0,ix0-1,ix0,  ix0+1,ix0-1, ix0-1, ix0+1]
        
        #true neighbors
        fpt=nonzero((xind>=0)*(xind<xm)*(yind>=0)*(yind<ym))[0]
        sind_true=ravel_multi_index([yind[fpt],xind[fpt]],ds); dem_true=self.dem.ravel()[sind_true]
        
        if ireturn==6:
            #get neighboring minimum dem and indices
            fp=dem_true==self.nodata;  dem_true[fp]=1e5
            dem=ones([8,slen])*1e5; dem.ravel()[fpt]=dem_true            
            iy_min=argmin(dem,axis=0); dem_min=dem[iy_min,arange(slen)]            
            
            #modify sind0's dir and dem
            dir0=array([128,64,32,16,8,4,2,1])
            self.dir.ravel()[sind0]=dir0[iy_min]
            self.dem.ravel()[sind0]=dem_min
            
        if ireturn==7:
            #get neighboring minimum dem and indices in seg=0
            fp=(self.seg.ravel()[sind_true]!=0)|(dem_true==self.nodata); dem_true[fp]=1e5;
            dem=ones([8,slen])*1e5; dem.ravel()[fpt]=dem_true
            
            sind=zeros([8,slen]).astype('int'); sind.ravel()[fpt]=sind_true
            iy_min=argmin(dem,axis=0);             
            sind_min=sind[iy_min,arange(slen)]; dem_min=dem[iy_min,arange(slen)];             
            return sind_min, dem_min
        
        if ireturn==8:
            #get neighboring minimum dem and indices in seg=0
            #dir=tile(array([128,64,32,16,8,4,2,1]),slen).reshape([slen,8]).T
            #dir=tile(array([8,4,2,1,128,64,32,16]),slen).reshape([slen,8]).T            
            dir=tile(array([16,4,1,64,8,2,128,32]),slen).reshape([slen,8]).T            
            seg_true=self.seg.ravel()[sind_true]; seg0=self.seg.ravel()[sind0]
            fp=(seg_true==tile(seg0,8)[fpt])|(dem_true==self.nodata); dem_true[fp]=1e5;
            dem=ones([8,slen])*1e5; dem.ravel()[fpt]=dem_true
            
            sind=zeros([8,slen]).astype('int'); sind.ravel()[fpt]=sind_true
            iy_min=argmin(dem,axis=0);             
            sind_min=sind[iy_min,arange(slen)]; 
            dem_min=dem[iy_min,arange(slen)];
            dir_min=dir[iy_min,arange(slen)]; 
            return sind_min, dem_min,dir_min
        
        if ireturn==2:
            fph=nonzero(dem_true>(tile(dem0,8)[fpt]+zlimit))[0]
            num=zeros(slen*8); num[fpt[fph]]=1;             
            isum=num.reshape([8,slen]).sum(axis=0); 
            fp=isum>0; sind_next=sind0[fp]
            return sind_next
        
        #neighbor with same elevation
        fps=nonzero(abs(dem_true-tile(dem0,8)[fpt])<=zlimit)[0]
        sind_next=sind_true[fps]; 
                
        #exclude low edges
        if ireturn==4:            
            fpl=nonzero(self.dir.ravel()[sind_next]==0)[0]; fpl_len=len(sind_next)
            sind_next=sind_next[fpl]
                    
        #unique indices
        sind_next,fpu,fpu_inverse=unique(sind_next,return_index=True,return_inverse=True)
        
        if ireturn==9:
            dir_flat=zeros(slen).astype('int')        
            dem0=int(1e5)*ones(nsize).astype('int'); dem0[sind0]=dem_flat;
            dem=int(1e5)*ones([8,slen]).astype('int'); dem.ravel()[fpt[fps]]=dem0[sind_true[fps]]
            iy_min=argmin(dem,axis=0); dem_min=dem[iy_min,arange(slen)]
            dir_min=tile(array([1,64,16,4,128,32,8,2]),[slen,1]).T[iy_min,arange(slen)]           
            fp=dem_flat>dem_min; dir_flat[fp]=dir_min[fp]
            return dir_flat
        
        if ireturn==0:
            return sind_next
        
        if ireturn==1:
            num=zeros(slen*8); num[fpt[fps]]=1
            isum=reshape(num,[8,slen]).sum(axis=0).astype('int')
            return isum
        
        if ireturn in [3,4,5]:
            if wlevel==0: #first level loop
                #init
                self.sind_next=sind0.copy()
                self.sind_list=[] 
                self.flag_search=True
                self.dem_flag=zeros(nsize).astype('int32')
                                
                if ireturn in [4,5]:
                    self.v0=1
                    self.vmax=None                    
                    self.dem_flat=zeros(nsize).astype('int32')                    
                if ireturn==4: self.sind_list.append(sind0)
                                 
                #1st round of search loop from outside to inside
                while self.flag_search:
                    sind_next=self.sind_next
                    if len(sind_next)!=0:
                        self.search_flat(sind_next,ireturn=ireturn,wlevel=1,level=0,level_max=level_max)
                        
                #2nd round of search loop from inside to outside
                if ireturn==4: 
                    self.dem_flag[:]=0
                    for i in arange(len(self.sind_list)):
                        self.vmax=self.search_flat(self.sind_list[-i-1],ireturn=ireturn,wlevel=2,level=0,level_max=level_max)
                
                #save results
                if ireturn==3: 
                    sind=array(self.sind_list)
                else:
                    dem_flat=self.dem_flat
                    
                #clean
                delattr(self,'sind_next'); delattr(self,'sind_list');delattr(self,'dem_flag'); delattr(self,'flag_search')
                if ireturn in [4,5]: delattr(self,'v0'); delattr(self,'vmax');delattr(self,'dem_flat')
                
                #return results
                if ireturn==3: 
                    return sind
                else:
                    return dem_flat

            elif wlevel==1: #2nd level search, for ireturn==3
                self.dem_flag[sind0]=1
                fpn=self.dem_flag[sind_next]==0; sind_next=sind_next[fpn]
                if ireturn==3: self.sind_list.extend(sind0);  
                x=ravel_multi_index([1437,2574],ds) in sind0
                
                if ireturn in [4,5]: self.dem_flat[sind0]=(ones(slen)*(level+self.v0)).astype('int32')
                
                if level!=level_max:                    
                    if sum(fpn)!=0: #continue                      
                        self.search_flat(sind_next,ireturn=ireturn,wlevel=1,level=level+1,level_max=level_max)
                    else:
                        self.flag_search=False #reach the end                        
                else:                
                    if sum(fpn)!=0:
                        self.sind_next=sind_next                         
                        if ireturn in [4,5]:
                            self.sind_list.append(sind_next)
                            self.v0=int(self.v0+level_max)                        
                                              
            elif wlevel==2: #2nd level search, for ireturn==4 
                self.dem_flag[sind0]=1
                fpn=nonzero(self.dem_flag[sind_next]==0)[0]; fpn_len=len(sind_next)
                sind_next=sind_next[fpn]                    
                vmax=self.dem_flat[sind0].copy()
                
                if level!=level_max:                                       
                    if sum(fpn)!=0: #continue                      
                        vmax_next=self.search_flat(sind_next,ireturn=ireturn,wlevel=2,level=level+1,level_max=level_max)                                                                                          
                    else:
                        self.flag_search=False #reach the end
                else:                
                    if sum(fpn)!=0:
                        vmax_next=self.vmax
                    
                #modify dem_flat
                if sum(fpn)!=0:
                    #replace vmax with vmax_next if vmax_next is larger                    
                    num_n=zeros(fpn_len); num_n[fpn]=vmax_next                                        
                    num=zeros(8*slen);  num[fpt[fps[fpl]]]=num_n[fpu_inverse] 
                    
                    nmax=reshape(num,[8,slen]).max(axis=0)                                       
                    fp=vmax<nmax; vmax[fp]=nmax[fp]
                    
                self.dem_flat[sind0]=vmax-self.dem_flat[sind0]                    
                
                return vmax
        
    def search_boundary(self,sind0,wlevel=0,level=0,level_max=500):
        #search pts of boundary pts of watershed segment, sind0 are the initial index of segments
        
        #pre-defind variables
        slen=len(sind0); ds=self.ds
        # print(sind0)
        iy0,ix0=unravel_index(sind0,ds)
        seg0=self.seg.ravel()[sind0]
        
        #construct matrix
        yind=r_[iy0-1,iy0-1,iy0-1,iy0,iy0+1,iy0+1,iy0+1,iy0]
        xind=r_[ix0+1,ix0,ix0-1,ix0-1,ix0-1,ix0,ix0+1,ix0+1]
        sind=zeros(8*slen).astype('int'); seg=sind.copy()
                
        #true neighbor
        fpt=nonzero((xind>=0)*(xind<ds[1])*(yind>=0)*(yind<ds[0]))[0]
        sind_true=ravel_multi_index([yind[fpt],xind[fpt]],ds); 
        sind[fpt]=sind_true; seg[fpt]=self.seg.ravel()[sind_true]
        
        #indices belong to the same segment
        fps=seg!=tile(seg0,8); sind[fps]=0
        
        #exclude previous cells
        # if hasattr(self,'pind'): fpn=sind==tile(self.pind,8); sind[fpn]=0
        if hasattr(self,'ppind'): fpn=sind==tile(self.ppind,8); sind[fpn]=0
        
        sind=sind.reshape([8,slen])
        
        iy0=zeros(slen).astype('int'); 
        if hasattr(self,'pind'):
            pind=self.pind
        else:
            pind=0
        #find the starting search pt
        for i in arange(8):
            fp=sind[i,:]==pind
            iy0[fp]=i
        
        #find the next pts
        sind_next=zeros(slen).astype('int')
        ix=arange(slen).astype('int')
        for i in arange(8):
            iy0=iy0+1
            iy1=mod(iy0,8); iy2=mod(iy0+1,8)
            fpn=(sind[iy1,ix]==0)*(sind[iy2,ix]!=0)*(sind_next==0)
            sind_next[fpn]=sind[iy2,ix][fpn]
        
        #if sind_next=0, search back
        fpb=sind_next==0; sind_next[fpb]=sind0[fpb]      
              
        #recursive search 
        if wlevel==0:
            #init
            self.sind_next=sind0.copy()
            self.sind0=sind0.copy()
            self.pind=sind0.copy() #previous pts
            self.ppind=sind0.copy() #previous pts
            self.flag_search=True
            self.sind_list=[[i] for i in sind0]
            self.snum=arange(slen).astype('int')
            
            #1-level search loop
            iflag=0
            while self.flag_search:                
                iflag=iflag+1;
                if len(self.sind_next)==0: break
                self.search_boundary(self.sind_next,wlevel=1,level=0,level_max=level_max)
                #print('search bnd: {}'.format(iflag*level_max))
                
            #save list info
            sind_list=array([array(i) for i in self.sind_list])
            
            #clean
            delattr(self,'sind_next'); delattr(self,'sind0'); delattr(self,'flag_search')
            delattr(self,'sind_list'); delattr(self,'snum'); delattr(self,'pind');delattr(self,'ppind');
            
            return sind_list
            
        elif wlevel==1:
            fpz=sind_next!=self.sind0
            [self.sind_list[i].append(j) for i,j in zip(self.snum[~fpz],sind_next[~fpz])]
            if sum(fpz)==0: self.flag_search=False
            if (level!=level_max) and sum(fpz)!=0: #not reach the end                                
                self.sind0=self.sind0[fpz]
                self.snum=self.snum[fpz]
                self.ppind=self.pind[fpz]
                self.pind=self.sind_next[fpz]
                self.sind_next=sind_next[fpz]
                
                #save pts to the list    
                [self.sind_list[i].append(j) for i,j in zip(self.snum,self.sind_next)]
                
                #continue search
                self.search_boundary(self.sind_next,wlevel=1,level=level+1,level_max=level_max)
                                     
    def search_upstream(self,sind0,ireturn=0,seg=None,acc_limit=0,wlevel=0,level=0,level_max=500,acc_calc=False):
        #ireturn=0: all the 1-level upstream index. if seg is not None, return seg number (seg_up) also
        #ireturn=1: all the 1-level upstream index, also with flags of true neighbor and incoming flow
        #ireturn=2: number of upstream indices 
        #ireturn=3: search all upstreams, compute acc and seg
        #ireturn=4: just one-level upstream index with largest acc 
        #ireturn=5: all the one-level upstream index except the one with largest acc, and numbering index
        #ireturn=6: uppermost upstream index along the largest river stem
        #ireturn=7: save all the cells along the largest upstream river
        #ireturn=8: return how many neighbors with same segment number after acc and seg are known
        #ireturn=9: return segment boundary indices
             
        #--pre-define variables        
        slen=len(sind0); ds=self.ds
        dir_in0=[8,4,2,1,128,64,32,16]
        
        #convert index
        iy0,ix0=unravel_index(sind0,ds)
                       
        #compute init neighbor index
        yind=r_[iy0-1,iy0-1,iy0-1,iy0,iy0+1,iy0+1,iy0+1,iy0]
        xind=r_[ix0+1,ix0,ix0-1,ix0-1,ix0-1,ix0,ix0+1,ix0+1]
        if seg is not None: seg0=tile(seg,8)
        
        #true neighbor index
        fpt=nonzero(((xind>=0)*(xind<ds[1])*(yind>=0)*(yind<ds[0])).ravel())[0]
        
        #init dir
        dir_in=tile(dir_in0,[len(sind0),1]).T.ravel()
        
        #true neighbor
        sind_true=ravel_multi_index([yind[fpt],xind[fpt]],ds);
        if seg is not None: seg_true=seg0[fpt]
        
        #dir diff for true neighbors  
        dir_diff=self.dir.ravel()[sind_true]-dir_in[fpt]
        
        #index of neighbors with incoming flow
        fpf=nonzero(dir_diff==0)[0]
        sind_up=sind_true[fpf] 
        if seg is not None: seg_up=seg_true[fpf]
        
        if ireturn==0: 
            if seg is None:
                return sind_up
            else:
                return sind_up, seg_up
        
        if ireturn==1:
            return sind_up, fpt[fpf]
        
        if ireturn==2:
            num=zeros(slen*8); num[fpt[fpf]]=1
            num_up=(reshape(num,[8,slen]).sum(axis=0)).astype('int')
            return num_up
        
        if ireturn==8:
            seg0=self.seg.ravel()[sind0]
            #method 0
            #fps=tile(seg0,8)[fpt]!=self.seg.ravel()[sind_true]
            #num=zeros(slen*8); fnum=ones(len(fpt)); fnum[fps]=0; num[fpt]=fnum
            #method 1
            fps=nonzero(tile(seg0,8)[fpt]==self.seg.ravel()[sind_true])[0]
            num=zeros(slen*8); num[fpt[fps]]=1
            
            num_up=(reshape(num,[8,slen]).sum(axis=0)).astype('int')
            return num_up
        
        if ireturn==3:
            if wlevel==0: #first level recursive search 
                
                #init
                if seg is None: 
                    seg0=arange(slen)+1
                    self.acc=zeros(prod(ds)).astype('int')                                                  
                else:
                    seg0=seg
                    
                self.seg=zeros(prod(ds)).astype('int')                
                self.sind_list=[]; 
                self.seg_list=[];
                self.flag_search=True
                                    
                self.sind_next=sind0.copy(); self.seg_next=seg0; 
                self.sind_list.append(sind0); self.seg_list.append(seg0) 
                      
                #1-level search loop, from downstream to upstream
                iflag=0
                while self.flag_search:
                    iflag=iflag+1
                    print('search upstream: {}, {}'.format(iflag,len(self.sind_next)))
                    if len(self.sind_next)==0: break
                    self.search_upstream(self.sind_next,ireturn=ireturn,seg=self.seg_next, wlevel=1,level=0,level_max=level_max)
                
                #search from upstream to downstream on the 1-level
                for i in arange(len(self.sind_list)):   
                    if seg is not None: continue
                    print('search downstream: {}, {}'.format(i,len(self.sind_list[-i-1])))
                    self.search_upstream(self.sind_list[-i-1],ireturn=ireturn,seg=self.seg_list[-i-1], wlevel=1,level=0,level_max=level_max,acc_calc=True)
                
                #clean 
                delattr(self,'sind_list');delattr(self,'seg_list'); delattr(self,'flag_search')
                delattr(self,'sind_next'); delattr(self,'seg_next')
                if seg is None: delattr(self,'seg')
                    
                #reshape
                if hasattr(self,'acc'): self.acc=self.acc.reshape(ds)                
                if seg is not None: self.seg=self.seg.reshape(ds)
                
                return
                         
            elif wlevel==1: #2-level recursive search    
                #assign seg number
                self.seg[sind0]=seg            
               
                acc_up=0; acc=0
                if level!=level_max:    
                    if len(sind_up)==0: #reach the end
                        self.flag_search=False                        
                    else:                                             
                        acc_up=self.search_upstream(sind_up,ireturn=ireturn,seg=seg_up,wlevel=1,level=level+1,level_max=level_max,acc_calc=acc_calc) 
                                         
                else:
                    if not acc_calc:
                        self.sind_next=sind_up
                        self.seg_next=seg_up
                        self.sind_list.append(sind_up)                    
                        self.seg_list.append(seg_up)
                    else:
                        acc_up=self.acc[sind_up]  
                
                if acc_calc:
                    #compute acc for sind0
                    #get flags of sind_up first      
                    sind_up,fptf=self.search_upstream(sind0,ireturn=1)
                    acc0=zeros(slen*8); acc0[fptf]=acc_up                     
                    acc=reshape(acc0,[8,slen]).sum(axis=0)+1
                    
                    self.acc[sind0]=acc
                
                return acc   
        
        #----------------------------------------------------------------------
        #code below works after self.acc is computed, except for ireturn=9
        #----------------------------------------------------------------------
        sind=-ones(8*slen).astype('int'); sind[fpt[fpf]]=sind_up
        
        if ireturn==9:                                                            
            sind_next=sind_up           
        else:
            #get corresponding acc
            acc=zeros(slen*8).astype('int'); acc[fpt[fpf]]=self.acc.ravel()[sind_up];            
            acc=acc.reshape([8,slen]); sind=sind.reshape([8,slen]) 
            
            #apply limit
            if acc_limit!=0:
                fpc=acc<acc_limit; acc[fpc]=0; sind[fpc]=-1;
                
            #get index for cells_next with maximum acc
            ix=arange(slen).astype('int'); iy=argmax(acc,axis=0)        
            sind_next=sind[iy,ix]
           
        #just one level
        if ireturn==4:
            return sind_next
        
        #one-level upstream, all the cells except the largest one
        if ireturn==5:            
            sind[iy,ix]=-1; fpn=sind!=-1;
            nind0=tile(arange(slen),8).reshape([8,slen]).astype('int') 
            sind_next=sind[fpn]; nind=nind0[fpn]
            return sind_next,nind
        
        #get to the uppermost pts, using 2-level recursive search. 
        #1-level recursive search will crash reaching the setrecursionlimit
        if ireturn in [6,7,9]:            
            if wlevel==0: #first level recursive
                #init
                #if not hasattr(self,'sind_next'):
                self.sind_next=sind0.copy()
                self.sind_flag=ones(slen)
                self.flag_search=True
                self.sind_list=[[i] for i in sind0]
                self.sind_bnd=[[] for i in sind0]
                self.sind_seg=[[] for i in sind0]
                self.seg_next=arange(slen)                                   
                
                #1-level search loop
                while self.flag_search:
                    if ireturn==9:
                        sind_next=self.sind_next
                    else:
                        fp=self.sind_flag>0; sind_next=self.sind_next[fp]                        
                    if len(sind_next)==0: break;                    
                    self.search_upstream(sind_next,ireturn=ireturn,seg=self.seg_next,acc_limit=acc_limit,wlevel=1,level=0,level_max=level_max,acc_calc=acc_calc) 
                
                #search finished at this point
                sind_next=self.sind_next                
                sind_list=array([array(i) for i in self.sind_list]) 
                sind_bnd=array([array(i) for i in self.sind_bnd]) 
                sind_seg=array([sort(array(i)) for i in self.sind_seg]) 
                
                #clean
                delattr(self,'sind_next'); delattr(self,'flag_search'); delattr(self,'sind_list')
                delattr(self,'sind_bnd'); delattr(self,'seg_next'); delattr(self,'sind_flag'); delattr(self,'sind_seg')
                
                if ireturn==6:
                    return sind_next
                elif ireturn==7:
                    return sind_list
                elif ireturn==9:
                    if acc_calc:
                        return sind_seg
                    else:
                        return sind_bnd
                                     
            elif wlevel==1: #second level recursive search                 
                fpnz=nonzero(sind_next!=-1)[0]; fpz=nonzero(sind_next==-1)[0]
                if level!=level_max:
                    if ireturn==9:
                        seg_next=seg_up
                        #collect boundary cells                              
                        isum=sind.reshape([8,slen]).sum(axis=0); fpb=isum==-8
                        [self.sind_bnd[i].append(j) for i,j in zip(seg[fpb],sind0[fpb])]
                        if acc_calc:
                            [self.sind_seg[i].append(j) for i,j in zip(seg,sind0)]
                        
                    if sum(fpnz)==0:
                        #reach the end
                        if ireturn!=9: fpn=self.sind_flag>0; self.sind_next[fpn]=sind0
                        if ireturn==7:
                            ind_list=nonzero(fpn)[0]; sind_list=sind0
                            [self.sind_list[i].append(j) for i,j in zip(ind_list,sind_list)]
                        self.flag_search=False
                    else:
                        if ireturn==9:
                            self.sind_next=sind_up
                            self.seg_next=seg_next
                        else:
                            #save the pts that reach the end
                            sind_next[fpz]=sind0[fpz]
                            fpn=self.sind_flag>0; self.sind_next[fpn]=sind_next; self.sind_flag[nonzero(fpn)[0][fpz]]=0
                        if ireturn==7:
                            ind_list=nonzero(fpn)[0]; sind_list=abs(sind_next)
                            [self.sind_list[i].append(j) for i,j in zip(ind_list,sind_list)]
                        
                        #continue search the rest pts 
                        self.search_upstream(sind_next[fpnz],ireturn=ireturn,seg=seg_next,acc_limit=acc_limit,wlevel=1,level=level+1,level_max=level_max,acc_calc=acc_calc)
                        
    def search_downstream(self,sind0,ireturn=0,wlevel=0,level=0,level_max=500):
        #ireturn=0: return 1-level downstream index
        #ireturn=1: downmost downstream index
        #ireturn=2: save all the index along the downstream
        
        #pre-define variables
        slen=len(sind0); ds=self.ds
        
        #variables for each direction
        dir_out=array([128,64,32,16,8,4,2,1]).astype('int')
        offset_y=array([-1,-1,-1,0,1,1,1,0]).astype('int')
        offset_x=array([1,0,-1,-1,-1,0,1,1]).astype('int')
        
        #check each dir
        sdir=self.dir.ravel()[sind0]
        iy0,ix0=unravel_index(sind0,ds)       
        yind=-ones(slen).astype('int'); xind=yind.copy(); sind_next=yind.copy()
        for i in arange(8):
            fp=sdir==dir_out[i]
            yind[fp]=iy0[fp]+offset_y[i]
            xind[fp]=ix0[fp]+offset_x[i]
        
        #true donwstream
        fpt=(xind>=0)*(xind<ds[1])*(yind>=0)*(yind<ds[0])
        if sum(fpt)!=0:
            sind_next[fpt]=ravel_multi_index([yind[fpt],xind[fpt]],ds)
            
        if ireturn==0:
            return sind_next
        
        if ireturn==1 or ireturn==2:
            if wlevel==0: #first level recrusive search
                #init
                #if not hasattr(self,'sind_next'):
                self.sind_next=sind0.copy()
                self.flag_search=True
                self.sind_list=[[i] for i in sind0]
                
                #1-level search loop
                while self.flag_search:
                     fp=self.sind_next>0; sind_next=self.sind_next[fp]
                     if sum(fp)==0: break;                    
                     self.search_downstream(sind_next,ireturn=ireturn,wlevel=1,level=0,level_max=level_max) 
                      
                #search finished at this point
                sind_next=-self.sind_next                
                sind_list=array([array(i) for i in self.sind_list])                
                
                #clean
                delattr(self,'sind_next'); delattr(self,'flag_search'); delattr(self,'sind_list')
                
                if ireturn==1:
                    return sind_next
                elif ireturn==2:
                    return sind_list  
                
            elif wlevel==1: #second level recursive search
                fpz=sind_next!=-1                
                if level!=level_max:
                    if sum(fpz)==0:
                        #reach the end
                        fpn=self.sind_next>0; self.sind_next[fpn]=-sind0
                        if ireturn==2:
                            ind_list=nonzero(fpn)[0]; sind_list=sind0
                            [self.sind_list[i].append(j) for i,j in zip(ind_list,sind_list)]
                        self.flag_search=False
                    else:
                        #save the pts that reach the end
                        sind_next[~fpz]=-sind0[~fpz]
                        fpn=self.sind_next>0; self.sind_next[fpn]=sind_next
                        if ireturn==2:
                            ind_list=nonzero(fpn)[0]; sind_list=abs(sind_next)
                            [self.sind_list[i].append(j) for i,j in zip(ind_list,sind_list)]
                        
                        #continue search the rest pts 
                        self.search_downstream(sind_next[fpz],ireturn=ireturn,wlevel=1,level=level+1,level_max=level_max)                              
            
    def compute_extent(self):
        #recalculate extent from affine        
        ym,xm=self.ds; dxy=self.affine[0]; ly0=self.affine[5]; lx0=self.affine[2]        
        self.extent=array([lx0,lx0+(xm-1)*dxy,ly0-(ym-1)*dxy,ly0])
        
    def remove_nan_nodata(self,nodata=-99999):
        fpn=isnan(self.dem)
        self.dem[fpn]=nodata
        self.nodata=nodata
    
    def get_coordinates(self,sind0,affine=None,ds=None,nodata=None):     
        #example: sind0=nonzero(self.dir.ravel()==0)[0], or sind0=nonzero(self.dir==0)
        
        sind=sind0.copy()
        #check parameters first
        if ds is None: ds=self.ds
        if affine is None: affine=self.affine
        if nodata is None: nodata=self.nodata
        
        #check index
        if len(sind)==2 and hasattr(sind[0],'__len__'):            
            indy,indx=sind  #sind0=[indy,indx]            
            fpn=(indy==nodata)|(indx==nodata)
        else:
            fpn=sind0==nodata; sind[fpn]=0;
            indy,indx=unravel_index(sind,ds)
            indy[fpn]=nodata; indx[fpn]=nodata
            
        #convert index to coordinates
        dxy=affine[0]; ly0=affine[5]; lx0=affine[2]
        
        xi=indx*dxy+lx0
        yi=-indy*dxy+ly0
        
        xi[fpn]=nodata; yi[fpn]=nodata
        
        return yi,xi
    
    def write_shapefile(self,data,sname,stype='POLYLINE',prjname='epsg:4326',attname=None,attvalue=None):
        #data: string name or data
        #sname: shapefile name
        
        #get data
        S=npz_data()
        if type(data)==type(''):
            exec('S.data=self.{}'.format(data))
        else:
            S.data=data
        
        #get type and prj
        S.type=stype
        S.prj=get_prj_file(prjname)
        
        #get attrs
        if (attname is not None) and (attvalue is not None):
            S.attname=attname
            S.attvalue=attvalue
        
        #get xy
        S.xy=[];
        for i in arange(len(S.data)):            
            yi,xi=self.get_coordinates(S.data[i])
            fp=(xi==self.nodata)|(yi==self.nodata);
            xi[fp]=nan; yi[fp]=nan;            
            S.xy.append(c_[xi,yi])
        S.xy=array(S.xy)
        
        #write shapefile
        write_shapefile_data(sname,S)
                       
         
if __name__=="__main__":    
    close('all')
    sys.setrecursionlimit(100000)

    import copy
    
    #S0=dem_dir();S0.read_data('S1.npz') 
    
    S=dem_dir(); S.read_data('S2_DEM.npz'); S.ds=S.dem.shape; ds=S.ds; ym,xm=ds
    dzm=min(diff(unique(S.dem)))*0.9; S.remove_nan_nodata();
        
    S.compute_dir(method=0);
    
    #resolve flats 
    sind0=nonzero(S.dir.ravel()==0)[0]; 
    isum=S.search_upstream(sind0,ireturn=2)
    S.resolve_flat(sind0[isum==0])    

    #resolve flats for single cells
    sind0=nonzero(S.dir.ravel()==0)[0]; 
    isum=S.search_upstream(sind0,ireturn=2)    
    S.search_flat(sind0[isum==0],ireturn=6)
      
    #identify depression catchment
    sind0=nonzero(S.dir.ravel()==0)[0]; iy,ix=unravel_index(sind0,ds)
    fp=(iy>0)*(iy<(ym-1))*(ix>0)*(ix<(xm-1)); sind0=sind0[fp]
    
    #search and mark each depression
    slen=len(sind0); seg0=arange(slen)+1; h0=S.dem.ravel()[sind0];  
    S.search_upstream(sind0,ireturn=3,seg=seg0,level_max=100)
    
    #get indices for each depression; here seg is numbering, not segment number
    sind_segs=S.search_upstream(sind0,ireturn=9,seg=arange(slen),acc_calc=True)
    ns0=array([len(i) for i in sind_segs]) 
    
    #get depression boundary
    S.compute_boundary(sind=sind0); 
    
    #loop to reverse dir along streams
    ids=arange(len(sind0))
    while len(sind0)!=0:
        A0=S.dir.copy(); B0=S.seg.copy();
        #-------------------------------------------------------------------------
        #seg0,     ns0,     h0,      sind0,       sind_bnd,   h_bnd
        #seg_min,  ns_min,  h0_min,  sind0_min,   sind_min,   h_min  h_bnd_min  dir_min
        #-------------------------------------------------------------------------        
        # #search and mark each depression
        # slen=len(sind0); seg0=arange(slen)+1; h0=S.dem.ravel()[sind0];  
        # S.search_upstream(sind0,ireturn=3,seg=seg0,level_max=100)
        
        # #get indices for each depression; here seg is numbering, not segment number
        # sind_segs=S.search_upstream(sind0,ireturn=9,seg=arange(slen),acc_calc=True)
        # ns0=array([len(i) for i in sind_segs])    
            
        # #get depression boundary
        # S.compute_boundary(sind=sind0); 
        
        #rearange the boundary index
        sind_bnd_all=[]; ids_bnd=[]; id1=0; id2=0;  
        for i in arange(len(sind0)):            
            sindi=unique(S.boundary[i]); id2=id1+len(sindi)
            sind_bnd_all.extend(sindi)        
            ids_bnd.append(arange(id1,id2).astype('int')); id1=id1+len(sindi)             
        sind_bnd_all=array(sind_bnd_all);ids_bnd=array(ids_bnd);
        
        #find all the neighboring indices with minimum depth
        sind_min_all,h_min_all,dir_min_all=S.search_flat(sind_bnd_all,ireturn=8)
                
        #minimum for each seg
        h_min=array([h_min_all[id].min() for id in ids_bnd])
        mind=array([nonzero(h_min_all[id]==hi)[0][0] for id,hi in zip(ids_bnd,h_min)])        
        sind_bnd=array([sind_bnd_all[id][mi] for id,mi in zip(ids_bnd,mind)]); h_bnd=S.dem.ravel()[sind_bnd]
        sind_min=array([sind_min_all[id][mi] for id,mi in zip(ids_bnd,mind)])
        dir_min=array([dir_min_all[id][mi] for id,mi in zip(ids_bnd,mind)])
        
        #get s0_min,h0_min,sind0_min
        seg_min=S.seg.ravel()[sind_min]; fps=nonzero(seg_min!=0)[0]; h_bnd_min=h_bnd[seg_min-1]
        ns_min=zeros(slen).astype('int'); h0_min=-1e5*ones(slen); sind0_min=-ones(slen).astype('int') 
        ns_min[fps]=ns0[seg_min[fps]-1];  h0_min[fps]=h0[seg_min[fps]-1]; sind0_min[fps]=sind0[seg_min[fps]-1];  
    
        #get stream from head to catchment
        fp1=seg_min==0; 
        fp2=(seg_min!=0)*(h0_min<h0); 
        fp3=(seg_min!=0)*(h0_min==h0)*(ns_min>ns0)
        fp4=(seg_min!=0)*(h0_min==h0)*(ns_min==ns0)*(h_bnd>h_bnd_min); 
        fp5=(seg_min!=0)*(h0_min==h0)*(ns_min==ns0)*(h_bnd==h_bnd_min)*(sind0>sind0_min);         
        fph=nonzero(fp1|fp2|fp3|fp4|fp5)[0]; sind_head=sind_bnd[fph]; 
        sind_streams=S.search_downstream(sind_head,ireturn=2)
        
        #get all stream index and its dir
        t0=time.time(); sind=[]; dir=[];
        for i in arange(len(fph)):
            id=fph[i];    
            #S.seg.ravel()[sind_segs[id]]=seg_min[id]; seg0[id]=seg_min[id]
            sind_stream=sind_streams[i][:-1]
            dir_stream=r_[dir_min[id],S.dir.ravel()[sind_stream][:-1]]            
            sind.extend(sind_stream); dir.extend(dir_stream)
        sind=array(sind); dir0=array(dir).copy(); dir=zeros(len(dir0)).astype('int')
            
        #reverse dir
        dir_0=[128,64,32,16,8,4,2,1]; dir_inv=[8,4,2,1,128,64,32,16]
        for i in arange(8):
            fpr=dir0==dir_0[i]; dir[fpr]=dir_inv[i]
        dt=time.time()-t0
        S.dir.ravel()[sind]=dir;
                
        #reset variables---------
        sind0=setdiff1d(sind0,sind0[fph]); delattr(S,'seg')
        #search and mark each depression
        slen=len(sind0); seg0=arange(slen)+1; h0=S.dem.ravel()[sind0];  
        S.search_upstream(sind0,ireturn=3,seg=seg0,level_max=100)
        
        #get indices for each depression; here seg is numbering, not segment number
        sind_segs=S.search_upstream(sind0,ireturn=9,seg=arange(slen),acc_calc=True)
        ns0=array([len(i) for i in sind_segs]) 
        
        #get depression boundary
        S.compute_boundary(sind=sind0); 
        
#------------------------------------------------------------------------------
    sys.exit()
    
    while True:
        #init dir: resolve cells with dir==0, but with no upstream cells
        
        S.compute_dir(method=1);
        sind0=nonzero(S.dir.ravel()==0)[0]
        isum=S.search_upstream(sind0,ireturn=2)
        S.search_flat(sind0[isum==0],ireturn=6)

        #identify depression catchment
        sind0=nonzero(S.dir.ravel()==0)[0]; iy,ix=unravel_index(sind0,ds)
        fp=(iy>0)*(iy<(ym-1))*(ix>0)*(ix<(xm-1)); sind0=sind0[fp]
        
        print(len(sind0)); dir=S.dir.copy(); dem=S.dem.copy()
        if len(sind0)==0: break
        #-------------------------------------------------------------------------
        #seg0,     ns0,     h0,      sind0,       sind_bnd,   h_bnd
        #seg_min,  ns_min,  h0_min,  sind0_min,   sind_min,   h_min  h_bnd_min
        #-------------------------------------------------------------------------
        #search and mark each depression
        slen=len(sind0); seg0=arange(slen)+1; h0=S.dem.ravel()[sind0];  
        S.search_upstream(sind0,ireturn=3,seg=seg0,level_max=100)
        
        #get indices for each depression; here seg is numbering, not segment number
        sind_segs=S.search_upstream(sind0,ireturn=9,seg=arange(slen),acc_calc=True)
        ns0=array([len(i) for i in sind_segs])    
            
        #get depression boundary
        S.compute_boundary(sind=sind0); 
        
        #rearange the boundary index
        sind_bnd_all=[]; ids_bnd=[]; id1=0; id2=0;  
        for i in arange(slen):
            sindi=unique(S.boundary[i]); id2=id1+len(sindi)
            sind_bnd_all.extend(sindi)        
            ids_bnd.append(arange(id1,id2).astype('int')); id1=id1+len(sindi)             
        sind_bnd_all=array(sind_bnd_all);ids_bnd=array(ids_bnd);
        
        #find all the neighboring indices with minimum depth
        sind_min_all,h_min_all=S.search_flat(sind_bnd_all,ireturn=8)
        
        #minimum for each seg
        h_min=array([h_min_all[id].min() for id in ids_bnd])
        mind=array([nonzero(h_min_all[id]==hi)[0][0] for id,hi in zip(ids_bnd,h_min)])        
        sind_bnd=array([sind_bnd_all[id][mi] for id,mi in zip(ids_bnd,mind)]); h_bnd=S.dem.ravel()[sind_bnd]
        sind_min=array([sind_min_all[id][mi] for id,mi in zip(ids_bnd,mind)])
        
        #get s0_min,h0_min,sind0_min
        seg_min=S.seg.ravel()[sind_min]; fps=nonzero(seg_min!=0)[0]; h_bnd_min=h_bnd[seg_min-1]
        ns_min=zeros(slen).astype('int'); h0_min=-1e5*ones(slen); sind0_min=-ones(slen).astype('int') 
        ns_min[fps]=ns0[seg_min[fps]-1];  h0_min[fps]=h0[seg_min[fps]-1]; sind0_min[fps]=sind0[seg_min[fps]-1];  
    
        #get stream
        fp1=seg_min==0; 
        fp2=(seg_min!=0)*(h0_min<h0); 
        fp3=(seg_min!=0)*(h0_min==h0)*(ns_min>ns0)
        fp4=(seg_min!=0)*(h0_min==h0)*(ns_min==ns0)*(h_bnd>h_bnd_min); 
        fp5=(seg_min!=0)*(h0_min==h0)*(ns_min==ns0)*(h_bnd==h_bnd_min)*(sind0>sind0_min);         
        fph=nonzero(fp1|fp2|fp3|fp4|fp5)[0]; sind_head=sind_bnd[fph]; sind_min_head=sind_min[fph]
        sind_streams=S.search_downstream(sind_head,ireturn=2)
        
        t0=time.time();  sind_tmp=[];
        for i in arange(len(fph)):
            id=fph[i]; 
            z0=h0[id]; z1=h_min[id]
            sind_seg=sind_segs[id]; sind_tmp.extend(sind_seg)
            sind_stream=sind_streams[i]
            h_seg=S.dem.ravel()[sind_seg]; fpz=nonzero(h_seg<=z1)[0]
            if z1<z0: sys.exit('impossible')
            if min(h_seg[fpz])<z0: sys.exit('impossible 2')
            if z0!=z1: S.dem.ravel()[sind_seg[fpz]]=z1+dzm*(h_seg[fpz]-z0)/(z1-z0)        
            S.dem.ravel()[sind_stream]=z1
        dt=time.time()-t0
                    
        sys.exit()
        #compute dir again
        dir_head=S.dir.ravel()[sind_min_head].copy(); fp=dir_head==0; S.dir.ravel()[sind_min_head[fp]]=-1;
        S.resolve_flat(sind_min_head)
        
        sys.exit();
    
     
    sys.exit()
    
    iy0,ix0=unravel_index(sind0,ds)
    #--compute again
    S.compute_dir(method=0);
    sind0=nonzero(S.dir.ravel()==0)[0]
    isum=S.search_upstream(sind0,ireturn=2)
    S.search_flat(sind0[isum==0],ireturn=6)  
    S.resolve_flat()
    
    #identify depression catchment
    sind0=nonzero(S.dir.ravel()==0)[0]; iy,ix=unravel_index(sind0,ds)
    fp=(iy>0)*(iy<(ym-1))*(ix>0)*(ix<(xm-1)); sind0=sind0[fp]
    
    slen=len(sind0); seg0=arange(slen)+1; h0=S.dem.ravel()[sind0];  
    S.search_upstream(sind0,ireturn=3,seg=seg0,level_max=100)
    
    # sys.exit()
    
    #plot
    sz=S.seg.astype('float'); sz[sz==0]=nan; imshow(sz);
    iy1,ix1=unravel_index(sind0,ds)
    for i in arange(len(S.boundary)):            
        iy,ix=unravel_index(S.boundary[i],ds)
        plot(ix,iy,'k-');
    plot(ix0,iy0,'ro',ix1,iy1,'b.')
    
    sys.exit()
        
    while len(sind0)>0:  
        slen=len(sind0); seg0=arange(slen)+1; 
        print('ZG:', slen)
        

        
        
        # #plot
        # figure()
        # imshow(sx);
        # iy0,ix0=unravel_index(sind0,ds)
        # iy2,ix2=unravel_index(sind_min[dem_min!=1e5],ds)
        # for i in arange(len(S.boundary)):            
        #     iy,ix=unravel_index(S.boundary[i],ds)
        #     plot(ix,iy,'k-');
        # plot(ix0,iy0,'r.',ix2,iy2,'b.')
        
        # sys.exit(0)
        
        #modify depth in each depression
        ids_head=[]; sind_head=[]; dem_head=[]; iseg=[]; id1=0; id2=0;  hmins=[];
        for i in arange(slen):
            id=ids_bnd[i]
            sind_bndi=sind_bnd[id];dem_bndi=dem_bnd[id]; sind_mini=sind_min[id]; dem_mini=dem_min[id]; 
            
            #sort
            indsort=argsort(dem_mini); 
            sind_bndi=sind_bndi[indsort]; sind_mini=sind_mini[indsort]
            dem_bndi=dem_bndi[indsort]; dem_mini=dem_mini[indsort];         
            
            #get minimum depths
            h1=dem_mini[0]; 
            if h1==100000: continue 
        
            #check results 
            #if h1<dem0[i]:
            #    sys.exit('impossible: neighboring minimum is lower than catchment')
        
            #get the 2nd minimum depths
            demi_unique=unique(dem_mini)        
            if len(demi_unique)>1: h2=demi_unique[1]                      
            if h2==100000: h2=h1
                    
            #save all the indices with minimum depth neighbors
            fpm=nonzero(dem_mini==h1)[0]
            id2=id1+len(fpm); ids_head.append(arange(id1,id2)); id1=id1+len(fpm)        
            sind_head.extend(sind_bndi[fpm])        
            dem_head.extend(dem_bndi[fpm]); 
            hmins.append([dem0[i],h1,h2])        
            iseg.append(i)
            
        ids_head=array(ids_head); sind_head=array(sind_head); dem_head=array(dem_head); 
        iseg=array(iseg); hmins=array(hmins)
             
        #search for downstream sind_link
        sind_stream=S.search_downstream(sind_head,ireturn=2)
        len_stream=array([len(i) for i in sind_stream])
        

        # #plot
        # figure()
        # imshow(sx);
        # iy0,ix0=unravel_index(sind0,ds)
        # iy1,ix1=unravel_index(sind_head,ds)        
        # for i in arange(len(S.boundary)):            
        #     iy,ix=unravel_index(S.boundary[i],ds)
        #     plot(ix,iy,'k-');
        # plot(ix0,iy0,'r.',ix1,iy1,'b.')
        sx=S.seg.astype('float'); sx[sx==0]=nan;
                        
        sind_stream_true=[]
        #modify depth in each depression
        for i in arange(len(ids_head)):
            id=ids_head[i]
            len_streami=len_stream[id]; sind_streami=sind_stream[id]
            
            sind_headi=sind_head[id]; dem_headi=dem_head[id];  
            h0,h1,h2=hmins[i]; isegi=iseg[i]
            
            #find the out channel of depression
            fps=nonzero(len_streami==len_streami.min())[0] #shortest route
            fph=nonzero(dem_headi[fps]==dem_headi[fps].min())[0] #lowest head
            idi=id[fps[fph]][0]
         
            # #modify depression depth        
            sind_segi=sind_seg[isegi]; dem_segi=S.dem.ravel()[sind_segi]
            
            if h1==h2:
                fph=nonzero(dem_segi<h1)[0]; S.dem.ravel()[sind_segi[fph]]=h1
            else:
                fph=nonzero((dem_segi>h0)*(dem_segi<h2))[0];             
                S.dem.ravel()[sind_segi[fph]]=h1+(h2-h1)*(dem_segi[fph]-h0)/(h2-h0)    
            
            S.dem.ravel()[sind_stream[idi]]=h1
            sind_stream_true.append(sind_stream[idi])
            
            #set seg num to zero
            S.seg.ravel()[sind_segi]=0
        
        #reset and continue work on other segs
        iseg_left=setdiff1d(arange(slen),iseg)
        
        #plot
        figure()
        imshow(sx);
        iy0,ix0=unravel_index(sind0,ds)
        iy1,ix1=unravel_index(sind_head,ds)        
        for i in arange(len(S.boundary)):            
            iy,ix=unravel_index(S.boundary[i],ds)
            plot(ix,iy,'k-');
        plot(ix1,iy1,'b.'); plot(ix0,iy0,'r.')
        for i in arange(len(sind_stream_true)):
            iy2,ix2=unravel_index(array(sind_stream_true[i]),ds)
            plot(ix2,iy2,'c-')        
        
        sind0=sind0[iseg_left]
        sind_seg=sind_seg[iseg_left]
        S.boundary=S.boundary[iseg_left]
        dem0=dem0[iseg_left]
    
        
    # A=S.dem.copy(); A2=A-A0;  
    # B=S.dir; C=S.seg;
    # #make sure every cell(dir=0) has upstream
    #resolve cells with dir==0, but with no upstream cells
    # S.compute_dir(method=0);
    # S.resolve_flat()
    # S.compute_watershed()
    
    # sind0=nonzero(S.dir.ravel()==0)[0]
    # isum=S.search_upstream(sind0,ireturn=2)
    
    
    sys.exit()
    
    #plot
    # seg=S.seg.astype('float'); fp=seg==0; seg[fp]=nan;
    # figure();
    # imshow(seg)
    
    # for i in arange(len(S.boundary)):
    #     iy,ix=unravel_index(S.boundary[i],ds)
    #     plot(ix,iy,'r.')
    
        
    
    # # # #compute dir

    
    # # iy,ix=unravel_index(sind,ds)
    # # fp=(iy>0)*(iy<(ym-1))*(ix>0)*(ix<(xm-1))
    
    # # siy=iy[fp]; six=ix[fp]; ssind=sind[fp]
    # # biy=iy[~fp]; bix=ix[~fp]; bsind=sind[~fp]
    
    # # S.search_flat(bsind,ireturn=6)
    
    # dt1=time.time()-t0
    
    # # S.resolve_flat()
    # S.compute_dir()
    # S.compute_watershed()
    # dt2=time.time()-t0
    
    # A=S.dem;
    # B=S.dir;
    
    # # sind0=nonzero(S.dir.ravel()==0)[0]
    # # isum=S.search_upstream(sind0,ireturn=2)
    
    # # sys.exit()
    
    # # # #fill depression
    # # # sind0=nonzero(S.dir.ravel()==0)[0]; 
    # # # isum=S.search_upstream(sind0,ireturn=2); sind=sind0[isum==8]
    # # # S.fill_depression(sind0,nodata=S.nodata)
    
    
    # # # #resolve flats
    # # # S.compute_dir('dem');
    # # # S.resolve_flat(zlimit=1e-5);
    
    # # # sys.exit()
    
    # # # #compute watershed
    # # # S.compute_dir('dem');
    # # # S.compute_watershed()
    
    # # # S1.compute_dir('dem');
    # # # S1.resolve_flat();
    # # # S1.compute_watershed()
    
    # figure();
    # subplot(2,1,1)
    # imshow(S0.acc,vmin=0,vmax=1e3)
    
    # subplot(2,1,2)
    # imshow(S.acc,vmin=0,vmax=1e3)
    
    # sys.exit()
#-----------------test resolve flats-------------------------------------------  
    # S=dem_dir(); #S.read_data('S2_DEM.npz')
    
    # ds=[20,20]
    # S.dem=arange(prod(ds)).reshape(ds); S.ds=ds; S.nodata=-9999;
    # iy,ix=nonzero(S.dem>-99)
    
    # fp=sqrt((iy-6)**2+(ix-6)**2)<5; S.dem.ravel()[fp]=-2; S.dem[14,13]=-8;
    
    # S.dem[13:18,13:18]=-2; S.dem[18,14]=-8;  A=S.dem
    
    # S.compute_dir(data='dem');
    # S.resolve_flat(); B=S.dir    

#------write shapefile---------------------------------------------------------
    # S=dem_dir(); S.read_data('S2.npz'); 
    # acc_limit=1e2; seg=arange(1,1e5)
    # S.compute_river(seg,acc_limit=acc_limit,apply_mask=True)
    # S.write_shapefile('rivers','B_rivers')
    # # S.compute_boundary(seg,acc_limit=acc_limit)            
    # # S.write_shapefile('boundary','B_boundary')    

#--------------------------precalculation--------------------------------------
    # # # # # # # # #--------------------------------------------------------------
    # # # compute dir and acc, then save them
    # # # fname='./13arcs/southern_louisiana_13_navd88_2010.asc'; sname='S3'
    # # # fname='ne_atl_crm_v1.asc'; sname='S2'
    # fname='GEBCO.asc'; sname='S1'
        
    # #read dem info
    # gds=read_deminfo(fname,'dem',size_subdomain=1e7); 
    # gds.add_deminfo(fname,'dem0',size_subdomain=1e9); 
    
    # #read all data first
    # gds.read_dem_data('dem0')  
    
    # flat_option=0

    # S=dem_dir();
    # #reconstrcut data from subdomains
    # id=-1; print(gds.dem_info.nsubdomain)
    # for i in arange(gds.dem_info.ny):
    #     for j in arange(gds.dem_info.nx):            
    #         id=id+1;  print(id)
            
    #         #read dem,            
    #         gds.read_dem_data('dem',[id],dem_data=gds.dem0_data[0].data);
    #         # gds.read_dem_data('dem',[id]);
    #         gd=gds.dem_data[0]

    #         #processing raw data
    #         if flat_option==0:
    #             try:                
    #                 gd.fill_depressions(gd.data, out_name='fdem')
    #                 try:
    #                     gd.resolve_flats(gd.fdem,'flats')
    #                 except:
    #                     gd.flats=gd.fdem
    #             except:
    #                 try:
    #                     gd.resolve_flats(gd.data,'flats')
    #                 except:
    #                     gd.flats=gd.data
    #         elif flat_option==1:
    #             #don't use resolve flats
    #             try:                
    #                 gd.fill_depressions(gd.data, out_name='fdem')    
    #                 gd.flats=gd.fdem
    #             except:      
    #                 gd.flats=gd.data
            
    #         #local data
    #         ixy=gds.dem_info.domains[id].rind            
    #         demii=array(gd.flats[ixy[0]:ixy[1],ixy[2]:ixy[3]])
            
    #         if j==0:
    #             demi=demii
    #         else:
    #             demi=c_[demi,demii]
    #     #vertical
    #     if i==0:
    #         dem=demi
    #     else:
    #         dem=r_[dem,demi]
                
    #     #vertical
    #     if i==0:
    #         S.dem=demi
    #     else:
    #         S.dem=r_[S.dem,demi]
    
    # #save information
    # S.nodata=gd.nodata; S.ds=S.dem.shape; S.affine=[*gds.dem_info.affine]
    # S.compute_extent();
        
    # #compute dir
    # t0=time.time();
    # S.compute_dir(S.dem,nodata=S.nodata);
    # S.compute_watershed()
    # dt=time.time()-t0
    
    # # S.save_data(sname,['dir','acc','seg','nodata','extent','ds','affine'])
    # #add mask
    # S.mask=(S.dem>=0)*(S.dem<=50)
    # S.save_data(sname,['dir','acc','seg','nodata','extent','ds','affine','mask'])
    
#------------------------plot rivers and seg boundary--------------------------
    # S=dem_dir(); S.read_data('S1.npz'); 
    
    # acc_limit=1e3; seg=arange(1,10)
    # S.compute_river(seg,acc_limit=acc_limit)
    # S.compute_boundary(seg,acc_limit=acc_limit)
    

    # figure();        
    # #plot acc
    # imshow(S.acc,vmin=0,vmax=5e3,extent=S.extent)
    
    # #plot rivers
    # for i in arange(len(S.rivers)):
    #     sindi=S.rivers[i].copy();
    #     yi,xi=S.get_coordinates(sindi)
        
    #     fp=(xi==S.nodata)|(yi==S.nodata);
    #     xi[fp]=nan; yi[fp]=nan;
        
    #     plot(xi,yi,'r-')

    # #plot river mouths
    # yi,xi=S.get_coordinates(S.info.rind)
    # plot(xi,yi,'b*',ms=12)
    
    # #plot seg boundary
    # colors='www'
    # for i in arange(len(S.boundary)):
    #     sindi=S.boundary[i];
    #     # sindi=array(S.sind_list[i])
    #     yi,xi=S.get_coordinates(array(sindi))
        
    #     fp=(xi==S.nodata)|(yi==S.nodata);
    #     xi[fp]=nan; yi[fp]=nan;
        
    #     plot(xi,yi,color=colors[mod(i,3)],lw=1)    
    
    
#---------------test method "search_downstream---------------------------------
    # S=dem_dir(); S.read_data('S1.npz'); S.compute_extent();
    
    # sind0=ravel_multi_index(nonzero((S.seg<=10)*(S.seg>0)*(S.dir==0)),S.ds);   
    
    # sind_up=S.search_upstream(sind0,ireturn=6)
    # list_down=S.search_downstream(sind_up,ireturn=2)
    # list_up=S.search_upstream(sind0,ireturn=7)
    
    # dyi,dxi=S.get_coordinates(sind0)
    # uyi,uxi=S.get_coordinates(sind_up)
    
    # figure();     
    # subplot(2,1,1)
    # imshow(S.acc,vmin=0,vmax=1e2,extent=S.extent)
    # plot(dxi,dyi,'b*',uxi,uyi,'w*')
    # for i in arange(len(sind0)):
    #     sindi=list_down[i];
    #     yi,xi=S.get_coordinates(sindi)
        
    #     fp=(xi==S.nodata)|(yi==S.nodata);
    #     xi[fp]=nan; yi[fp]=nan;
        
    #     plot(xi,yi,'r-')        
    
        
    # subplot(2,1,2)
    # imshow(S.acc,vmin=0,vmax=1e2,extent=S.extent)    
    # plot(dxi,dyi,'b*',uxi,uyi,'w*')
    # for i in arange(len(sind0)):
    #     sindi=list_up[i];
    #     yi,xi=S.get_coordinates(sindi)
        
    #     fp=(xi==S.nodata)|(yi==S.nodata);
    #     xi[fp]=nan; yi[fp]=nan;
        
    #     plot(xi,yi,'r-')    
    