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
        
    def extract_river(self,sind0,acc_limit=1e4,nodata=None):
                        
        #pre-define varibles, will update in the loop
        if nodata is None: nodata=self.nodata
        sind=sind0; slen=len(sind); 
        num=arange(slen).astype('int'); pind0=None; 
        
        #each pts cooresponding to a river network
        self.rivers=[[] for i in arange(len(sind0))]; 
                
        while len(sind)!=0:
            #search the largest stream
            sind_list=self.search_upstream(sind,ireturn=6,acc_limit=acc_limit)
                                    
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
            sind,pnum_ind=self.search_upstream(pind,ireturn=4,acc_limit=acc_limit)
            slen=len(sind); num=pnum[pnum_ind];  pind0=pind[pnum_ind]
                    
        #format            
        for i in arange(len(sind0)):
            S.rivers[i]=array(S.rivers[i]).astype('int')        
        
    def compute_acc(self):        
        #acc_max: compute acc and segments
        if not hasattr(self,'dir'): sys.exit('dir not exist')
        
        #initialize acc
        sind=nonzero(self.dir.ravel()==0)[0]; 
        seg=arange(len(sind)).astype('int')+1
        self.search_watershed(sind,seg=seg);  
        
        #reorder segment number
        delattr(self,'seg')
        acc=self.acc[sind]; ind=flipud(argsort(acc)); sind=sind[ind]; 
        self.search_watershed(sind,seg=seg,segsort=True); 
        
        #reshape
        self.acc=self.acc.reshape(self.ds)
        self.seg=self.seg.reshape(self.ds)
        
        return
                                
    def compute_dir(self,data='dem',nodata=-9999,subdomain_size=1e6): 
        #dem data
        if type(data)==str: 
            exect('self._temp=self.{}'.format(data));
            dem=self._temp; delattr(self,'_temp')
        else:
            dem=data
            
        #pre-calculation                
        ds=dem.shape; ym,xm=ds; nsize=prod(ds); nlim=subdomain_size
        self.dir=zeros(nsize).astype('int32')
        nsubdomain=round(nsize/nlim)
        
        offsets_all=array([1,-xm+1,-xm,-xm-1,-1,xm-1,xm,xm+1])
        ndirs_all=array([1,128,64,32,16,8,4,2])
        pdirs_all=array([16,8,4,2,1,128,64,32])
           
        #calculate dir for each subdomain
        dem=dem.ravel();
        
        for it in arange(4):
            if it==0: 
                nsub=nsubdomain #for subdomain
            elif it==1 or it==2:
                nsub=4 #for 4 sides and 4 corners
            elif it==3:
                nsub=6 #fill_depression
                
                #for cells with dir==0 but no upstream cells
                sind00=nonzero(self.dir==0)[0]
                sind_up=self.search_upstream(sind00,ireturn=2)
                fp=sind_up==0; sind00=sind00[fp]; 
                yind,xind=unravel_index(sind00,ds)
        
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
                        flag=array([0,4,5,6,7])
                    elif i==1:
                        sind0=arange(1,xm-1)+(ym-1)*xm; 
                        flag=array([0,1,2,3,4])
                    elif i==2:
                        sind0=arange(1,ym-1)*xm
                        flag=array([0,1,2,6,7])
                    elif i==3:
                        sind0=arange(1,ym-1)*xm+(xm-2)
                        flag=array([2,3,4,5,6])
                    self.dir[sind0]=0
                elif it==2:
                    #get index for each corner
                    if i==0:
                        sind0=array([0])
                        flag=array([0,6,7])
                    elif i==1:
                        sind0=array([xm-1])
                        flag=array([4,5,6])
                    elif i==2:
                        sind0=array([(ym-1)*xm])
                        flag=array([0,1,2])
                    elif i==3:
                        sind0=array([nsize-1])
                        flag=array([2,3,4])
                    self.dir[sind0]=0
                elif it==3:                    
                        if i==0: #for cells inside
                            fp=(xind>0)*(xind<(xm-1))*(yind>0)*(yind<(ym-1))
                            sind0=sind00[fp]                                                        
                            flag=arange(8)
                        elif i==1: #for side
                            fp=(xind>0)*(xind<(xm-1))*(yind==0)
                            sind0=sind00[fp] 
                            flag=array([0,4,5,6,7])
                        elif i==2:
                            fp=(xind>0)*(xind<(xm-1))*(yind==(ym-1))
                            sind0=sind00[fp] 
                            flag=array([0,1,2,3,4])
                        elif i==3:
                            fp=(yind>0)*(yind<(ym-1))*(xind==0)
                            sind0=sind00[fp] 
                            flag=array([0,1,2,6,7])
                        elif i==4:
                            fp=(yind>0)*(yind<(ym-1))*(xind==(xm-1))
                            sind0=sind00[fp]
                            flag=array([2,3,4,5,6])
                        elif i==5: #for 4 corners
                            sind0=array([0,xm-1,(ym-1)*xm,nsize-1])
                            dir0=array([2,8,128,32])
                            for ii in arange(4):
                                if sind0[ii] in sind00: self.dir[sind0[ii]]=dir0[ii] 
                            continue
                        
                        if len(sind0)==0: continue
                                                
                #begin compute dir-----------------------------------------------
                #define offsets for each sid-e
                offsets=offsets_all[flag]
                ndirs=ndirs_all[flag]
                pdirs=pdirs_all[flag]
                               
                #exclude nodata
                fp=array(dem[sind0]!=nodata);  self.dir[sind0[~fp]]=-1
                sind0=sind0[fp]; dem0=dem[sind0]; 
                if len(sind0)==0: continue
                    
                #find the neighbor with min. depth
                nind=sind0.copy(); ndem=ones(len(sind0))*nodata; ndir=zeros(len(sind0)); pdir=zeros(len(sind0)) 
                for m in arange(len(offsets)):
                    sindi=sind0+offsets[m]; 
                    fp=sindi>=nsize; sindi[fp]=nsize-1
                    demi=dem[sindi]
                    
                    #check values
                    fp1=(demi!=nodata)*(ndem==nodata)
                    fp2=(demi!=nodata)*(ndem!=nodata)*(demi<ndem)            
                    fp=fp1|fp2
                    
                    #assign new ndem if ndem==nodata and demi!=nodata, and if demi<ndem            
                    nind[fp]=sindi[fp]
                    ndem[fp]=demi[fp]
                    ndir[fp]=ndirs[m]
                    pdir[fp]=pdirs[m]
                                
                if it==3:
                    fp=(ndem!=nodata); self.dir[sind0[fp]]=ndir[fp]
                else:
                    zlim=0
                    #compare with dem0
                    fp1=(ndem!=nodata)*(dem0>ndem)*(dem0>(ndem+zlim)); 
                    fp2=(ndem!=nodata)*(abs(dem0-ndem)<=zlim)*(nind>sind0)
                    fp=fp1|fp2; self.dir[sind0[fp]]=ndir[fp]
                    
                    fp2=(ndem!=nodata)*(abs(dem0-ndem)<=zlim)*(nind<sind0)
                    # fp=fp2; self.dir[sind0[fp]]=pdir[fp]
                    fp=fp2; self.dir[sind0[fp]]=0
        
        #reshape
        self.dir=self.dir.reshape(ds)
        
    def search_watershed(self,sind0,seg=1,level_max=100,wlevel=0,segsort=False):
        #sind0: index of catchment pts. 
        #example: sind0=nonzero(self.dir.ravel()==0)[0], or sind0=nonzero(self.dir==0)
        #algorithm: two layer of recursive, to avoid errors, and improve efficiency
        
        #assign variables
        ds=self.dir.shape   
        if wlevel==0:
            if not hasattr(self,'seg'): self.seg=zeros(self.dir.size).astype('int')
            if not hasattr(self,'acc'): self.acc=zeros(self.dir.size).astype('int')
            self.sinds=[]; self.sinds.append([ds,seg,level_max])
            #self.segs=[]; 
        
            #convert index
            if len(sind0)==2 and hasattr(sind0[0],'__len__'):
                sind=ravel_multi_index(sind0,ds)  #sind0=[indy,indx]
            else:
                sind=sind0
                
            if not hasattr(sind,'__len__'): sind=array([sind])
            
            #format segs
            if hasattr(seg,'__len__') and len(sind)!=len(seg): sys.exit('sind0 and seg: not same length')
            if hasattr(seg,'__len__'):
                segs=seg
            else:
                segs=ones(len(sind))*seg     
                
            #exclude watershed that was already claimed
            fp=nonzero(self.seg[sind]==0)[0]
            sind=sind[fp]; segs=segs[fp]
        
        #define sind and segs
        if wlevel!=0: 
            sind=sind0; segs=seg
        
        if len(sind)==0: return
        
        print('search upstream: {}, {}'.format(wlevel,len(sind)))
        #save sind and segs    
        self.sinds.append(sind)
        
        #searach watershed from downstream to upstream
        self.get_neighbor_cells(sind,seg=segs,level=0,level_max=level_max)        
        
        #figure();imshow(C.seg.reshape(ds))
        #continue search if level_max is reached
        sind_next=nonzero(self.seg<0)[0]; seg_next=-self.seg[sind_next]
        if len(sind_next)!=0:
            self.search_watershed(sind_next,seg=seg_next,level_max=level_max,wlevel=wlevel+1)
            
        #search from upstream to downstream 
        if wlevel==0: 
            ds,seg,level_max=self.sinds[0]
            ns=len(self.sinds)-1
            for i in arange(ns):          
                if segsort: continue
                print('search downstream: {}, {}'.format(i,len(self.sinds[-i-1])))
                self.get_neighbor_cells(self.sinds[-i-1],seg=0,level=0,level_max=level_max,acc_calc=True)        
            delattr(self,'sinds')
    
    def get_neighbor_cells(self,sind0,seg=0,level=0,level_max=100,acc_calc=False):
        #determine watershed area first
        if acc_calc==False:
            if level!=level_max:
                #get index_next                      
                ind_next,seg_next=self.search_upstream(sind0,seg=seg)
                            
                #search index_next
                self.get_neighbor_cells(ind_next,seg=seg_next,level=level+1,level_max=level_max,acc_calc=acc_calc)
                
                #assign seg number
                self.seg[sind0]=seg
            else:
                #assign seg number
                self.seg[sind0]=-seg
         
        #compuate accumlation
        elif acc_calc==True:        
            if level!=level_max: 
                fpn=self.acc[sind0]==0; sind00=sind0[fpn]
                
                #get index_next                     
                ind_next,fpt,fpf=self.search_upstream(sind00,ireturn=1)
                
                #get acc of index_next
                if len(ind_next)!=0: 
                    acc_next=self.get_neighbor_cells(ind_next,seg=None,level=level+1,level_max=level_max,acc_calc=acc_calc)                    
                else:
                    acc_next=0                
                
                #compute and assign acc in sind00
                cc=zeros(len(sind00)*8); fcc=cc[fpt]; fcc[fpf]=acc_next; cc[fpt]=fcc
                acc0=(reshape(cc,[8,len(sind00)]).sum(axis=0))+1                
                self.acc[sind00]=acc0
             
            acc=self.acc[sind0]
                     
            return acc    

    def search_watershed_bnd(self,sind0=None):
        #search pts of boundary pts of watershed segment, sind0 are the initial index of segments
        
        #check index
        if sind0 is None:
            sind=sind0
        else:
            sind=nonzero(self.dir.ravel()==0)[0]
            
        #check closed segments
        sind_next=self.search_upstream(sind,self.ds,ireturn=2)
        
        #search boundary index 
        self.search_segment_bnd(sind)

    def search_upstream(self,sind0,ireturn=0,seg=None,acc_limit=0,wlevel=0,level=0,level_max=500):
        #ireturn=0: all the 1-level upstream index. if seg is not None, return seg number (seg_up) also
        #ireturn=1: all the 1-level upstream index, also with flags of true neighbor and incoming flow
        #ireturn=2: number of upstream indices        
        #ireturn=3: just one-level upstream index with largest acc 
        #ireturn=4: all the one-level upstream index except the one with largest acc, and numbering index
        #ireturn=5: uppermost upstream index along the largest river stem
        #ireturn=6: save all the cells along the largest upstream river
             
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
            return sind_up,fpt,fpf
        
        if ireturn==2:
            num=zeros(slen*8); fnum=num[fpt]; fnum[fpf]=1; num[fpt]=fnum; 
            num_up=(reshape(num,[8,slen]).sum(axis=0)).astype('int')
            return num_up
        
        #----------------------------------------------------------------------
        #code below works after self.acc is computed
        #----------------------------------------------------------------------
        #get corresponding acc
        acc=zeros(slen*8).astype('int'); sind=acc.copy()
        facc=acc[fpt]; facc[fpf]=self.acc.ravel()[sind_up]; acc[fpt]=facc; 
        fsind=sind[fpt]; fsind[fpf]=sind_up; sind[fpt]=fsind
        acc=acc.reshape([8,slen]); sind=sind.reshape([8,slen]) 
        
        #apply limit
        if acc_limit!=0:
            fpc=acc<acc_limit; acc[fpc]=0; sind[fpc]=0;
            
        #get index for cells_next with maximum acc
        ix=arange(slen).astype('int'); iy=argmax(acc,axis=0)        
        sind_next=sind[iy,ix]
       
        #just one level
        if ireturn==3:
            return sind_next
        
        #one-level upstream, all the cells except the largest one
        if ireturn==4:            
            sind[iy,ix]=0; fpn=sind!=0;
            nind0=tile(arange(slen),8).reshape([8,slen]).astype('int') 
            sind_next=sind[fpn]; nind=nind0[fpn]
            return sind_next,nind
        
        #get to the uppermost pts, using 2-level recursive search. 
        #1-level recursive search will crash reaching the setrecursionlimit
        if ireturn==5 or ireturn==6:            
            if wlevel==0: #first level recursive
                #init
                if not hasattr(self,'sind_next'):
                    self.sind_next=sind0.copy()
                    self.flag_search=True
                    self.sind_list=[[i] for i in sind0]                    
                
                #1-level search loop
                while self.flag_search:
                    fp=self.sind_next>0; sind_next=self.sind_next[fp]
                    if sum(fp)==0: break;                    
                    self.search_upstream(sind_next,ireturn=ireturn,acc_limit=acc_limit,wlevel=1,level=0,level_max=level_max) 
                
                #search finished at this point
                sind_next=-self.sind_next                
                sind_list=array([array(i) for i in self.sind_list])                
                
                #clean
                delattr(self,'sind_next'); delattr(self,'flag_search'); delattr(self,'sind_list')
                
                if ireturn==5:
                    return sind_next
                elif ireturn==6:
                    return sind_list
                                     
            elif wlevel==1: #second level recursive search
                fpz=sind_next!=0                
                if level!=level_max:
                    if sum(fpz)==0:
                        #reach the end
                        fpn=self.sind_next>0; self.sind_next[fpn]=-sind0
                        if ireturn==6:
                            ind_list=nonzero(fpn)[0]; sind_list=sind0
                            [self.sind_list[i].append(j) for i,j in zip(ind_list,sind_list)]
                        self.flag_search=False
                    else:
                        #save the pts that reach the end
                        sind_next[~fpz]=-sind0[~fpz]
                        fpn=self.sind_next>0; self.sind_next[fpn]=sind_next
                        if ireturn==6:
                            ind_list=nonzero(fpn)[0]; sind_list=abs(sind_next)
                            [self.sind_list[i].append(j) for i,j in zip(ind_list,sind_list)]
                        
                        #continue search the rest pts 
                        self.search_upstream(sind_next[fpz],ireturn=ireturn,acc_limit=acc_limit,wlevel=1,level=level+1,level_max=level_max)
                        
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
                #delattr(self,'sind_next'); delattr(self,'flag_search'); delattr(self,'sind_list')
                
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
    
    def get_coordinates(self,sind0,affine=None,ds=None,nodata=None):     
        #example: sind0=nonzero(self.dir.ravel()==0)[0], or sind0=nonzero(self.dir==0)
        
        #check parameters first
        if ds is None: ds=self.ds
        if affine is None: affine=self.affine
        if nodata is None: nodata=self.nodata
        
        #check index
        if len(sind0)==2 and hasattr(sind0[0],'__len__'):            
            indy,indx=sind0  #sind0=[indy,indx]            
            fpn=(indy==nodata)|(indx==nodata)
        else:
            fpn=sind0==nodata; sind0[fpn]=0;
            indy,indx=unravel_index(sind0,ds)
            indy[fpn]=nodata; indx[fpn]=nodata
            
        #convert index to coordinates
        dxy=affine[0]; ly0=affine[5]; lx0=affine[2]
        
        xi=indx*dxy+lx0
        yi=-indy*dxy+ly0
        
        xi[fpn]=nodata; yi[fpn]=nodata
        
        return yi,xi
                    
         
if __name__=="__main__":    
    close('all')
    sys.setrecursionlimit(100000)


#-----------------------------------
    # S=dem_dir(); S.read_data('S1.npz'); S.compute_extent();
   
    # #check closed segments
    # # sind=nonzero(S.dir.ravel()==0)[0]; acc=S.acc.ravel()[sind];
    
    # sind0=ravel_multi_index(nonzero((S.seg<=10)*(S.seg>0)*(S.dir==0)),S.ds);   
    
    # sind_up=S.search_upstream(sind0,ireturn=5)
    # list_down=S.search_downstream(sind_up,ireturn=2)
    # list_up=S.search_upstream(sind0,ireturn=6)
    
    # dyi,dxi=S.get_coordinates(sind0)
    # uyi,uxi=S.get_coordinates(sind_up)
    
    # figure(); 
    
    # subplot(2,1,1)
    # imshow(S.acc,vmin=0,vmax=1e2,extent=S.extent)
    
    # for i in arange(len(sind0)):
    #     sindi=list_down[i];
    #     yi,xi=S.get_coordinates(sindi)
        
    #     fp=(xi==S.nodata)|(yi==S.nodata);
    #     xi[fp]=nan; yi[fp]=nan;
        
    #     plot(xi,yi,'r-')
        
    # plot(dxi,dyi,'b*',uxi,uyi,'w*')
        
    # subplot(2,1,2)
    # imshow(S.acc,vmin=0,vmax=1e2,extent=S.extent)
    
    # colors='rgb'
    # for i in arange(len(sind0)):
    #     sindi=list_up[i];
    #     yi,xi=S.get_coordinates(sindi)
        
    #     fp=(xi==S.nodata)|(yi==S.nodata);
    #     xi[fp]=nan; yi[fp]=nan;
        
    #     plot(xi,yi,'r-')
    
    # plot(dxi,dyi,'b*',uxi,uyi,'w*')

    
    
    
     
#------------------------plot watershed bnd------------------------------------
    # S=dem_dir(); S.read_data('S1_0.npz'); S.compute_extent();
   
    # #check closed segments
    # # sind=nonzero(S.dir.ravel()==0)[0]; acc=S.acc.ravel()[sind];
    
    # sind0=ravel_multi_index(nonzero((S.seg<=10)*(S.seg>0)*(S.dir==0)),S.ds); 
    
    # S.extract_river(sind0,acc_limit=1e2)
    
    # figure(); 
    # imshow(S.acc,vmin=0,vmax=1e2,extent=S.extent)
    
    # for i in arange(len(sind0)):
    #     sindi=S.rivers[i];
    #     yi,xi=S.get_coordinates(sindi)
        
    #     fp=(xi==S.nodata)|(yi==S.nodata);
    #     xi[fp]=nan; yi[fp]=nan;
        
    #     plot(xi,yi,'r-')
         
    
    # sind0=ravel_multi_index(nonzero((S.seg==1)),S.ds);
    # yi,xi=S.get_coordinates(sind0)
    # plot(xi,yi,'g.')
 
#--------------------------plot river with catchment---------------------------
    # S=dem_dir(); S.read_data('S2.npz')
    
    # #get xy
    # sind=nonzero(S.dir.ravel()==0)[0]; acc=S.acc.ravel()[sind]; fpt=acc>1e6;
    # yi,xi=S.get_coordinates(sind[fpt])
    
    # #plot
    # figure(); 
    # imshow(S.acc,vmin=0,vmax=1e5,extent=S.extent)
    # plot(xi,yi,'r.')
    
 
#--------------------------precalculation--------------------------------------
    # # # # # # #--------------------------------------------------------------
    # # compute dir and acc, then save them
    # # fname='./13arcs/southern_louisiana_13_navd88_2010.asc';
    # # fname='ne_atl_crm_v1.asc'; sname='S2_0'
    # fname='GEBCO.asc'; sname='S1_0'
        
    # #read dem info
    # gds=read_deminfo(fname,'dem',size_subdomain=1e7); 
    # gds.add_deminfo(fname,'dem0',size_subdomain=1e9); 
    
    # #read all data first
    # gds.read_dem_data('dem0')    

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
    #         try:                
    #             gd.fill_depressions(gd.data, out_name='fdem')
    #             try:
    #                 gd.resolve_flats(gd.fdem,'flats')
    #             except:
    #                 gd.flats=gd.fdem
    #         except:
    #             try:
    #                 gd.resolve_flats(gd.data,'flats')
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
    # S.nodata=gd.nodata; S.extent=gds.dem_info.extent; S.ds=S.dem.shape;  S.affine=[*gds.dem_info.affine]
    
    # #compute dir
    # t0=time.time();
    # S.compute_dir(S.dem,nodata=S.nodata);
    # S.compute_acc()
    # dt=time.time()-t0
    
    # S.save_data(sname,['dir','acc','seg','nodata','extent','ds','affine'])
    
    # #--------------------------------------------------------------------------
    
