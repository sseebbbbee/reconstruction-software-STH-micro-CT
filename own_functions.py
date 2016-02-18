import numpy as np
import astra
import scipy.io as sio
import glob
import pylab
#from skimage.transform import resize,rotate #Better, keeps real values if we want to rotate.
from scipy import misc

def find_char(s, ch):
    return [i for i, ltr in enumerate(s) if ltr == ch]

def get_mean_lf_and_dc(pixels,beam_profile_files,dark_current_files,lag_corr,a,b):
    # reads the bin files and creates a mean dark current and beam profile signal. corrections of defective elements are also done. 
    
    dc_len= len(dark_current_files)
    dc=np.zeros(pixels*pixels,dtype= np.float16)
    for i in range(0 ,dc_len):
        fname=dark_current_files[i]  
        B=np.fromfile(fname,dtype='int16')
        dc=dc + B[256:pixels*pixels+256]
    dc=dc/float(dc_len)
    max_pixels= np.where(dc>(np.mean(dc)+4*np.std(dc))) # finds pixels values that are larger than 4 standard deviations from the mean.
    max_row=max_pixels[0]
    for ii in max_row:
        dc[ii]=(dc[ii-10]+dc[ii-5]+dc[ii+5]+dc[ii+10])/4 # defective elements are interpolated.

    lf_len= len(beam_profile_files)
    lf=np.zeros(pixels*pixels,dtype= np.float16)
    if lag_corr==True: # if lag correction is chosen the beam profile also has to be corrected
        S1=np.zeros((pixels*pixels),dtype= np.float32)
        S2=np.zeros((pixels*pixels),dtype= np.float32)
        S3=np.zeros((pixels*pixels),dtype= np.float32)
        fname=beam_profile_files[0]
        B=np.fromfile(fname,dtype='int16')
        B=B[256:pixels*pixels+256]
        X=B-dc 
           

        S1=X+S1*np.exp(-a[0])
        S2=X+S2*np.exp(-a[1])
        S3=X+S3*np.exp(-a[2])
        
    else:
        fname=beam_profile_files[0]
        B=np.fromfile(fname,dtype='int16')
        B=B[256:pixels*pixels+256]
        X=B-dc        
    lf= lf +  X   
    for i in range(1 ,lf_len):
        fname=beam_profile_files[i]  
        B=np.fromfile(fname,dtype='int16')
        B=B[256:pixels*pixels+256]
        Y=B-dc  
        
        if lag_corr==True:
            X=Y-S1*b[0]*np.exp(-a[0])- S2*b[1]*np.exp(-a[1])- S3*b[2]*np.exp(-a[2])
            S1=X+S1*np.exp(-a[0])
            S2=X+S2*np.exp(-a[1])
            S3=X+S3*np.exp(-a[2])
        #X=X*(X>0)
        else:
            X=Y
                
        lf= lf +  X 
        
        
    lf=lf/float(lf_len) # getting the average beam profile
    hist, bin_edges = np.histogram(lf,bins=20, density=False)
    #correcting defective pixels.
    try: 
        hist_ind=np.where(hist < 100)  #Maybe good to change this number
        print hist
        print bin_edges
        hist_ind=hist_ind[0]
        hist_ind=hist_ind[0]
        threshold=bin_edges[hist_ind]

        max_pixels= np.where(lf>threshold) 
        max_row=max_pixels[0]
        for ii in max_row:
            lf[ii]=(lf[ii-10]+lf[ii-5]+lf[ii+5]+lf[ii+10])/4
    except:
        pass        
        
    return lf,dc

def correct_dead_pixels(projection,dead_row,dead_col):
    
    for ii in range(0,len(dead_row)): # might be good to have some safety procedure to look if index +-1 is in the image.
        projection[dead_row[ii],dead_col[ii]]=(projection[dead_row[ii]-1,dead_col[ii]]+projection[dead_row[ii]+1,dead_col[ii]]+projection[dead_row[ii],dead_col[ii]-1]+projection[dead_row[ii],dead_col[ii]+1])/4
    return projection

def correction_of_misalignment(Projection,pixels,pixels_to_crop,alignment_in_mm):
    ####This should be done for the whole matrix at once but now it's projections per projection#
    alignment_col_in_mm=-alignment_in_mm[0] #minus as Lorenzos coordinates system is that way.
    alignment_row_in_mm=alignment_in_mm[1]   #
    rotation=alignment_in_mm
    alignment_row_in_pixel=np.round(alignment_row_in_mm/0.05)
    alignment_col_in_pixel=np.round(alignment_col_in_mm/0.05)
    
#     if np.abs(alignment_row_in_pixel) > np.abs(alignment_col_in_pixel):   This can be used instead of just saying 50
#         crop_after_move=np.abs(alignment_row_in_pixel)+2
#     else:
#         crop_after_move=np.abs(alignment_col_in_pixel)+2
    crop_after_move=50 # used to insure that the dead pixels in the edges of the detector isn't included in the image.
    A=np.zeros((pixels,pixels),dtype= np.float16) #default is  numpy.float64
    #Projection=rotate(Projection,rotation,resize=False, order=5, preserve_range=True) rotation should probably be after translation
    A[pixels_to_crop+alignment_row_in_pixel:pixels-pixels_to_crop+alignment_row_in_pixel,pixels_to_crop+alignment_col_in_pixel:pixels-pixels_to_crop+alignment_col_in_pixel] = Projection[pixels_to_crop:pixels-pixels_to_crop,pixels_to_crop:pixels-pixels_to_crop]
    Projection=A[pixels_to_crop+crop_after_move:pixels-pixels_to_crop-crop_after_move,pixels_to_crop+crop_after_move:pixels-pixels_to_crop-crop_after_move] 
    return Projection


def making_sino_from_bin(filenames,lag_corr,number_of_files, number_of_projections,ang_space, pixels,beam_profile_file,dark_current_file,new_size,det_low,pixels_to_crop,alignment_in_mm,a,b):
    A=np.zeros((new_size[0],number_of_projections,new_size[1]),dtype= np.float16) # matrix with the highest resolution projection values
    D=np.zeros((det_low,number_of_projections,det_low),dtype= np.float16) # matrix with low resolution projection values

    ii=0;

    beam_profile_files=glob.glob(beam_profile_file) #get the filenames in the director specified in beam_profile_file
    dark_current_files=glob.glob(dark_current_file) #get the filenames in the director specified in beam_profile_file
    print
    if len(beam_profile_files)==1: # this part needed when I switched between using Lorenzos data and mine.
        Beam_profile=np.fromfile(beam_profile_file,dtype='int16')
        Dark_current=np.fromfile(dark_current_file,dtype='int16')
    else:
        Beam_profile,Dark_current=get_mean_lf_and_dc(pixels,beam_profile_files,dark_current_files,lag_corr,a,b)

    Beam_profile=np.reshape(Beam_profile,(pixels,pixels))
    Beam_profile=correction_of_misalignment(Beam_profile,pixels,pixels_to_crop,alignment_in_mm)

    ##### finding dead dead pixels in beam profile
    hist, bin_edges = np.histogram(Beam_profile,bins=20, density=False)
    hist_max=np.max(hist)
    ind_mean=np.where(hist==hist_max)
    ind_mean=ind_mean[0]
    beam_mean=bin_edges[ind_mean[0]-1] #the mean value of the beam profile.

    dead_pixels= np.where(Beam_profile<beam_mean/2) # another value than 2 is maybe better
    dead_row=dead_pixels[0]
    dead_col=dead_pixels[1]
    Beam_profile_max=np.max(Beam_profile)
    Beam_profile=correct_dead_pixels(Beam_profile,dead_row,dead_col)

    Beam_profile=misc.imresize(Beam_profile, new_size, mode='F') 
    Beam_profile_low=misc.imresize(Beam_profile, (det_low,det_low), mode='F') 

    #### start of getting the projection data
    fname=filenames[0]  
    B=np.fromfile(fname,dtype='int16')
    B=B[256:pixels*pixels+256] # here B is a vector 2400*2400, the first 256 values are not projection values.
    B=np.reshape(B-Dark_current,(pixels,pixels)) # makes data to 2D projection data
    X=correction_of_misalignment(B,pixels,pixels_to_crop,alignment_in_mm)
    X=correct_dead_pixels(X,dead_row,dead_col) #correcting pixels
    max_pixels= np.where(X>Beam_profile_max)
    max_row=max_pixels[0]
    max_col=max_pixels[1]
    X=correct_dead_pixels(X,max_row,max_col) #correcting hot pixels
    
    if lag_corr==True:
        numrows = len(X)
        numcols = len(X[0]) 
        S1=np.zeros((numrows,numcols),dtype= np.float32)
        S2=np.zeros((numrows,numcols),dtype= np.float32)
        S3=np.zeros((numrows,numcols),dtype= np.float32)
    else:
        print 'no lag correction'
    
    B=misc.imresize(X, new_size,mode='F') 
    E=misc.imresize(X, (det_low,det_low),mode='F') 
    #B=resize(X, new_size,order=1, preserve_range=True) 

    A[:,0,:]=-1*np.log(B/Beam_profile) ## getting attenuation coeffs
    D[:,0,:]=-1*np.log(E/Beam_profile_low) ## getting attenuation coeffs

    if lag_corr==True:
        S1=X+S1*np.exp(-a[0])
        S2=X+S2*np.exp(-a[1])
        S3=X+S3*np.exp(-a[2])

    ii=1;
    for i in range(1 ,number_of_files):
        fname=filenames[i]  
        B=np.fromfile(fname,dtype='int16')
        B=B[256:pixels*pixels+256]
        B=np.reshape(B-Dark_current,(pixels,pixels))
        Y=correction_of_misalignment(B,pixels,pixels_to_crop,alignment_in_mm)
        Y=correct_dead_pixels(Y,dead_row,dead_col)
        ### lag effect correction
        if lag_corr==True:
            X=Y-S1*b[0]*np.exp(-a[0])-S2*b[1]*np.exp(-a[1])-S3*b[2]*np.exp(-a[2])
            X=X*(X>0) # to make sure no negative values.
        else:
            X=Y
            X=X*(X>0)
        max_pixels= np.where(X>Beam_profile_max)
        max_row=max_pixels[0]
        max_col=max_pixels[1]
        X=correct_dead_pixels(X,max_row,max_col)
        if lag_corr==True:
            S1=X+S1*np.exp(-a[0])
            S2=X+S2*np.exp(-a[1])
            S3=X+S3*np.exp(-a[2])
        ##########
        
        if i%ang_space==0: # if you have taken more projections than you want to reconstruct from then a ang space > 1 can be put.
#             pylab.figure()
#             pylab.imshow(X)
#             pylab.show()
            B=misc.imresize(X, new_size,mode='F')
            hist, bin_edges = np.histogram(B,bins=20, density=False)

            E=misc.imresize(X, (det_low,det_low),mode='F') 
            A[:,ii,:]=-1*np.log(B/Beam_profile) ## getting attenuation coeffs
            D[:,ii,:]=-1*np.log(E/Beam_profile_low) ## getting attenuation coeffs

            ii=ii+1
            print i,ii

     
    return A,D

def make_reconsruction(proj_data,number_of_projections,vol_size_ratio,det_size,rec_volume,source_to_origin_pixels,origin_to_detector_pixels,nr_iterations,algo,relaxation,initializer):

    vol_geom = astra.create_vol_geom(rec_volume)
    
    angles = np.linspace(0, 2*np.pi, number_of_projections, False)
    det_spacing=np.round((float(det_size[0])/det_size[1])/vol_size_ratio,decimals=3)

    proj_geom = astra.create_proj_geom('cone', det_spacing, det_spacing, det_size[1], det_size[1], angles, source_to_origin_pixels,origin_to_detector_pixels)

    #np.round(float(rec_volume[0])/new_size[0],decimals=2)
    if algo[0:4] != 'EGEN':
        # Create projection data from this
        proj_id=4

        # Create a data object for the reconstruction
        rec_id = astra.data3d.create('-vol', vol_geom,initializer)
        sinogram_id = astra.data3d.create('-proj3d', proj_geom, proj_data)
        #Set up the parameters for a reconstruction algorithm using the GPU
        cfg = astra.astra_dict(algo)
        cfg['ReconstructionDataId'] = rec_id
        cfg['ProjectionDataId'] = sinogram_id
        cfg['option']={}
        cfg['option']['MinConstraint'] = 0 #attenuation coefficient can't be negative
 
        # Create the algorithm object from the configuration structure
        alg_id = astra.algorithm.create(cfg)

    print "startar algo"
    if algo == 'SIRT3D_CUDA' or algo=='CGLS3D_CUDA':
        astra.algorithm.run(alg_id,nr_iterations)
    elif algo=='FDK_CUDA':
        astra.algorithm.run(alg_id)
    elif algo=='EGEN_SIRT':
        if isinstance( initializer, ( int, long ) ):
            rec_volume=np.zeros(rec_volume,dtype= np.float16)
        else:
            rec_volume=initializer
        for ii in range(0,nr_iterations):

            sinogram_id, proj_it = astra.create_sino3d_gpu(rec_volume, proj_geom,vol_geom)
            residual=proj_it-proj_data

            h=relaxation*float(1)/(rec_volume.shape[0]*number_of_projections)
            [id, rec_volume_it] = astra.create_backprojection3d_gpu(residual, proj_geom, vol_geom)
            rec_volume= rec_volume -h*rec_volume_it
            rec_volume=rec_volume*(rec_volume>0)
            astra.data3d.delete(id)
            astra.data3d.delete(sinogram_id)
            
    elif algo=='EGEN_SART':
        if isinstance( initializer, ( int, long ) ):
            rec_volume=np.zeros(rec_volume,dtype= np.float16)
        else:
            rec_volume=initializer
        #print rec_volume.shape[0]
        
        for ii in range(0,nr_iterations):
            angles_ind = np.linspace(0, number_of_projections, number_of_projections,False)
            np.random.shuffle(angles_ind)

            for jj in angles_ind:
                proj_geom = astra.create_proj_geom('cone', det_spacing, det_spacing, det_size[1], det_size[1], angles[jj], source_to_origin_pixels,origin_to_detector_pixels)
            #print np.min(sub_volume_high)
                sinogram_id, proj_it = astra.create_sino3d_gpu(rec_volume, proj_geom,vol_geom)
                proj_it=proj_it[:,0,:]
                residual=np.zeros((det_size[1],1,det_size[1]))
                residual[:,0,:]=proj_it-proj_data[:,jj,:]
                #print np.shape(proj_it),np.shape(residual),np.shape(proj_data[:,jj,:])
                #astra.data3d.store(sinogram_id, residual)
                
                h=relaxation*float(1)/(rec_volume.shape[0])#0.00001
                [id, rec_volume_it] = astra.create_backprojection3d_gpu(residual, proj_geom, vol_geom)
                rec_volume= rec_volume -h*rec_volume_it
                rec_volume=rec_volume*(rec_volume>0)
                astra.data3d.delete(id)
                astra.data3d.delete(sinogram_id)    
#             pylab.figure()
#             pylab.gray()
#             pylab.imshow(rec_volume[:,:,128])
            

    else:
        print"algorithm is not supported"
        
        
    print "slut algo"
    if algo[0:4] != 'EGEN':
        # Get the result
        rec = astra.data3d.get(rec_id)
    
        # Clean up. Note that GPU memory is tied up in the algorithm object,
        # and main RAM in the data objects.
        astra.algorithm.delete(alg_id)
        astra.data3d.delete(rec_id)
        astra.data3d.delete(proj_id)
        astra.data3d.delete(sinogram_id)
        return rec
    else:
        return rec_volume





def create_boundaries(nr_subvol_x,nr_subvol_y,nr_subvol_z,vol_size,overlap):
    # Creating the boundaries for the multireosolution technique
    x_sub_boundaries=np.zeros((nr_subvol_x+1,3)) # ind 0 = subvolume boundarier, 1 overlap to the volume below, 2 overlap to volume above
    x_subvol_width=vol_size/nr_subvol_x
    x_overlap=np.round(vol_size*overlap)
    for ii in range(0,nr_subvol_x+1):
        if ii == nr_subvol_x:
            x_sub_boundaries[ii,0]=vol_size
        else:
            x_sub_boundaries[ii,0]=ii*x_subvol_width
    for ii in range(0,nr_subvol_x+1):
        if x_sub_boundaries[ii,0]- x_overlap > 0:
            x_sub_boundaries[ii,1]=x_sub_boundaries[ii,0]- x_overlap
    for ii in range(0,nr_subvol_x+1):
        if x_sub_boundaries[ii,0]+ x_overlap < vol_size and ii != nr_subvol_x:
            x_sub_boundaries[ii,2]=x_sub_boundaries[ii,0]+ x_overlap
        else:
            x_sub_boundaries[ii,2]=vol_size
    print x_sub_boundaries
    y_sub_boundaries=np.zeros((nr_subvol_y+1,3)) # ind 0 = subvolume boundarier, 1 overlap to the volume below, 2 overlap to volume above
    y_subvol_width=vol_size/nr_subvol_y
    y_overlap=np.round(vol_size*overlap)
    for ii in range(0,nr_subvol_y+1):
        if ii == nr_subvol_y:
            y_sub_boundaries[ii,0]=vol_size
        else:
            y_sub_boundaries[ii,0]=ii*y_subvol_width
    for ii in range(0,nr_subvol_y+1):
        if y_sub_boundaries[ii,0]- y_overlap > 0:
            y_sub_boundaries[ii,1]=y_sub_boundaries[ii,0]- y_overlap
    for ii in range(0,nr_subvol_y+1):
        if y_sub_boundaries[ii,0]+ y_overlap < vol_size and ii != nr_subvol_y:
            y_sub_boundaries[ii,2]=y_sub_boundaries[ii,0]+ y_overlap
        else:
            y_sub_boundaries[ii,2]=vol_size
    print y_sub_boundaries
    z_sub_boundaries=np.zeros((nr_subvol_z+1,3)) # ind 0 = subvolume boundarier, 1 overlap to the volume below, 2 overlap to volume above
    z_subvol_width=vol_size/nr_subvol_z
    z_overlap=np.round(vol_size*overlap)
    for ii in range(0,nr_subvol_z+1):
        if ii == nr_subvol_z:
            z_sub_boundaries[ii,0]=vol_size
        else:
            z_sub_boundaries[ii,0]=ii*z_subvol_width
    for ii in range(0,nr_subvol_z+1):
        if z_sub_boundaries[ii,0]- z_overlap > 0:
            z_sub_boundaries[ii,1]=z_sub_boundaries[ii,0]- z_overlap
    for ii in range(0,nr_subvol_z+1):
        if z_sub_boundaries[ii,0]+ z_overlap < vol_size and ii != nr_subvol_z:
            z_sub_boundaries[ii,2]=z_sub_boundaries[ii,0]+ z_overlap
        else:
            z_sub_boundaries[ii,2]=vol_size
    print z_sub_boundaries
    return x_sub_boundaries,y_sub_boundaries,z_sub_boundaries,x_subvol_width,y_subvol_width,z_subvol_width

def create_parts(nr_subvol_x,nr_subvol_y,nr_subvol_z):
    #This funtion is not used in the GUI version
    if nr_subvol_x >1:
        con=True
        while con==True:
            x_answer = raw_input("What x part do you want, all or number for part: ")
            if x_answer=="all":
                x_parts=np.linspace(0,nr_subvol_x-1,nr_subvol_x)
                con=False
            elif len(str(x_answer))==1 and int(x_answer[0]) < nr_subvol_x:
                x_parts=np.linspace(int(x_answer),int(x_answer),1)
                con=False
            elif len(str(x_answer))==2 and int(x_answer) < nr_subvol_x:
                x_parts=np.linspace(int(x_answer),int(x_answer),1)
                con=False
            elif len(x_answer)==3 and x_answer[0].isdigit() and x_answer[2].isdigit() and int(x_answer[0]) < nr_subvol_x and int(x_answer[2]) < nr_subvol_x:
                first_x=int(x_answer[0])
                last_x=int(x_answer[2])
                x_parts=np.linspace(first_x,last_x,last_x-first_x+1)
                con=False
            else:
                print 'Incorrect input, Try again'
    else:
        x_parts=np.linspace(0,0,1)
    
    print x_parts 
    if nr_subvol_y >1:
        con=True
        while con==True:
            y_answer = raw_input("What y part do you want, all or number for part: ")
            if y_answer=="all":
                y_parts=np.linspace(0,nr_subvol_y-1,nr_subvol_y)
                con=False
            elif len(str(y_answer))==1 and int(y_answer[0]) < nr_subvol_y:
                y_parts=np.linspace(int(y_answer),int(y_answer),1)
                con=False
            elif len(y_answer)==3 and y_answer[0].isdigit() and y_answer[2].isdigit() and int(y_answer[0]) < nr_subvol_y and int(y_answer[2]) < nr_subvol_y:
                first_y=int(y_answer[0])
                last_y=int(y_answer[2])
                y_parts=np.linspace(first_y,last_y,last_y-first_y+1)
                con=False
            else:
                print 'Incorrect input, Try again'
    else:
        y_parts=np.linspace(0,0,1)    
    print y_parts 

    if nr_subvol_z >1:
        con=True
        while con==True:
            z_answer = raw_input("What z part do you want, all or number for part: ")
            if z_answer=="all":
                z_parts=np.linspace(0,nr_subvol_z-1,nr_subvol_z)
                con=False
            elif len(str(z_answer))==1 and int(z_answer[0]) < nr_subvol_z:
                z_parts=np.linspace(int(z_answer),int(z_answer),1)
                con=False
            elif len(z_answer)==3 and z_answer[0].isdigit() and z_answer[2].isdigit() and int(z_answer[0]) < nr_subvol_z and int(z_answer[2]) < nr_subvol_z:
                first_z=int(z_answer[0])
                last_z=int(z_answer[2])
                z_parts=np.linspace(first_z,last_z,last_z-first_z+1)
                con=False
            else:
                print 'Incorrect input, Try again'
    else:
        z_parts=np.linspace(0,0,1)      
    print z_parts
    return x_parts,y_parts, z_parts

def how_to_move(x_sub_boundaries,y_sub_boundaries,z_sub_boundaries,vol_size,x_subvol_width,y_subvol_width,z_subvol_width,n_x,n_y,n_z):
    # in Astra 1.6 the reconstruction volume is centered around the origin. But the detector and source is moveable
    # due to this limitations we must calculate how the source and x-ray should be moved to correspond to
    # the geometry as the projection data is obtained from. 
    if z_sub_boundaries[n_z,0] != 0 and z_sub_boundaries[n_z+1,2] != vol_size:
        move_x = float(vol_size) / 2 - z_sub_boundaries[n_z,0] - float(z_subvol_width) / 2
    elif z_sub_boundaries[n_z,0] == 0:
        move_x = float(vol_size) / 2 - z_sub_boundaries[n_z,0] - float(z_sub_boundaries[n_z+1,2] - z_sub_boundaries[n_z,1]) / 2
    elif z_sub_boundaries[n_z+1,2] == vol_size:
        move_x = float(vol_size) / 2 - z_sub_boundaries[n_z,1] - float(z_sub_boundaries[n_z+1,2] - z_sub_boundaries[n_z,1]) / 2       
    else:
        move_x=0
    if y_sub_boundaries[n_y,0] != 0 and y_sub_boundaries[n_y+1,2] != vol_size:
        move_y = float(vol_size) / 2 - y_sub_boundaries[n_y,0] - float(y_subvol_width) / 2
    elif y_sub_boundaries[n_y,0] == 0:
        move_y = float(vol_size) / 2 - y_sub_boundaries[n_y,0] - float(y_sub_boundaries[n_y+1,2] - y_sub_boundaries[n_y,1]) / 2
    elif y_sub_boundaries[n_y+1,2] == vol_size:
        move_y = float(vol_size) / 2 - y_sub_boundaries[n_y,1] - float(y_sub_boundaries[n_y+1,2] - y_sub_boundaries[n_y,1]) / 2
    else:
        move_y=0

    if x_sub_boundaries[n_x,0] != 0 and x_sub_boundaries[n_x+1,2] != vol_size:
        move_z = float(vol_size) / 2 - x_sub_boundaries[n_x,0] - float(x_subvol_width) / 2
    elif x_sub_boundaries[n_x,0] == 0:
        move_z = float(vol_size) / 2 - x_sub_boundaries[n_x,0] - float(x_sub_boundaries[n_x+1,2] - x_sub_boundaries[n_x,1]) / 2
    elif x_sub_boundaries[n_x+1,2] == vol_size:
        move_z = float(vol_size) / 2 - x_sub_boundaries[n_x,1] - float(x_sub_boundaries[n_x+1,2] - x_sub_boundaries[n_x,1]) / 2
    else:
        move_z=0

    return move_x,move_y,move_z

def read_and_show(x_parts,y_parts,z_parts,x_sub_boundaries,y_sub_boundaries,z_sub_boundaries,det):
    # puts the subvolumes at the right place and form a full reconstruction.
    nr_x_parts=len(x_parts)
    nr_y_parts=len(y_parts)
    nr_z_parts=len(z_parts)
    print 'x_parts', x_parts
    x_rec=0
    for i in x_parts:
        x_rec=x_rec + x_sub_boundaries[i+1,0] - x_sub_boundaries[i,0]   
    y_rec=0
    for i in y_parts:
        y_rec=y_rec + y_sub_boundaries[i+1,0] - y_sub_boundaries[i,0]  
    z_rec=0
    for i in z_parts:
        z_rec=z_rec + z_sub_boundaries[i+1,0] - z_sub_boundaries[i,0]  

        
    print z_rec, y_rec
    full_rec=np.zeros((x_rec,y_rec,z_rec),np.float16)
    print np.shape(full_rec)
    ii=1
    for n_x in range(0,len(x_parts)):
        for n_y in range(0,len(y_parts)):
            for n_z in range(0,len(z_parts)):
                #full_rec[x_sub_boundaries[n_x,0]:x_sub_boundaries[n_x+1,0],y_sub_boundaries[n_y,0]:y_sub_boundaries[n_y+1,0],z_sub_boundaries[n_z,0]:z_sub_boundaries[n_z+1,0]]=rec_ROI[:,:,:,ii]
                outfile='C:/Users/Sebastian/Documents/CTdata/Saved Data/rec_ROI_'+str(ii)+'det_'+str(det)+'.mat'
                P = sio.loadmat(outfile)['rec_sub_roi']
                print np.shape(P)
                if (nr_x_parts*nr_y_parts*nr_z_parts)==1:
                    full_rec=P
                else:
                    full_rec[x_sub_boundaries[n_x,0]:x_sub_boundaries[n_x+1,0],y_sub_boundaries[n_y,0]:y_sub_boundaries[n_y+1,0],z_sub_boundaries[n_z,0]:z_sub_boundaries[n_z+1,0]]=P
                ii=ii+1
             

    return full_rec

def recon_multi_part(algon,proj_data,proj_data_low,det_move,angles,x_sub_boundaries,y_sub_boundaries,z_sub_boundaries,vol_size,det_spacing_x,det_spacing_y,x_subvol_width,y_subvol_width,z_subvol_width,n_x,n_y,n_z,nr_iterations_high,source_origin, origin_det,relaxation):
    # this part performs the reconstruction of the subvolume.
    pro_proj = proj_data - proj_data_low  # we only want to have the projection data that comes from the sub volume
    det_row_count= pro_proj.shape[0]
    det_col_count=pro_proj.shape[2]

    vectors_ROI = np.zeros((len(angles), 12)) # this vector will contain the information of the source and detector position.
    move_x,move_y,move_z=how_to_move(x_sub_boundaries,y_sub_boundaries,z_sub_boundaries,vol_size,x_subvol_width,y_subvol_width,z_subvol_width,n_x,n_y,n_z)
    for i in range(0,len(angles)):
        # source
        vectors_ROI[i,0] = np.sin(angles[i]) * (source_origin) +move_x 
        vectors_ROI[i,1] = -np.cos(angles[i]) * (source_origin) + move_y
        vectors_ROI[i,2] = move_z
        # center of detector
        vectors_ROI[i,3] = -np.sin(angles[i]) * (origin_det) + move_x 
        vectors_ROI[i,4] = np.cos(angles[i]) * (origin_det)  + move_y
        vectors_ROI[i,5] = move_z +det_move
        # vector from detector pixel 0 to 1
        vectors_ROI[i,6] = np.cos(angles[i]) *det_spacing_x    
        vectors_ROI[i,7] = np.sin(angles[i])*det_spacing_x
        vectors_ROI[i,8] = 0

        # vector from detector pixel (0,0) to (1,0)
        vectors_ROI[i,9] = 0
        vectors_ROI[i,10] = 0
        vectors_ROI[i,11] = det_spacing_y
    
            #ROI=cube[:,:,z_sub_boundaries[n_z,0]:z_sub_boundaries[n_z+1,0]]
            # mlab.figure()
            # mlab.contour3d(ROI)
    width_z=z_sub_boundaries[n_z+1,2]-z_sub_boundaries[n_z,1]
    width_y=y_sub_boundaries[n_y+1,2]-y_sub_boundaries[n_y,1]
    width_x=x_sub_boundaries[n_x+1,2]-x_sub_boundaries[n_x,1]
    
    vol_geom_ROI= astra.create_vol_geom(width_y.astype(int), width_z.astype(int), width_x.astype(int))
    proj_geom_ROI = astra.create_proj_geom('cone_vec', det_row_count, det_col_count, vectors_ROI)

    if algon[0:4] != 'EGEN':

        # Create a data object for the reconstruction


        rec_id_ROI = astra.data3d.create('-vol', vol_geom_ROI)
        sinogram_id_ROI = astra.data3d.create('-proj3d', proj_geom_ROI, pro_proj)
        if algon=='SIRT':
            algon='SIRT3D_CUDA'
        elif algon=='CGLS':
            algon='CGLS3D_CUDA'

        # Set up the parameters for a reconstruction algorithm using the GPU
        cfg = astra.astra_dict(algon)   
        cfg['ReconstructionDataId'] = rec_id_ROI
        cfg['ProjectionDataId'] = sinogram_id_ROI
        cfg['option']={}
        cfg['option']['MinConstraint'] = 0 #attenuation coefficient can't be negative
        #cfg['option']['VoxelSuperSampling'] = np.round(float(det_spacing_x)/vol_size,decimals=2)
        alg_id = astra.algorithm.create(cfg)

    print "startar algo"
    if algon == 'SIRT3D_CUDA' or algon=='CGLS3D_CUDA':
        astra.algorithm.run(alg_id,nr_iterations_high)

            
    elif algon=='EGEN_SART':

        rec_volume=(width_x.astype(int), width_y.astype(int), width_z.astype(int))
        rec_volume=np.zeros(rec_volume,dtype= np.float16)
        
        sides=[rec_volume.shape[0],rec_volume.shape[1],rec_volume.shape[2]]
        max_side=np.max(sides)
        print max_side
        h=relaxation*float(1)/(max_side)#0.00001 change to max size of y o z


        vectors=np.matrix('0 0 0 0 0 0 0 0 0 0 0 0',dtype= np.float16)
        for ii in range(0,nr_iterations_high):
            angles_ind = np.linspace(0, len(angles), len(angles),False)
            np.random.shuffle(angles_ind) #shuffle the projection angles to get faster convergence
            
            for jj in angles_ind:
                
                #proj_geom = astra.create_proj_geom('cone', det_spacing, det_spacing, det_size[1], det_size[1], angles[jj], source_to_origin_pixels,origin_to_detector_pixels)
                vectors[:,:]=vectors_ROI[jj,:]
                proj_geom_SART = astra.create_proj_geom('cone_vec', det_row_count, det_col_count, vectors)
                sinogram_id, proj_it = astra.create_sino3d_gpu(rec_volume, proj_geom_SART,vol_geom_ROI)

                proj_it=proj_it[:,0,:]


                residual=np.zeros((det_row_count,1,det_col_count))
                residual[:,0,:]=proj_it-pro_proj[:,jj,:]

                astra.data3d.store(sinogram_id, residual)

                #h=0.001
                [id, rec_volume_it] = astra.create_backprojection3d_gpu(residual, proj_geom_SART, vol_geom_ROI)
                rec_volume= rec_volume -h*rec_volume_it
                rec_volume=rec_volume*(rec_volume>0)
                astra.data3d.delete(id)
                astra.data3d.delete(sinogram_id)  
                
            

    else:
        print"algorithm is not supported"

    print "slut algo"
    if algon[0:4] != 'EGEN':
        # Get the result
        
        rec_sub_vol=astra.data3d.get(rec_id_ROI)
    
        # Clean up. Note that GPU memory is tied up in the algorithm object,
        # and main RAM in the data objects.
        astra.algorithm.delete(alg_id)
        astra.data3d.delete(rec_id_ROI)
        astra.data3d.delete(sinogram_id_ROI)
        return rec_sub_vol
    else:
        return rec_volume



    
                        
