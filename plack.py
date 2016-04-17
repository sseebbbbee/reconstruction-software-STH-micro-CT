import astra
import numpy as np
#from mayavi import mlab
import nibabel as nib
import scipy.io as sio
import glob
from own_functions import *
#from draw_test import *
from scipy import misc



# det=int(raw_input("Number of detectors: "))
# vol_size_ratio=int(raw_input("Vol size ratio: "))
# ang_space=int(raw_input("Angular spacing: "))

relaxation=0.05
lag_corr=True
det=int(raw_input("det size: "))
vol_size_ratio=int(raw_input("vol size ratio: "))
ang_space=1
# if det <= 750:
#     det_low=det
# else:
#     det_low=750
det_low=det
# vol_size_ratio=2
# ang_space=1
new_size=(int(det),int(det))

#plack='QRM_10_12_curr800_vol30_pulse1500'#raw_input("folder: ")
plack='plaque36kV800mA1500ms'
KTH=True
if KTH==True:
    filenames=glob.glob('CTdata/'+plack+'/*.bin') # makes filenames a list
    beam_profile_file='CTdata/lf36kV800mA1500ms/*.bin'
    dark_current_file='CTdata/dc36kV/*.bin'
    outfile='CTdata/Saved Data/proj_data_det_'+plack+'_det'+str(det)+'ang'+str(ang_space)
    outfile_low='CTdata/Saved Data/proj_data_det_'+plack+'_det'+str(det_low)+'ang'+str(ang_space)
else:
    import pylab
    filenames=glob.glob('C:/Users/Sebastian/Documents/CTdata/'+plack+'/*.bin') # makes filenames a list
    
    beam_profile_file="C:/Users/Sebastian/Documents/CTdata/lf36kV800mA1500ms/*.bin"
    dark_current_file="C:/Users/Sebastian/Documents/CTdata/dc36kV/*.bin"
    outfile='C:/Users/Sebastian/Documents/CTdata/Saved Data/proj_data_det_'+plack+'_det'+str(det)+'ang'+str(ang_space)
    outfile_low='C:/Users/Sebastian/Documents/CTdata/Saved Data/proj_data_det_'+plack+'_det'+str(det_low)+'ang'+str(ang_space)


###The geometry of the CT-system
with open('allignment parameters.txt') as f:
    lines = f.readlines()
param=np.zeros(5)
ii=0
for line in lines:
    inde=line.index('=')
    x=line[inde+1:len(line)-1]
    if x[0]=='-':
        param[ii]=-1*float(x[1:len(x)])
    else:
        param[ii]=float(x[0:len(x)])
    ii=ii+1
print param    
pixels=2400
source_to_origin_mm=param[0]
source_to_origin_pixels=np.round(source_to_origin_mm/0.05);
origin_to_detector_mm=param[1]-source_to_origin_mm
origin_to_detector_pixels=np.round(origin_to_detector_mm/0.05);
# changes the distances to fit the new sampling of the detector
origin_det=np.round(origin_to_detector_pixels/np.round(pixels/new_size[1]))
source_origin=np.round(source_to_origin_pixels/np.round(pixels/new_size[1]))
                        
alignment_in_mm=param[2:4] #dx, dz translation, third input is rotation, this part should be read from a txt file instead.
###
            
filenames.sort() # sorts the filenames so the come in the order they are obtained
number_of_files= len(filenames)
print number_of_files
number_of_projections=number_of_files/ang_space
print number_of_projections, "number of projections "
pixels_to_crop=150 # the number of pixels to crop on every side of the projection.
            
with open('Lag correction parameters.txt') as corr:
    lines = corr.readlines()
    
kV=45

a=np.zeros(3)
b=np.zeros(3)
for line in lines:
    if line==lines[0]:
        pass
    else:
        if int(line[0:2])==kV:
            indexes= find_char(line,',')
            a[0]=float(line[indexes[0]+1:indexes[1]])
            a[1]=float(line[indexes[1]+1:indexes[2]])
            a[2]=float(line[indexes[2]+1:indexes[3]])
            b[0]=float(line[indexes[3]+1:indexes[4]])
            b[1]=float(line[indexes[4]+1:indexes[5]])
            b[2]=float(line[indexes[5]+1:len(line)])

            print a, b


try:
    proj_data = sio.loadmat(outfile+'.mat')['proj_data']
    proj_data_low = sio.loadmat(outfile_low+'.mat')['proj_data']
except:
    proj_data,proj_data_low= making_sino_from_bin(filenames,lag_corr,number_of_files, number_of_projections,ang_space, pixels,beam_profile_file,dark_current_file,new_size,det_low,pixels_to_crop,alignment_in_mm,a,b)
    #sio.savemat(outfile+'.mat', {'proj_data':proj_data})
    #sio.savemat(outfile_low+'.mat', {'proj_data':proj_data_low})
print np.shape(proj_data[:,0,:])
origin_det=np.round(origin_to_detector_pixels/np.round(pixels/new_size[1]))
source_origin=np.round(source_to_origin_pixels/np.round(pixels/new_size[1]))
print origin_det
print source_origin
FOV= 2*int((det/2)/(origin_det + source_origin) * source_origin ) 
print FOV
nr_iterations_low=int(raw_input("Number of iterations low: "))
nr_iterations_high=int(raw_input("Number of iterations high: "))

vol_size=FOV

#vol_size_ratio=4
vol_size_low=vol_size/vol_size_ratio


overlap= 0.1 # overlap between regions needed because of  bad results in the border area.
nr_subvol_x=int(raw_input("Nr subvol in x dir.: "))
nr_subvol_y=int(raw_input("Nr subvol in y dir.: "))
nr_subvol_z=int(raw_input("Nr subvol in z dir.: "))
# nr_subvol_x=1
# nr_subvol_y=2
# nr_subvol_z=2
x_sub_boundaries,y_sub_boundaries,z_sub_boundaries,x_subvol_width,y_subvol_width,z_subvol_width=create_boundaries(nr_subvol_x,nr_subvol_y,nr_subvol_z,vol_size,overlap)


angles = np.linspace(0, 2*np.pi, number_of_projections,False)

# det_spacing_x=np.round(float(vol_size)/new_size[0],decimals=2) #change here if you want maginification of rec volume.
# det_spacing_y=np.round(float(vol_size)/new_size[0],decimals=2)
det_spacing_x=1
det_spacing_y=1
det_row_count=new_size[0]
det_col_count=new_size[0]

det_spacing_x_low=np.round(float(det_spacing_x)/vol_size_ratio,decimals=2) #0.5
det_spacing_y_low=np.round(float(det_spacing_y)/vol_size_ratio,decimals=2)

#algon='FDK_CGLS'
algon=raw_input("Algo low: ")
initializer=0
rec_volume=(vol_size_low,vol_size_low,vol_size_low)
det_size=(det,det_low)
if algon=='FDK':
    algo='FDK_CUDA'
    nr_iterations_low=1
    rec_low=make_reconsruction(proj_data_low,number_of_projections,vol_size_ratio,det_size,rec_volume,source_origin,origin_det,nr_iterations_low,algo,relaxation,initializer)
    rec_low=rec_low*(rec_low>0)
elif algon=='SIRT':
    algo='SIRT3D_CUDA'
    rec_low=make_reconsruction(proj_data_low,number_of_projections,vol_size_ratio,det_size,rec_volume,source_origin,origin_det,nr_iterations_low,algo,relaxation,initializer)
elif algon=='CGLS':
    algo='CGLS3D_CUDA'
    rec_low=make_reconsruction(proj_data_low,number_of_projections,vol_size_ratio,det_size,rec_volume,source_origin,origin_det,nr_iterations_low,algo,relaxation,initializer)
elif algon=='SART':
    algo='EGEN_SART'
    rec_low=make_reconsruction(proj_data,number_of_projections,vol_size_ratio,det_size,rec_volume,source_origin,origin_det,nr_iterations_low,algo,relaxation,initializer)
elif algon=='FDK_SIRT':
    algo='FDK_CUDA'
    rec_low=make_reconsruction(proj_data_low,number_of_projections,vol_size_ratio,det_size,rec_volume,source_origin,origin_det,nr_iterations_low,algo,relaxation,initializer)
    rec_low=rec_low*(rec_low>0)
    algo='SIRT3D_CUDA'
    rec_low=make_reconsruction(proj_data_low,number_of_projections,vol_size_ratio,det_size,rec_volume,source_origin,origin_det,nr_iterations_low,algo,relaxation,rec_low)
elif algon=='FDK_CGLS':
    algo='FDK_CUDA'
    rec_low=make_reconsruction(proj_data_low,number_of_projections,vol_size_ratio,det_size,rec_volume,source_origin,origin_det,nr_iterations_low,algo,relaxation,initializer)
    rec_low=rec_low*(rec_low>0)
    algo='CGLS3D_CUDA'
    rec_low=make_reconsruction(proj_data_low,number_of_projections,vol_size_ratio,det_size,rec_volume,source_origin,origin_det,nr_iterations_low,algo,relaxation,rec_low)
else:
    print'Algorithm not found, following algorithms exists: FDK,SIRT,SART,CGLS,FDK_SIRT and FDK_CGLS'
    quit()

proj_data_low=[]
save_as_low=plack+'_it'+str(nr_iterations_low)+'_'+algon
affine=np.eye(4) # affine transformation matrix when nothing is moved.
array_img=nib.Nifti1Image(rec_low, affine)
nib.save(array_img, save_as_low+'_low.nii')

##############################
################

# draw_answer = raw_input("Do you want to draw ROI (y/n): ")
# if draw_answer=='y':
#     im=misc.imresize(rec_low[0.5*vol_size_low,:,:], (vol_size,vol_size)) ## decide what interpolation
#     start_z,start_y,last_z,last_y=draw_ROI(im)
#  
#     z_overlap=np.round(overlap*(last_z-start_z))
#     z_subvol_width=last_z-start_z
#  
#     z_sub_boundaries=np.array([[start_z,start_z-z_overlap,start_z+z_overlap], [last_z, last_z-z_overlap, last_z+z_overlap]])
#  
#     nr_subvol_z=len(z_sub_boundaries)-1
#     y_overlap=np.round(overlap*(last_y-start_y))
#     y_subvol_width=last_y-start_y
#  
#     y_sub_boundaries=np.array([[start_y,start_y-y_overlap,start_y+y_overlap], [last_y, last_y-y_overlap, last_y+y_overlap]])
#  
#     nr_subvol_y=len(y_sub_boundaries)-1

################
if KTH==False:
    pylab.figure()
    pylab.gray()
    pylab.imshow(rec_low[0.50*vol_size_low,:,:])
    pylab.figure()
    pylab.imshow(rec_low[:,0.50*vol_size_low,:]) 
    pylab.show()
#

rec_low_temp=np.zeros((vol_size_low,vol_size_low,vol_size_low),dtype= np.float16)


##################################HERE Multi resolution starts
ii=0

x_parts,y_parts, z_parts=create_parts(nr_subvol_x,nr_subvol_y,nr_subvol_z) 

vol_geom_low = astra.create_vol_geom(vol_size_low, vol_size_low, vol_size_low)
proj_geom_low = astra.create_proj_geom('cone',  det_spacing_x_low, det_spacing_y_low, det_row_count, det_col_count, angles, source_origin, origin_det)
algon_high=raw_input("Algo high: ")
for n_x in x_parts:
    for n_y in y_parts:
        for n_z in z_parts: ## nr_subvol_z
            rec_low_temp[:,:,:]=rec_low[:,:,:]
            rec_low_temp[x_sub_boundaries[n_x,1] / vol_size_ratio:x_sub_boundaries[n_x+1,2] / vol_size_ratio, y_sub_boundaries[n_y,1] / vol_size_ratio:y_sub_boundaries[n_y+1,2] / vol_size_ratio, z_sub_boundaries[n_z,1] / vol_size_ratio:z_sub_boundaries[n_z+1,2] / vol_size_ratio] = 0

            proj_id_low, proj_data_low = astra.create_sino3d_gpu(rec_low_temp, proj_geom_low, vol_geom_low)
            astra.data3d.delete(proj_id_low)
            x_start=x_sub_boundaries[n_x,1]
            x_stop=x_sub_boundaries[n_x+1,2]
            h_start=FOV/2 - x_start
            h_stop=FOV/2 - x_stop  #height
            if x_start==0:
                sin_angle_start= float(FOV/2 - x_start + FOV/2*(det/2)/(origin_det + source_origin))/(source_origin - FOV/2) #not needed
                pro_start=0
            else:
                if h_start < 0:
                    sin_angle_start= float(FOV/2 - x_start)/(source_origin + FOV/2)
                else:
                    sin_angle_start= float(FOV/2 - x_start)/(source_origin - FOV/2)
                pro_start=np.floor(det/2 -(source_origin+origin_det)*sin_angle_start)
            if x_stop == vol_size:
                sin_angle_stop=float(FOV/2 - x_stop - FOV/2*(det/2)/(origin_det + source_origin))/(source_origin - FOV/2)
                pro_stop=det
            else:
                if h_stop < 0:
                    sin_angle_stop= float(FOV/2 - x_stop)/(source_origin - FOV/2)
                else:
                    sin_angle_stop= float(FOV/2 - x_stop)/(source_origin + FOV/2)
                pro_stop=np.ceil(det/2 -(source_origin+origin_det)*sin_angle_stop)
            det_move=det/2 - pro_stop + float(pro_stop -pro_start)/2
            det_move=-det_move
#             print x_start, x_stop,h_start,h_stop, sin_angle_start, sin_angle_stop, pro_start, pro_stop, det_move
#             pylab.figure()
#             pylab.imshow(proj_data_low[:, 0, :])

            #pylab.show()
            #rec_sub_vol=np.zeros((x_sub_boundaries[n_x+1,2]-x_sub_boundaries[n_x,1],y_sub_boundaries[n_y+1,2]-y_sub_boundaries[n_y,1],z_sub_boundaries[n_z+1,2]-z_sub_boundaries[n_z,1]))
            proj_data_part=proj_data[pro_start:pro_stop,:,:]
            proj_data_low_part=proj_data_low[pro_start:pro_stop,:,:]
            print np.shape(proj_data_part), 'shape  data part'
            print np.shape(proj_data_low_part),'shape data low part'
#             pylab.figure()
#             pylab.imshow(proj_data_part[:, 0, :])
#             pylab.figure()
#             pylab.imshow(proj_data_low_part[:, 0, :])

            rec_sub_vol =recon_multi_part(algon_high,proj_data_part,proj_data_low_part,det_move,angles,x_sub_boundaries,y_sub_boundaries,z_sub_boundaries,vol_size,det_spacing_x,det_spacing_y,x_subvol_width,y_subvol_width,z_subvol_width,n_x,n_y,n_z,nr_iterations_high,source_origin, origin_det,relaxation)
#                 pylab.figure()
#                 pylab.imshow(rec_sub_vol[0.65*vol_size,:,:])
#             pylab.figure()
#             pylab.imshow(rec_sub_vol[:,0.5*vol_size,:])
#             pylab.show()
#             print 'rec_sub shape', np.shape(rec_sub_vol)
#             print 'z start:', (z_sub_boundaries[n_z,0]-z_sub_boundaries[n_z,1])
#             print 'z slut:', (z_sub_boundaries[n_z,0]-z_sub_boundaries[n_z,1]+z_subvol_width)
            
            rec_sub_roi= rec_sub_vol[(x_sub_boundaries[n_x,0]-x_sub_boundaries[n_x,1]):(x_sub_boundaries[n_x,0]-x_sub_boundaries[n_x,1]+x_subvol_width),(y_sub_boundaries[n_y,0]-y_sub_boundaries[n_y,1]):(y_sub_boundaries[n_y,0]-y_sub_boundaries[n_y,1]+y_subvol_width),(z_sub_boundaries[n_z,0]-z_sub_boundaries[n_z,1]):(z_sub_boundaries[n_z,0]-z_sub_boundaries[n_z,1]+z_subvol_width)]
            #rec_ROI[:,:,:,ii]=rec_sub_roi
            print ii
            ii=ii+1


            
            #sio.savemat('rec_slice_'+algo+str(nr_iterations)+'_det'++'.mat', {'rec_slice':rec[270,:,:]})
            if KTH==False:
                outfile='C:/Users/Sebastian/Documents/CTdata/Saved Data/rec_ROI_'+str(ii)+'det_'+str(det)
            else:
                outfile='CTdata/Saved Data/rec_ROI_'+str(ii)+'det_'+str(det)+str(number_of_projections)
            sio.savemat(outfile+'.mat', {'rec_sub_roi':rec_sub_roi})
            #np.savez(outfile, rec_ROI=rec_ROI, width=width)
    
if KTH==False:
    full_rec=read_and_show(x_parts,y_parts,z_parts,x_sub_boundaries,y_sub_boundaries,z_sub_boundaries,det)
    save_as_low=raw_input("Save low rec woth name: ")
    affine=np.eye(4) # affine transformation matrix when nothing is moved.
    array_img=nib.Nifti1Image(rec_sub_roi, affine)
    nib.save(array_img, save_as_low +'_full.nii')
    pylab.figure()
    x_rec=len(x_parts)*x_subvol_width
    y_rec=len(y_parts)*y_subvol_width
    pylab.imshow(full_rec[0.5*x_rec,:,:])
    pylab.figure()
    pylab.imshow(full_rec[:,0.5*y_rec,:])
    pylab.show()

