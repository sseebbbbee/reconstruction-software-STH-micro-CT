from Tkinter import *
import scipy.io as sio
import numpy as np
import tkFileDialog
import ttk
import pylab
from PIL import Image, ImageTk
import glob
from own_functions import *
import nibabel as nib

#After uninstalling pillow and PIL:
#pip install image
#Fixed this issue.

#The class draw is used for drawing a region of interest in the
#low reconstrruction that will be used in the multiresolution technique.
class draw(Toplevel):
    def __init__(self,im,ROI,vol_size,vol_size_ratio,*rec_vol):
        # im: is the image that you draw a region on.
        # ROI: can be true or false depending on if you have the image of the parallel projection in a window already
        # rec_vol: is the low reconstruction 
        self.vol_size=vol_size
        self.vol_size_ratio=vol_size_ratio
        Toplevel.__init__(self)

        self.x = self.y = 0
        self.ROI=ROI
        try:
            self.rec_vol=rec_vol[0]
        except:
            pass
        w, h = im.size
     
        self.canvas = Canvas(self, width=w, height=h, cursor="cross")
        self.canvas.pack()
        self.canvas.bind("<ButtonPress-1>", self.on_button_press)
        self.canvas.bind("<B1-Motion>", self.on_move_press)
        self.canvas.bind("<ButtonRelease-1>", self.on_button_release)
        
        self.im=ImageTk.PhotoImage(im,master = self)
        self.rect = None
        self.start_x = None
        self.start_y = None
        

        self.draw_image()
        if self.ROI==False: #it's only posible to draw the region viewing from two direction
            self.draw_button=Button(self, text="Draw ROI")
            self.draw_button.pack()
            self.draw_button.bind("<ButtonPress-1>", self.draw_yz_ROI)
        
        self.quit_button=Button(self, text="Quit")
        self.quit_button.pack()
        self.quit_button.bind("<ButtonPress-1>", self.quitting)

    def draw_image(self):
        #puts either the reconstruction slice or parallel projection as background image
        self.canvas.create_image(0,0,anchor="nw",image=self.im)
        
    def on_button_press(self, event):
        # saves mouse start position
        self.start_x = event.x
        self.start_y = event.y
        # create rectangle if not yet exist
        if not self.rect:
            self.rect = self.canvas.create_rectangle(self.x, self.y, 1, 1, outline="red")

    def on_move_press(self, event):
         
        self.last_x, self.last_y = (event.x, event.y)
        #change the variables so it fits the reconstruction system
        if self.start_y < self.last_y:
            self.rec_start_x=self.start_y
            self.rec_last_x=self.last_y
        else:
            self.rec_start_x=self.last_y
            self.rec_last_x=self.start_y
        # expand rectangle as you drag the mouse
        self.canvas.coords(self.rect, self.start_x, self.start_y, self.last_x, self.last_y)

    def on_button_release(self, event):
        pass
    def draw_yz_ROI(self,event):
    # code to exit
        #self.rec_vol[self.last_y/self.vol_size_ratio,:,:] is the slice shown in the second drawing window
        P=misc.imresize(self.rec_vol[self.last_y/self.vol_size_ratio,:,:], (self.vol_size,self.vol_size)) # makes the images same size as the reconstruction volume
        max_val=float(np.max(P))
        phan=Image.fromarray(P/(max_val/255)) # For visualisation
        ROI=True #so the draw class knows that Region of interest is drwan in x-y plane.
        draw_ROI = draw(phan,ROI,self.vol_size,self.vol_size_ratio)
        
        draw_ROI.mainloop()
        
        # sorting the values so they have the lowest value at start, so you can draw the rectangle in any way you want.
        if draw_ROI.start_x < draw_ROI.last_x:
            self.rec_start_z=draw_ROI.start_x
            self.rec_last_z=draw_ROI.last_x
        else:
            self.rec_start_z=draw_ROI.last_x
            self.rec_last_z=draw_ROI.start_x
        if draw_ROI.start_y < draw_ROI.last_y:
            self.rec_start_y=draw_ROI.start_y
            self.rec_last_y=draw_ROI.last_y
        else:
            self.rec_start_y=draw_ROI.last_y
            self.rec_last_y=draw_ROI.start_y
        
    def quitting(self,event):
        #closes the window where you can draw
        self.quit()
        self.destroy()

class reconstruction_GUI:

    def __init__(self,top):

        self.overlap=0.10 # the overlap used in the multi resolution technique.
        
        #### creating the variables and the input boxes in the GUI that has to do with projection data ######
        self.dirname_data=StringVar()
        self.dirname_lf=StringVar()
        self.dirname_dc=StringVar()
        
        self.L_data = Label(top, textvariable=self.dirname_data)
        self.L_data.grid(row=0,column=1,columnspan=7) 
        self.button_dir_data=Button(top, text='Load projection data    ',command= lambda: self.askdirectory('proj'))
        self.button_dir_data.grid(row=0,column=0,sticky=W)#Sticky E=east, W=west, S=south and N=North
        
        self.L_lf = Label(top, textvariable=self.dirname_lf)
        self.L_lf.grid(row=1,column=1,columnspan=7) #Sticky E=east, W=west, S=south and N=North
        self.button_dir_lf=Button(top, text='Load beam profile data',command= lambda: self.askdirectory('beam')) 
        self.button_dir_lf.grid(row=1,column=0,sticky=W)
        
        self.L_dc = Label(top, textvariable=self.dirname_dc)
        self.L_dc.grid(row=2,column=1,columnspan=7) #Sticky E=east, W=west, S=south and N=North
        self.button_dir_dc=Button(top, text='Load dark current data ',command= lambda: self.askdirectory('dark')) 
        self.button_dir_dc.grid(row=2,column=0,sticky=W)

        self.lag_corr=BooleanVar()
        self.lag_corr_check= Checkbutton(top, text='Lag effect coorrection',variable=self.lag_corr)
        self.lag_corr_check.grid(row=3,column=0, columnspan=2,sticky=W)
        
        self.nr_det_label = Label(top, text="Nr. of detectors (max 2000):")
        self.nr_det_label.grid(row=3,column=1,columnspan=2,sticky=E)
        self.nr_det_entry = Entry(top, width=5)
        self.nr_det_entry.grid(row=3,column=3,sticky=W)
        
 
        self.vol_ratiobox = ttk.Combobox(top, values=('1','2', '3', '4','6','8'),state='readonly')
        self.vol_ratiobox.grid(row=3,column=5,columnspan=3,sticky=W)
        self.vol_ratio_label = Label(top, text="Vol. ratio:")
        self.vol_ratio_label.grid(row=3,column=4)
        
        self.button_get_proj=Button(top, text='Process data')
        self.button_get_proj.bind("<ButtonPress-1>",lambda event: self.proj_status_text(event,'process'))
        self.button_get_proj.bind("<ButtonRelease-1>",self.get_proj)
        self.button_get_proj.grid(row=4,column=0)
        self.proj_status=StringVar()
        self.label_get_proj_status=Label(top,textvariable=self.proj_status)
        self.label_get_proj_status.grid(row=4,column=1,columnspan=4,sticky=W)
        
        #####################
        
        #### creating the variables and the input boxes in the GUI that has to do with the low resolution reconstruction######
        
        self.algo = StringVar()
        self.algobox = ttk.Combobox(top, values=('FDK', 'SIRT', 'SART', 'CGLS','FDK_SIRT','FDK_CGLS'),state='readonly')
        self.algobox.bind("<<ComboboxSelected>>", self.set_comboparam)
        self.algobox.grid(row=5,column=1,columnspan=3,sticky=W)
        self.algo_label = Label(top, text="Algorithm:")
        self.algo_label.grid(row=5,column=0)

        self.nr_low_label = Label(top, text="Nr. of iterations low resolution:")
        self.nr_low_label.grid(row=6,column=0)
        self.nr_low_entry = Entry(top, width=5)
        self.nr_low_entry.grid(row=6,column=1,sticky=W)

        
        self.button_rec_low=Button(top, text='Start low resolution reconstruction')
        self.button_rec_low.bind("<ButtonPress-1>",lambda event: self.proj_status_text(event,'rec_low'))
        self.button_rec_low.bind("<ButtonRelease-1>",self.start_low_rec)
        self.button_rec_low.grid(row=7,column=0)
        self.rec_low_status=StringVar()
        self.label_rec_low_status=Label(top,textvariable=self.rec_low_status)
        self.label_rec_low_status.grid(row=7,column=1,columnspan=4,sticky=W)
        
        ####################
        
        #### creating the variables and the input boxes in the GUI that has to do with projection data ######
        
        self.button_choose_ROI=Button(top, text='Choose ROI',command=self.show)
        self.button_choose_ROI.grid(row=8,column=0)
        
        self.nr_x_label = Label(top, text="Nr. of subvolumes in x-direction:")
        self.nr_x_label.grid(row=8,column=5,columnspan=3,sticky=W)
        self.nr_x_entry = Entry(top, width=5)
        self.nr_x_entry.grid(row=8,column=8,sticky=W)
        
        self.nr_y_label = Label(top, text="Nr. of subvolumes in y-direction:")
        self.nr_y_label.grid(row=9,column=5,columnspan=3,sticky=W)
        self.nr_y_entry = Entry(top, width=5)
        self.nr_y_entry.grid(row=9,column=8,sticky=W)
        
        self.nr_z_label = Label(top, text="Nr. of subvolumes in z-direction:")
        self.nr_z_label.grid(row=10,column=5,columnspan=3,sticky=W)
        self.nr_z_entry = Entry(top, width=5)
        self.nr_z_entry.grid(row=10,column=8,sticky=W)
        ###############
        #### creating the variables and the input boxes in the GUI that has to do the multi resolution technique ######
        
        self.method_draw=False
        self.method_parts=False
        self.var = StringVar()
        self.radio_draw = Radiobutton(top, text="Drawn ROI", variable=self.var, value='drawn',command=self.ROI_method)
        self.radio_draw.grid( row=11,column=0 )
        self.radio_parts = Radiobutton(top, text="Choose parts", variable=self.var, value='parts',command=self.ROI_method)
        self.radio_parts.grid( row=11,column=1 )
        self.radio_parts_status=StringVar()
        self.label_radio_parts_status=Label(top,textvariable=self.radio_parts_status)
        self.label_radio_parts_status.grid(row=11,column=2,columnspan=4,sticky=W)
        
        self.algo_high = StringVar()
        self.algobox_high = ttk.Combobox(top, values=('SIRT', 'SART', 'CGLS'),state='readonly')
        self.algobox_high.bind("<<ComboboxSelected>>", self.set_comboparam)
        self.algobox_high.grid(row=21,column=1,columnspan=3,sticky=W)
        self.algo_high_label = Label(top, text="Algorithm:")
        self.algo_high_label.grid(row=21,column=0)
        
        self.nr_high_label = Label(top, text="Nr. of iterations high resolution:")
        self.nr_high_label.grid(row=22,column=0)
        self.nr_high_entry = Entry(top, width=5)
        self.nr_high_entry.grid(row=22,column=1,sticky=W)
        
        self.button_rec_high=Button(top, text='Start high resolution reconstruction')
        self.button_rec_high.bind("<ButtonPress-1>",lambda event: self.proj_status_text(event,'rec_high'))
        self.button_rec_high.bind("<ButtonRelease-1>",self.start_high_rec)
        self.button_rec_high.grid(row=23,column=0)
        self.rec_high_status=StringVar()
        self.label_rec_high_status=Label(top,textvariable=self.rec_high_status)
        self.label_rec_high_status.grid(row=23,column=1,columnspan=4,sticky=W)

        ###############

    def start_high_rec(self,event):
        if len(self.var.get())==0:
            self.radio_parts_status.set('Please choose method')
            self.rec_high_status.set('Something went wrong, please check your input')  
        else:
            if self.overlap>0:   #change to try.
                overlap=self.overlap
                if self.var.get()=='drawn':
                    ##### if the draw method was chosen new boundariy matrices has to be created
                    z_overlap=np.round(overlap*(self.last_z-self.start_z))
                    self.z_subvol_width=self.last_z-self.start_z
                    start_minus_overlap_z=self.start_z-z_overlap
                    start_plus_overlap_z=self.start_z+z_overlap
                    last_minus_overlap_z=self.last_z-z_overlap
                    last_plus_overlap_z=self.last_z+z_overlap
                    if start_minus_overlap_z < 0:
                        start_minus_overlap_z=0
                    if start_plus_overlap_z > self.vol_size:
                        start_plus_overlap_z=self.vol_size
                    if last_minus_overlap_z < 0:
                        last_minus_overlap_z=0
                    if last_plus_overlap_z > self.vol_size:
                        last_plus_overlap_z=self.vol_size                                                  
                    self.z_sub_boundaries=np.array([[self.start_z,start_minus_overlap_z,start_plus_overlap_z], [self.last_z, last_minus_overlap_z, last_plus_overlap_z]])
                    nr_subvol_z=len(self.z_sub_boundaries)-1
                    z_parts=np.linspace(0,0,1) # will only have one part in each direction
            
                    y_overlap=np.round(overlap*(self.last_y-self.start_y))
                    self.y_subvol_width=self.last_y-self.start_y 
                    start_minus_overlap_y=self.start_y-y_overlap
                    start_plus_overlap_y=self.start_y+y_overlap
                    last_minus_overlap_y=self.last_y-y_overlap
                    last_plus_overlap_y=self.last_y+y_overlap
                    if start_minus_overlap_y < 0:
                        start_minus_overlap_y=0
                    if start_plus_overlap_y > self.vol_size:
                        start_plus_overlap_y=self.vol_size
                    if last_minus_overlap_y < 0:
                        last_minus_overlap_y=0
                    if last_plus_overlap_y > self.vol_size:
                        last_plus_overlap_y=self.vol_size                                                  
                    self.y_sub_boundaries=np.array([[self.start_y,start_minus_overlap_y,start_plus_overlap_y], [self.last_y, last_minus_overlap_y, last_plus_overlap_y]])
                    nr_subvol_y=len(self.y_sub_boundaries)-1
                    y_parts=np.linspace(0,0,1) # will only have one part in each direction
            
                    x_overlap=np.round(overlap*(self.last_x-self.start_x))
                    self.x_subvol_width=self.last_x-self.start_x
                    start_minus_overlap_x=self.start_x-x_overlap
                    start_plus_overlap_x=self.start_x+x_overlap
                    last_minus_overlap_x=self.last_x-x_overlap
                    last_plus_overlap_x=self.last_x+x_overlap
                    if start_minus_overlap_x < 0:
                        start_minus_overlap_x=0
                    if start_plus_overlap_x > self.vol_size:
                        start_plus_overlap_x=self.vol_size
                    if last_minus_overlap_x < 0:
                        last_minus_overlap_x=0
                    if last_plus_overlap_x > self.vol_size:
                        last_plus_overlap_x=self.vol_size                                                  
                    self.x_sub_boundaries=np.array([[self.start_x,start_minus_overlap_x,start_plus_overlap_x], [self.last_x, last_minus_overlap_x, last_plus_overlap_x]])
                    nr_subvol_x=len(self.x_sub_boundaries)-1
                    x_parts=np.linspace(0,0,1) # will only have one part in each direction
                    ######
                else:
                    x_parts,y_parts,z_parts=self.parts_list()

                rec_low_temp=np.zeros((self.vol_size_low,self.vol_size_low,self.vol_size_low),dtype= np.float16)


                ##################################HERE Multi resolution starts
                ii=0

                nr_iterations_high=int(self.nr_high_entry.get())
                vol_geom_low = astra.create_vol_geom(self.vol_size_low, self.vol_size_low, self.vol_size_low)
                #creating the CT geometry
                proj_geom_low = astra.create_proj_geom('cone',  self.det_spacing_x_low, self.det_spacing_y_low, self.det_row_count, self.det_col_count, self.angles, self.source_origin, self.origin_det)
                algon_high=self.algobox_high.get()

                if str(algon_high)=='SART':
                    algon_high='EGEN_SART'

                for n_x in x_parts:
                    for n_y in y_parts:
                        for n_z in z_parts: ## nr_subvol_z
                            rec_low_temp[:,:,:]=self.rec_low[:,:,:]
                            #put subvolume to zero
                            rec_low_temp[self.x_sub_boundaries[n_x,1] / self.vol_size_ratio:self.x_sub_boundaries[n_x+1,2] / self.vol_size_ratio, self.y_sub_boundaries[n_y,1] / self.vol_size_ratio:self.y_sub_boundaries[n_y+1,2] / self.vol_size_ratio, self.z_sub_boundaries[n_z,1] / self.vol_size_ratio:self.z_sub_boundaries[n_z+1,2] / self.vol_size_ratio] = 0
                            # make projection data when subvolume is put to zero.
                            proj_id_low, proj_data_low = astra.create_sino3d_gpu(rec_low_temp, proj_geom_low, vol_geom_low)
                            astra.data3d.delete(proj_id_low)
                            ### When splitting in the axis of rotation the whole projection matrix is no needed, here we calculate what part that is actually needed
                            x_start=self.x_sub_boundaries[n_x,1]
                            x_stop=self.x_sub_boundaries[n_x+1,2]
                            h_start=self.FOV/2 - x_start
                            h_stop=self.FOV/2 - x_stop  #height
                            if x_start==0:
                                sin_angle_start= float(self.FOV/2 - x_start + self.FOV/2*(self.det/2)/(self.origin_det + self.source_origin))/(self.source_origin - self.FOV/2) #not needed
                                pro_start=0
                            else:
                                if h_start < 0:
                                    sin_angle_start= float(self.FOV/2 - x_start)/(self.source_origin + self.FOV/2)
                                else:
                                    sin_angle_start= float(self.FOV/2 - x_start)/(self.source_origin - self.FOV/2)
                                pro_start=np.floor(self.det/2 -(self.source_origin+self.origin_det)*sin_angle_start)
                            if x_stop == self.vol_size:
                                sin_angle_stop=float(self.FOV/2 - x_stop - self.FOV/2*(self.det/2)/(self.origin_det + self.source_origin))/(self.source_origin - self.FOV/2)
                                pro_stop=self.det
                            else:
                                if h_stop < 0:
                                    sin_angle_stop= float(self.FOV/2 - x_stop)/(self.source_origin - self.FOV/2)
                                else:
                                    sin_angle_stop= float(self.FOV/2 - x_stop)/(self.source_origin + self.FOV/2)
                                pro_stop=np.ceil(self.det/2 -(self.source_origin+self.origin_det)*sin_angle_stop)
                            self.det_move=self.det/2 - pro_stop + float(pro_stop -pro_start)/2
                            self.det_move=-self.det_move
                            proj_data_part=self.proj_data[pro_start:pro_stop,:,:]
                            proj_data_low_part=proj_data_low[pro_start:pro_stop,:,:]
                            #########
                            
                            relaxation=0.05 # relaxation for the self-made SART aalgorithm, this value  needs to be changed though.
                            rec_sub_vol =recon_multi_part(algon_high,proj_data_part,proj_data_low_part,self.det_move,self.angles,self.x_sub_boundaries,self.y_sub_boundaries,self.z_sub_boundaries,self.vol_size,self.det_spacing_x,self.det_spacing_y,self.x_subvol_width,self.y_subvol_width,self.z_subvol_width,n_x,n_y,n_z,nr_iterations_high,self.source_origin, self.origin_det,relaxation)
                            sub_shape=np.shape(rec_sub_vol)
                            rec_sub_roi= rec_sub_vol[(self.x_sub_boundaries[n_x,0]-self.x_sub_boundaries[n_x,1]):sub_shape[0]-(self.x_sub_boundaries[n_x+1,2]-self.x_sub_boundaries[n_x+1,0]),(self.y_sub_boundaries[n_y,0]-self.y_sub_boundaries[n_y,1]):sub_shape[1]-(self.y_sub_boundaries[n_y+1,2]-self.y_sub_boundaries[n_y+1,0]),(self.z_sub_boundaries[n_z,0]-self.z_sub_boundaries[n_z,1]):sub_shape[2]-(self.z_sub_boundaries[n_z+1,2]-self.z_sub_boundaries[n_z+1,0])]
                            #rec_sub_roi= rec_sub_vol[(self.x_sub_boundaries[n_x,0]-self.x_sub_boundaries[n_x,1]):(self.x_sub_boundaries[n_x,0]-self.x_sub_boundaries[n_x,1]+self.x_subvol_width),(self.y_sub_boundaries[n_y,0]-self.y_sub_boundaries[n_y,1]):(self.y_sub_boundaries[n_y,0]-self.y_sub_boundaries[n_y,1]+self.y_subvol_width),(self.z_sub_boundaries[n_z,0]-self.z_sub_boundaries[n_z,1]):(self.z_sub_boundaries[n_z,0]-self.z_sub_boundaries[n_z,1]+self.z_subvol_width)]
                            
                            print ii
                            ii=ii+1


                            outfile='C:/Users/Sebastian/Documents/CTdata/Saved Data/rec_ROI_'+str(ii)+'det_'+str(self.det)
                            sio.savemat(outfile+'.mat', {'rec_sub_roi':rec_sub_roi})
                            save_as_low='part'+ str(self.x_sub_boundaries[n_x,0]) + str(self.y_sub_boundaries[n_y,0] ) +str(self.z_sub_boundaries[n_z,0])
                            affine=np.eye(4)# affine transformation matrix when nothing is moved or scaled.
                            affine[0,3]=self.x_sub_boundaries[n_x,0] #translation in x-direction 
                            affine[1,3]=self.y_sub_boundaries[n_y,0] #translation in y-direction
                            affine[2,3]=self.z_sub_boundaries[n_z,0] #translation in z-direction
                            array_img=nib.Nifti1Image(rec_sub_roi, affine)
                            nib.save(array_img, save_as_low +'_full.nii')
                            #np.savez(outfile, rec_ROI=rec_ROI, width=width)
                self.rec_high_status.set('High resolution reconstruction done') 
                ### this part can be removed, better to watch the reconstructions in paraview.            
                full_rec=read_and_show(x_parts,y_parts,z_parts,self.x_sub_boundaries,self.y_sub_boundaries,self.z_sub_boundaries,self.det)

                pylab.figure()
                x_rec=0
                for i in x_parts:
                    x_rec=x_rec + self.x_sub_boundaries[i+1,0] - self.x_sub_boundaries[i,0]   
                y_rec=0
                for i in y_parts:
                    y_rec=y_rec + self.y_sub_boundaries[i+1,0] - self.y_sub_boundaries[i,0]  
                z_rec=0
                for i in z_parts:
                    z_rec=z_rec + self.z_sub_boundaries[i+1,0] - self.z_sub_boundaries[i,0] 
                pylab.imshow(full_rec[0.5*x_rec,:,:])
                pylab.figure()
                pylab.imshow(full_rec[:,0.5*y_rec,:])
                pylab.show()
                #####
                 
            else: # change to except:
                self.rec_high_status.set('Something went wrong, please check your input')  
                 
    def ROI_method(self):
        # creates the input box to choose method to create the subvolume(s)
        
        if self.var.get()=='draw':
            self.radio_parts_status.set(' ')
            self.button_choose_ROI=Button(top, text='Choose ROI',command=self.show)
            self.button_choose_ROI.grid(row=8,column=0)
            self.method_draw=True
            self.method_parts=False
        elif self.var.get()=='parts':
            self.radio_parts_status.set(' ')
            # create the boundary matrices which are used in the multi resolution technique
            self.x_sub_boundaries,self.y_sub_boundaries,self.z_sub_boundaries,self.x_subvol_width,self.y_subvol_width,self.z_subvol_width=create_boundaries(int(self.nr_x_entry.get()),int(self.nr_y_entry.get()),int(self.nr_z_entry.get()),self.vol_size,self.overlap)
            
            if self.method_draw==True:
                self.button_choose_ROI.destroy()
                self.method_draw=False
            if self.method_parts==True:
                self.label_x_bound.destroy()
                self.label_y_bound.destroy()
                self.label_z_bound.destroy()
                    
            #####BOXES for input#########    
            self.parts_x_label = Label(top, text='Part(s) in x-direction')
            self.parts_x_label.grid(row=18,column=0)
            self.parts_x_entry = Entry(top, width=5)
            self.parts_x_entry.grid(row=18,column=1,sticky=E)
            self.parts_streck_label = Label(top, text='-')
            self.parts_streck_label.grid(row=18,column=2)
            self.parts_x_entry2 = Entry(top, width=5)
            self.parts_x_entry2.grid(row=18,column=3,sticky=W)
            
            self.parts_y_label = Label(top, text='Part(s) in y-direction')
            self.parts_y_label.grid(row=19,column=0)
            self.parts_y_entry = Entry(top, width=5)
            self.parts_y_entry.grid(row=19,column=1,sticky=E)
            self.parts_streck_label = Label(top, text='-')
            self.parts_streck_label.grid(row=19,column=2)
            self.parts_y_entry2 = Entry(top, width=5)
            self.parts_y_entry2.grid(row=19,column=3,sticky=W)
            
            self.parts_z_label = Label(top, text='Part(s) in z-direction')
            self.parts_z_label.grid(row=20,column=0)
            self.parts_z_entry = Entry(top, width=5)
            self.parts_z_entry.grid(row=20,column=1,sticky=E)
            self.parts_streck_label = Label(top, text='-')
            self.parts_streck_label.grid(row=20,column=2)
            self.parts_z_entry2 = Entry(top, width=5)
            self.parts_z_entry2.grid(row=20,column=3,sticky=W)
            ###############################
            
            #####showing the parts##########
            self.label_head_bound=Label(top,text='Part:')
            self.label_head_bound.grid(row=17,column=4)
            if self.x_sub_boundaries.shape[0] >= self.y_sub_boundaries.shape[0]:
                if self.x_sub_boundaries.shape[0] >= self.z_sub_boundaries.shape[0]:
                    max_parts=self.x_sub_boundaries.shape[0]
                else:
                    max_parts=self.z_sub_boundaries.shape[0]
            elif self.z_sub_boundaries.shape[0] >= self.y_sub_boundaries.shape[0]:
                max_parts=self.z_sub_boundaries.shape[0]
            else:
                max_parts=self.y_sub_boundaries.shape[0]
                
            for ii in range(0,max_parts-1):
                self.label_head_bound=Label(top,text=str(ii))
                self.label_head_bound.grid(row=17,column=ii+5) 
            
            for ii in range(0,self.x_sub_boundaries.shape[0]-1):
                self.label_x_bound=Label(top,text=str(self.x_sub_boundaries[ii,0])+'-'+str(self.x_sub_boundaries[ii+1,0]))
                self.label_x_bound.grid(row=18,column=ii+5)
            for ii in range(0,self.y_sub_boundaries.shape[0]-1):
                self.label_y_bound=Label(top,text=str(self.y_sub_boundaries[ii,0])+'-'+str(self.y_sub_boundaries[ii+1,0]))
                self.label_y_bound.grid(row=19,column=ii+5)
            for ii in range(0,self.z_sub_boundaries.shape[0]-1):
                self.label_z_bound=Label(top,text=str(self.z_sub_boundaries[ii,0])+'-'+str(self.z_sub_boundaries[ii+1,0]))
                self.label_z_bound.grid(row=20,column=ii+5)
            self.method_parts=True
    def parts_list(self):
        #creating a list on what region every part equals in the reconstruction volume
        nr_subvol_x=int(self.nr_x_entry.get())
        nr_subvol_y=int(self.nr_y_entry.get())
        nr_subvol_z=int(self.nr_z_entry.get())
        if nr_subvol_x >1:
            if len(self.parts_x_entry2.get())==0 and int(self.parts_x_entry.get()) < nr_subvol_x:
                x_parts=np.linspace(int(self.parts_x_entry.get()),int(self.parts_x_entry.get()),1)
            elif len(self.parts_x_entry2.get())>0  and  int(len(self.parts_x_entry2.get())) < nr_subvol_x:
                first_x=int(self.parts_x_entry.get())
                last_x=int(self.parts_x_entry2.get())
                x_parts=np.linspace(first_x,last_x,last_x-first_x+1)
        else:
            x_parts=np.linspace(0,0,1)
    
        print x_parts 
        
        if nr_subvol_y >1:
            if len(self.parts_y_entry2.get())==0 and int(self.parts_y_entry.get()) < nr_subvol_y:
                y_parts=np.linspace(int(self.parts_y_entry.get()),int(self.parts_y_entry.get()),1)
            elif len(self.parts_y_entry2.get())>0  and  int(len(self.parts_y_entry2.get())) < nr_subvol_y:
                first_y=int(self.parts_y_entry.get())
                last_y=int(self.parts_y_entry2.get())
                y_parts=np.linspace(first_y,last_y,last_y-first_y+1)
        else:
            y_parts=np.linspace(0,0,1)
    
        print y_parts 
        
        if nr_subvol_z >1:
            if len(self.parts_z_entry2.get())==0 and int(self.parts_z_entry.get()) < nr_subvol_z:
                z_parts=np.linspace(int(self.parts_z_entry.get()),int(self.parts_z_entry.get()),1)
            elif len(self.parts_z_entry2.get())>0  and  int(len(self.parts_z_entry2.get())) < nr_subvol_z:
                first_z=int(self.parts_z_entry.get())
                last_z=int(self.parts_z_entry2.get())
                z_parts=np.linspace(first_z,last_z,last_z-first_z+1)
        else:
            z_parts=np.linspace(0,0,1)
    
        print z_parts 

        
        return x_parts,y_parts, z_parts   

    def start_low_rec(self,event):
        if 1==1: #try: #  better to use try and show the error for the user when the program is fully tested
            self.vol_size=self.FOV

            self.vol_size_low=self.vol_size/self.vol_size_ratio # the size of each side in the low rec. volume.

            self.angles = np.linspace(0, 2*np.pi, self.number_of_projections,False)
            nr_iterations_low=int(self.nr_low_entry.get())
            # self.det_spacing_x=np.round(float(self.vol_size)/new_size[0],decimals=2) #change here if you want maginification of rec volume.
            # self.det_spacing_y=np.round(float(self.vol_size)/new_size[0],decimals=2)
            self.det_spacing_x=1
            self.det_spacing_y=1
            self.det_row_count=self.new_size[0]
            self.det_col_count=self.new_size[0]

            self.det_spacing_x_low=np.round(float(self.det_spacing_x)/self.vol_size_ratio,decimals=2) #detector spacing compared to voxel size
            self.det_spacing_y_low=np.round(float(self.det_spacing_y)/self.vol_size_ratio,decimals=2)

            algon=self.algobox.get()
            relaxation=0.05 # the relaxation parameter, should maybe be a input box instead
            initializer=0 # start reconstruction with a volume of zeros
            rec_volume=(self.vol_size_low,self.vol_size_low,self.vol_size_low)
            self.det_size=(self.det,self.det_low)
            # the reconstruction function is called.
            if algon=='FDK':
                algo='FDK_CUDA'
                nr_iterations_low=1
                self.rec_low=make_reconsruction(self.proj_data_low,self.number_of_projections,self.vol_size_ratio,self.det_size,rec_volume,self.source_origin,self.origin_det,nr_iterations_low,algo,relaxation,initializer)
                self.rec_low=self.rec_low*(self.rec_low>0)
            elif algon=='SIRT':
                algo='SIRT3D_CUDA'
                self.rec_low=make_reconsruction(self.proj_data_low,self.number_of_projections,self.vol_size_ratio,self.det_size,rec_volume,self.source_origin,self.origin_det,nr_iterations_low,algo,relaxation,initializer)
            elif algon=='CGLS':
                algo='CGLS3D_CUDA'
                self.rec_low=make_reconsruction(self.proj_data_low,self.number_of_projections,self.vol_size_ratio,self.det_size,rec_volume,self.source_origin,self.origin_det,nr_iterations_low,algo,relaxation,initializer)
            elif algon=='SART': # can use the whole projection data
                algo='EGEN_SART'
                rec_low=make_reconsruction(self.proj_data,self.number_of_projections,self.vol_size_ratio,self.det_size,rec_volume,self.source_origin,self.origin_det,nr_iterations_low,algo,relaxation,initializer)
            elif algon=='Landweber': # can use the whole projection data
                algo='EGEN_SIRT'
                rec_low=make_reconsruction(self.proj_data_low,self.number_of_projections,self.vol_size_ratio,self.det_size,rec_volume,self.source_origin,self.origin_det,nr_iterations_low,algo,relaxation,initializer)
            elif algon=='FDK_SIRT':
                algo='FDK_CUDA'
                self.rec_low=make_reconsruction(self.proj_data_low,self.number_of_projections,self.vol_size_ratio,self.det_size,rec_volume,self.source_origin,self.origin_det,nr_iterations_low,algo,relaxation,initializer)
                self.rec_low=self.rec_low*(self.rec_low>0)
                algo='SIRT3D_CUDA'
                self.rec_low=make_reconsruction(self.proj_data_low,self.number_of_projections,self.vol_size_ratio,self.det_size,rec_volume,self.source_origin,self.origin_det,nr_iterations_low,algo,relaxation,self.rec_low)
            elif algon=='FDK_CGLS':
                algo='FDK_CUDA'
                self.rec_low=make_reconsruction(self.proj_data_low,self.number_of_projections,self.vol_size_ratio,self.det_size,rec_volume,self.source_origin,self.origin_det,nr_iterations_low,algo,relaxation,initializer)
                self.rec_low=self.rec_low*(self.rec_low>0)
                algo='CGLS3D_CUDA'
                self.rec_low=make_reconsruction(self.proj_data_low,self.number_of_projections,self.vol_size_ratio,self.det_size,rec_volume,self.source_origin,self.origin_det,nr_iterations_low,algo,relaxation,self.rec_low)
            else:
                print'Algorithm not found, following algorithms exists: FDK,SIRT,CGLS,FDK_SIRT and FDK_CGLS'
                quit()

            proj_data_low=[]
            ###Get the projection image with parallel rays so it can be used to for drawing region of interest
            vol_geom = astra.create_vol_geom(rec_volume)
            det_spacing=np.round((float(self.det_size[0])/self.det_size[1])/self.vol_size_ratio,decimals=3)
            proj_geom = astra.create_proj_geom('parallel3d', 1, 1, self.vol_size_low, self.vol_size_low, self.angles)
            sinogram_id, proj_parallel = astra.create_sino3d_gpu(self.rec_low, proj_geom,vol_geom)
            self.par_proj_im=proj_parallel[:,0,:] # only one angle is needed for the draw.
            vol_geom=[]
            proj_parallel=[]
            astra.data3d.delete(sinogram_id)
            ####
            
            save_as_low='test_it'+str(nr_iterations_low)+'_'+algon # the name of the nifti file that is saved.
            affine=np.eye(4) # affine transformation matrix when nothing is moved.
            array_img=nib.Nifti1Image(self.rec_low, affine)
            nib.save(array_img, save_as_low+'_low.nii')
            self.rec_low_status.set('Low resolution reconstruction done') 
        else:#except:
            self.rec_low_status.set('Something went wrong, please check your input') 
        
    def proj_status_text(self,event,invar):
        if invar=='process':
            self.proj_status.set('Processing data, please wait')
        elif invar=='rec_low':
            self.rec_low_status.set('Low resolution reconstruction started, please wait')  
        elif invar=='rec_high':
            self.rec_high_status.set('High resolution reconstruction started, please wait')       
    def get_proj(self,event):
        #self.proj_status.set('Processing data, please wait')
        if 1==1: #try: # use try when you want to show error in GUI
            self.det=int(self.nr_det_entry.get())
            self.vol_size_ratio=int(self.vol_ratiobox.get())
            self.ang_space=1
            if self.det <= 750: 
                self.det_low=self.det
            else: #uses a smaller detector size for low resolution reconstruction to save memory and time, the threshold depends on the computer
                self.det_low=750
                    # vol_size_ratio=2
                    # ang_space=1
            self.new_size=(int(self.det),int(self.det)) # new size is the number of detector in every direction. if 2000 it's equal to do nothing.
            print self.dirname_data.get()
            filenames=glob.glob(self.dirname_data.get()+'/*.bin') # makes a list of filenames to the projection data
            print filenames
            beam_profile_file=self.dirname_lf.get()+'/*.bin' # the place of beam profile files
            dark_current_file=self.dirname_dc.get()+'/*.bin' # the place of beam profile files
            
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
            self.origin_det=np.round(origin_to_detector_pixels/np.round(pixels/self.new_size[1]))
            self.source_origin=np.round(source_to_origin_pixels/np.round(pixels/self.new_size[1]))
                        
            alignment_in_mm=param[2:4] #dx, dz translation, third input is rotation, this part should be read from a txt file instead.
            ###
            
            filenames.sort() # sorts the filenames so the come in the order they are obtained
            number_of_files= len(filenames)
            print number_of_files
            self.number_of_projections=number_of_files/self.ang_space
            print self.number_of_projections, "number of projections "
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
            #reads the bin files and makes the transformation attenuation values
            self.proj_data,self.proj_data_low= making_sino_from_bin(filenames,self.lag_corr.get(),number_of_files, self.number_of_projections,self.ang_space, pixels,beam_profile_file,dark_current_file,self.new_size,self.det_low,pixels_to_crop,alignment_in_mm,a,b)
            

            # calculates the FOV so only this part is needed for the reconstruction
            self.FOV= 2*int((self.det/2)/(self.origin_det + self.source_origin) * self.source_origin ) 
            
            self.proj_status.set('Data processed!')
        else: #except: #use except to show error in GUI
            self.proj_status.set('Something went wrong, please try again!')

        
    def askdirectory(self,button_name):
        #gets the directory name
        if button_name=='proj':
            self.dirname_data.set(tkFileDialog.askdirectory())
        elif button_name=='beam':
            self.dirname_lf.set(tkFileDialog.askdirectory())
        elif button_name=='dark':
            self.dirname_dc.set(tkFileDialog.askdirectory())
        
                
    def set_comboparam(self,event):
        self.algo.set(self.algobox.get())

        
    def show(self):

        #P is resized to the same size as the x- plane in the reconstruction volume
        P=misc.imresize(self.par_proj_im, (self.vol_size,self.vol_size))

        
        max_val=float(np.max(P))
        phan=Image.fromarray(P/(max_val/255)) #to get values from 0-255, so they are properly displayed
        ROI=False #so the draw class knows that Region of interest is not drawn yet.
        app = draw(phan,ROI,self.vol_size,self.vol_size_ratio,self.rec_low)
        
        app.mainloop()
        
        # getting the boundaries for the region of interest
        self.start_z=app.rec_start_z
        self.last_z=app.rec_last_z
        self.start_y=app.rec_start_y
        self.last_y=app.rec_last_y
        self.start_x=app.rec_start_x
        self.last_x=app.rec_last_x
        
        print app.rec_start_x,app.rec_last_x,app.rec_start_y, app.rec_last_y,app.rec_start_z, app.rec_last_z
        
#creating a Tkinter object        
top = Tk() 
top.title('Reconstruction')
top.geometry('1000x700')

a=reconstruction_GUI(top)


top.mainloop()