#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 17 14:49:34 2018

@author: brh
"""

class FileReadWrite:
    
#==============================================================================
    def absolute_path(self, output_directory, working_directory):
        '''
        create full path for local output directory with error checking
        '''
        import os
        import sys
        
        if os.path.isdir(output_directory) != True:
            output_directory = os.path.join \
                (working_directory, output_directory)
            if os.path.isdir(output_directory) != True:
                os.makedirs(output_directory)
                print('warning: NOT FOUND output_directory of ', \
                      output_directory, '; creating directory')
                #sys.exit()
            #endif
        #endif
        
        return output_directory
    
#==============================================================================
    def add_timestamp_to_filename(self, thisfile, timestamp):
        '''
        Adds a time stamp to a file name
        If the filename has a . extention, the time stamp is before the .
        '''
        import datetime
        
        if timestamp == []:
            timestamp = datetime.datetime.now().strftime('%Y%m%d_%H%M')
        #endif
        
        temp = thisfile.rsplit('.')
        fout = temp[0]+'_'+timestamp
        
        if len(temp) > 1:
            fout = fout+'.'+temp[1]
        #endif
        
        return fout
    
#==============================================================================
    def get_oneline_data_matching_time(self,time_to_match, openfile, isnoisy):
        '''
        Looks through file to find time that is equal to or less than the 
        time_to_match
        '''
        
        import sys
        import numpy as np
        
        # HACK - the error checking on this needs work.  If it doesn't find
        # a time it seems to start with zero values instead of an error
        
        usenextdata = False
        dontquit = True
        dataout = None
        idxL = 0

        #----------------------------------------------------------------------                       
        while dontquit == True:
            rln = openfile.readline()
            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            if len(rln) == 0:
                if isnoisy == True:
                    print('end of file')
                #endif    
                break
            #endif
            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            rln = rln.rstrip('\n')
            pln = rln.partition(" ")
            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            if pln[0] == 'HEADER:':
                if isnoisy == True:
                    print(pln[2])
                #endif    
            elif pln[0] == 'FILETYPE:':
                if isnoisy == True:
                    print(pln[2])
                #endif    
            elif pln[0] == 'DATE:':
                if isnoisy == True:
                    print(pln[2])
                #endif    
            elif pln[0] == 'DATATYPE:':
                if isnoisy == True:
                    print(pln[2])
                #endif    
            elif pln[0] == 'UNITS:':
                if isnoisy == True:
                    print(pln[2])
                #endif    
            elif pln[0] == 'NX:':
                if isnoisy == True:
                    print(pln[2])
                #endif    
                NXlocal = int(pln[2])               
                dataout = np.zeros(NXlocal,dtype=np.float64)
                if isnoisy == True:
                    print(pln[2])
                #endif    
            elif pln[0] == 'END:':
                if isnoisy == True:
                    print(pln[2])
                #endif    
            elif pln[0] == 'ITEMS_PER_LINE:':
                if isnoisy == True:
                    print(pln[2])                                           
                #endif    
            elif pln[0] == 'TIME:':
                if usenextdata == True:
                    #sys.exit()
                    usenextdata = False
                    dontquit = False
                #endif    
                if isnoisy == True:
                    print(pln[2])
                #endif    
                linetime = float(pln[2])
                if linetime >= time_to_match:
                    usenextdata = True                       
                #endif    
            elif pln[0] == 'DATA:':
                if isnoisy == True:
                    print(pln[2])
                #endif    
            else:        
                if usenextdata == True:
                    rln = rln.lstrip(' ')
                    rln = rln.rstrip('\n')
                    rln = rln.rstrip(' ')
                    rln = rln.rstrip(']')
                    rln = rln.lstrip('[')
                    if isnoisy == True:
                        print('reading data')
                        print(rln)
                    #endif    
                    pln = rln.split()
                    dtemp = [float(s) for s in pln]
                    dataout[idxL:idxL + len(dtemp)] = dtemp
                    idxL = idxL + len(dtemp)
                #endif
            #endif
            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        #endwhile
        #----------------------------------------------------------------------                       
       
        return [linetime, dataout]

#==============================================================================
    def geometry_txout_write(self, geo, fa, setting, trun, NX):
        '''
        writes the text files of all the geometry
        '''
        import numpy as np
        import math
        
        if setting.txtout_writedata == True:

            setting.txtout_file_zbottom.write \
                (''.join(['TIME: ',str(trun.present_time),'\n']))
            setting.txtout_file_zbottom.write('DATA: NEXTLINE\n')

            setting.txtout_file_zbottom_face.write \
                (''.join(['TIME: ',str(trun.present_time),'\n']))
            setting.txtout_file_zbottom_face.write('DATA: NEXTLINE\n')
           
            setting.txtout_file_length.write \
                (''.join(['TIME: ',str(trun.present_time),'\n']))
            setting.txtout_file_length.write('DATA: NEXTLINE\n')
    
            setting.txtout_file_xvalue.write \
                (''.join(['TIME: ',str(trun.present_time),'\n']))
            setting.txtout_file_xvalue.write('DATA: NEXTLINE\n')
           
            setting.txtout_file_xvalue_face.write \
                (''.join(['TIME: ',str(trun.present_time),'\n']))
            setting.txtout_file_xvalue_face.write('DATA: NEXTLINE\n')
           
            setting.txtout_file_manningsn.write \
                (''.join(['TIME: ',str(trun.present_time),'\n']))
            setting.txtout_file_manningsn.write('DATA: NEXTLINE\n')
           
            # number of lines to be written. This is needed because
            # a large system might need to be written over multiple lines
            nline = math.ceil(NX / setting.txtout_items_on_line_maximum)
            
            idxL = 0
            idxH = setting.txtout_items_on_line_maximum
            
            #print('here ',nline, NX, setting.txtout_items_on_line_maximum)
            
            for ii in range(0,nline):            
                ts1 = geo['zbottom'][idxL:idxH]
                ts2 = geo['length'][idxL:idxH]
                ts3 = geo['xvalue'][idxL:idxH]
                ts4 = geo['manningsn'][idxL:idxH]
                ts5 = geo['breadth'][idxL:idxH]
                ts6 = fa['zbottom'][idxL:idxH]
                ts7 = fa['xvalue'][idxL:idxH]
                
                ts1 = np.array_str(ts1,precision = 16)
                ts2 = np.array_str(ts2,precision = 16)
                ts3 = np.array_str(ts3,precision = 16)
                ts4 = np.array_str(ts4,precision = 16)
                ts5 = np.array_str(ts5,precision = 16)
                ts6 = np.array_str(ts6,precision = 16)
                ts7 = np.array_str(ts7,precision = 16)
                
                setting.txtout_file_zbottom.write(ts1)
                setting.txtout_file_zbottom.write('\n')
                
                setting.txtout_file_length.write(ts2)
                setting.txtout_file_length.write('\n')
                
                setting.txtout_file_xvalue.write(ts3)
                setting.txtout_file_xvalue.write('\n')
                
                setting.txtout_file_manningsn.write(ts4)
                setting.txtout_file_manningsn.write('\n')
                
                setting.txtout_file_breadth.write(ts5)
                setting.txtout_file_breadth.write('\n')

                setting.txtout_file_zbottom_face.write(ts6)
                setting.txtout_file_zbottom_face.write('\n')

                setting.txtout_file_xvalue_face.write(ts7)
                setting.txtout_file_xvalue_face.write('\n')
                
                idxL = idxH
                idxH = idxL + setting.txtout_items_on_line_maximum 
            #endfor
            
            setting.txtout_file_zbottom.close()
            setting.txtout_file_length.close()
            setting.txtout_file_xvalue.close()
            setting.txtout_file_manningsn.close()
            setting.txtout_file_breadth.close()

            setting.txtout_file_zbottom_face.close()
            setting.txtout_file_xvalue_face.close()
            
        #endif
        
        return
    
#==============================================================================
    def initialize_all_writefile_names(self,setting):
        '''
        Opens and writes initial data to output text files
        These files can be used for restart
        ''' 
        import datetime
        import os
        import sys
        
        # HACK - this should be re-written with a flexible dictionary of 
        # filenames that are iterated over
        
        #Define a time stamp for this run
        setting.thisrun_timestamp  \
            = datetime.datetime.now().strftime('%Y%m%d_%H%M') 
 
        # text file output directory
        setting.txtout_directory = self.absolute_path \
            (setting.txtout_directory, setting.working_directory)   
         
        # filename for saving all simulation settings.    
        setting.setting_filename = self.add_timestamp_to_filename \
            (setting.setting_filename, setting.thisrun_timestamp)
            
        setting.setting_filename = os.path.join \
            (setting.txtout_directory, setting.setting_filename) 
     
        #----------------------------------------------------------------------                       
        # create full path for local output and datasave directories 
        if setting.txtout_writedata == True:
             
            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            # DYNAMIC SOLUTION FILES
            setting.txtout_volume_filename \
                = self.add_timestamp_to_filename \
                (setting.txtout_volume_filename, setting.thisrun_timestamp)
                
            setting.txtout_volume_filename = os.path.join \
                  (setting.txtout_directory, setting.txtout_volume_filename) 
                  
            
            setting.txtout_flowrate_filename \
                = self.add_timestamp_to_filename \
                (setting.txtout_flowrate_filename, setting.thisrun_timestamp)
                
            setting.txtout_flowrate_filename = os.path.join \
                  (setting.txtout_directory, setting.txtout_flowrate_filename) 


            setting.txtout_flowrate_face_filename \
                = self.add_timestamp_to_filename \
                (setting.txtout_flowrate_face_filename, setting.thisrun_timestamp)
                
            setting.txtout_flowrate_face_filename = os.path.join \
                  (setting.txtout_directory, setting.txtout_flowrate_face_filename) 
                  
    
            setting.txtout_hyddepth_filename \
                = self.add_timestamp_to_filename \
                (setting.txtout_hyddepth_filename, setting.thisrun_timestamp)
                
            setting.txtout_hyddepth_filename = os.path.join \
                  (setting.txtout_directory, setting.txtout_hyddepth_filename) 
                  
    
            setting.txtout_elevation_filename \
                = self.add_timestamp_to_filename \
                (setting.txtout_elevation_filename, setting.thisrun_timestamp)
                
            setting.txtout_elevation_filename = os.path.join \
                  (setting.txtout_directory, setting.txtout_elevation_filename) 
                  
                  
            setting.txtout_elevation_faceM_filename \
                = self.add_timestamp_to_filename \
                (setting.txtout_elevation_faceM_filename, setting.thisrun_timestamp)
                
            setting.txtout_elevation_faceM_filename = os.path.join \
                  (setting.txtout_directory, setting.txtout_elevation_faceM_filename) 


            setting.txtout_elevation_faceP_filename \
                = self.add_timestamp_to_filename \
                (setting.txtout_elevation_faceP_filename, setting.thisrun_timestamp)
                
            setting.txtout_elevation_faceP_filename = os.path.join \
                  (setting.txtout_directory, setting.txtout_elevation_faceP_filename) 


            setting.txtout_area_filename \
                = self.add_timestamp_to_filename \
                (setting.txtout_area_filename, setting.thisrun_timestamp)
                
            setting.txtout_area_filename = os.path.join \
                  (setting.txtout_directory, setting.txtout_area_filename) 
                  

            setting.txtout_timestep_filename \
                = self.add_timestamp_to_filename \
                (setting.txtout_timestep_filename, setting.thisrun_timestamp)
                
            setting.txtout_timestep_filename = os.path.join \
                  (setting.txtout_directory, setting.txtout_timestep_filename) 

            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            # GEOMETRY FILES
            setting.txtout_zbottom_filename \
                = self.add_timestamp_to_filename \
                (setting.txtout_zbottom_filename, setting.thisrun_timestamp)
                
            setting.txtout_zbottom_filename = os.path.join \
                  (setting.txtout_directory, setting.txtout_zbottom_filename) 

                  
            setting.txtout_zbottom_face_filename \
                = self.add_timestamp_to_filename \
                (setting.txtout_zbottom_face_filename, setting.thisrun_timestamp)
                
            setting.txtout_zbottom_face_filename = os.path.join \
                  (setting.txtout_directory, setting.txtout_zbottom_face_filename) 
                  
 
            setting.txtout_length_filename \
                = self.add_timestamp_to_filename \
                (setting.txtout_length_filename, setting.thisrun_timestamp)
                
            setting.txtout_length_filename = os.path.join \
                  (setting.txtout_directory, setting.txtout_length_filename) 
                  

            setting.txtout_xvalue_filename \
                = self.add_timestamp_to_filename \
                (setting.txtout_xvalue_filename, setting.thisrun_timestamp)
                
            setting.txtout_xvalue_filename = os.path.join \
                  (setting.txtout_directory, setting.txtout_xvalue_filename) 


            setting.txtout_xvalue_face_filename \
                = self.add_timestamp_to_filename \
                (setting.txtout_xvalue_face_filename, setting.thisrun_timestamp)
                
            setting.txtout_xvalue_face_filename = os.path.join \
                  (setting.txtout_directory, setting.txtout_xvalue_face_filename) 


            setting.txtout_manningsn_filename \
                = self.add_timestamp_to_filename \
                (setting.txtout_manningsn_filename, setting.thisrun_timestamp)
                
            setting.txtout_manningsn_filename = os.path.join \
                  (setting.txtout_directory, setting.txtout_manningsn_filename) 
                  

            setting.txtout_breadth_filename \
                = self.add_timestamp_to_filename \
                (setting.txtout_breadth_filename, setting.thisrun_timestamp)
                
            setting.txtout_breadth_filename = os.path.join \
                  (setting.txtout_directory, setting.txtout_breadth_filename) 
        #endif
        
        #----------------------------------------------------------------------                               
        # BINARY FILES          
        if setting.binsave_writedata == True:
            setting.binsave_directory = self.absolute_path \
                (setting.binsave_directory, setting.working_directory)

            setting.binsave_data_filename \
                =  self.add_timestamp_to_filename \
                (setting.binsave_data_filename, setting.thisrun_timestamp)
                
            setting.binsave_time_filename \
                =  self.add_timestamp_to_filename \
                (setting.binsave_time_filename, setting.thisrun_timestamp)
                                 
            setting.binsave_data_filename = os.path.join \
                 (setting.binsave_directory, setting.binsave_data_filename) 
                 
            setting.binsave_time_filename = os.path.join \
                 (setting.binsave_directory, setting.binsave_time_filename)         
        #endif
             
        return setting
    
#==============================================================================
    def open_txtout_file(self, setting, NX):
        '''
        Opens all the ouptut textfiles
        '''
        import sys
        
        # HACK - this should be rewritten for a more flexible dictionary of 
        # filenames with units, nomenclature, etc.
        
        # open the files for writing text (restart) output
        setting.txtout_file_volume          = open(setting.txtout_volume_filename,'w')
        setting.txtout_file_flowrate        = open(setting.txtout_flowrate_filename,'w')
        setting.txtout_file_flowrate_face   = open(setting.txtout_flowrate_face_filename,'w')
        setting.txtout_file_hyddepth        = open(setting.txtout_hyddepth_filename,'w')
        setting.txtout_file_elevation       = open(setting.txtout_elevation_filename,'w')
        setting.txtout_file_elevation_faceM = open(setting.txtout_elevation_faceM_filename,'w')
        setting.txtout_file_elevation_faceP = open(setting.txtout_elevation_faceP_filename,'w')
        setting.txtout_file_area            = open(setting.txtout_area_filename,'w')
        setting.txtout_file_timestep        = open(setting.txtout_timestep_filename,'w')

        setting.txtout_file_zbottom        = open(setting.txtout_zbottom_filename,'w')
        setting.txtout_file_zbottom_face   = open(setting.txtout_zbottom_face_filename,'w')
        setting.txtout_file_length         = open(setting.txtout_length_filename,'w')
        setting.txtout_file_xvalue         = open(setting.txtout_xvalue_filename,'w')
        setting.txtout_file_xvalue_face    = open(setting.txtout_xvalue_face_filename,'w')
        setting.txtout_file_manningsn      = open(setting.txtout_manningsn_filename,'w')
        setting.txtout_file_breadth        = open(setting.txtout_breadth_filename,'w')
        
        # write the headers
        self.write_header_for_text_output \
            (NX, setting.txtout_items_on_line_maximum, \
             setting.txtout_file_volume, 'volume', 'm^3', \
             setting.txtout_writedata_header, setting.thisrun_timestamp)
            
        self.write_header_for_text_output \
            (NX, setting.txtout_items_on_line_maximum, \
             setting.txtout_file_flowrate, 'flowrate', 'm^2/s',  \
             setting.txtout_writedata_header, setting.thisrun_timestamp)

        self.write_header_for_text_output \
            (NX+1, setting.txtout_items_on_line_maximum, \
             setting.txtout_file_flowrate_face, 'flowrate_face', 'm^2/s',  \
             setting.txtout_writedata_header, setting.thisrun_timestamp)
          
        self.write_header_for_text_output \
            (NX, setting.txtout_items_on_line_maximum, \
             setting.txtout_file_hyddepth, 'hyddepth', 'm',  \
             setting.txtout_writedata_header, setting.thisrun_timestamp)
            
        self.write_header_for_text_output \
            (NX, setting.txtout_items_on_line_maximum, \
             setting.txtout_file_elevation, 'elevation', 'm',  \
             setting.txtout_writedata_header, setting.thisrun_timestamp)

        self.write_header_for_text_output \
            (NX+1, setting.txtout_items_on_line_maximum, \
             setting.txtout_file_elevation_faceM, 'elevation_faceM', 'm',  \
             setting.txtout_writedata_header, setting.thisrun_timestamp)

        self.write_header_for_text_output \
            (NX+1, setting.txtout_items_on_line_maximum, \
             setting.txtout_file_elevation_faceP, 'elevation_faceP', 'm',  \
             setting.txtout_writedata_header, setting.thisrun_timestamp)

        self.write_header_for_text_output \
            (NX, setting.txtout_items_on_line_maximum, \
             setting.txtout_file_area, 'area', 'm^2',  \
             setting.txtout_writedata_header, setting.thisrun_timestamp)

        self.write_header_for_text_output \
            (1, setting.txtout_items_on_line_maximum, \
             setting.txtout_file_timestep, 'timestep', 's',  \
             setting.txtout_writedata_header, setting.thisrun_timestamp)

        self.write_header_for_text_output \
            (NX, setting.txtout_items_on_line_maximum, \
             setting.txtout_file_zbottom, 'zbottom', 'm',  \
             setting.txtout_writedata_header, setting.thisrun_timestamp)

        self.write_header_for_text_output \
            (NX+1, setting.txtout_items_on_line_maximum, \
             setting.txtout_file_zbottom_face, 'zbottom_face', 'm',  \
             setting.txtout_writedata_header, setting.thisrun_timestamp)

        self.write_header_for_text_output \
            (NX, setting.txtout_items_on_line_maximum, \
             setting.txtout_file_length, 'length', 'm',  \
             setting.txtout_writedata_header, setting.thisrun_timestamp)
            
        self.write_header_for_text_output \
            (NX, setting.txtout_items_on_line_maximum, \
             setting.txtout_file_xvalue, 'xvalue', 'm',  \
             setting.txtout_writedata_header, setting.thisrun_timestamp)

        self.write_header_for_text_output \
            (NX+1, setting.txtout_items_on_line_maximum, \
             setting.txtout_file_xvalue_face, 'xvalue_face', 'm',  \
             setting.txtout_writedata_header, setting.thisrun_timestamp)

        self.write_header_for_text_output \
            (NX, setting.txtout_items_on_line_maximum, \
             setting.txtout_file_manningsn, 'manningsn', 'n/a',  \
             setting.txtout_writedata_header, setting.thisrun_timestamp)

        self.write_header_for_text_output \
            (NX, setting.txtout_items_on_line_maximum, \
             setting.txtout_file_breadth, 'breadth', 'm',  \
             setting.txtout_writedata_header, setting.thisrun_timestamp)

        return setting

#==============================================================================
    def read_datapairs(self, openfile, isnoisy):
        '''
        Reads data pairs from a text file
        '''
        import numpy as np
        
        usenextdata = False
        dontquit = True
        
        ii = 0
        while dontquit == True:
            rln = openfile.readline()
            if len(rln) == 0:
                if isnoisy == True:
                    print('end of file')
                break
            rln = rln.rstrip('\n')
            if usenextdata == False:   
                pln = rln.partition(" ")
                if pln[0] == 'HEADER:':
                    if isnoisy == True:
                        print(pln[2])
                elif pln[0] == 'FILETYPE:':
                    if isnoisy == True:
                        print(pln[2])
                elif pln[0] == 'DATE:':
                    if isnoisy == True:
                        print(pln[2])
                elif pln[0] == 'DATATYPE:':
                    if isnoisy == True:
                        print(pln[2])
                elif pln[0] == 'UNITS:':
                    if isnoisy == True:
                        print(pln[2])
                elif pln[0] == 'NX:':  
                    if isnoisy == True:
                        print(pln[2])                       
                    NX = int(pln[2])
                    data1 = np.zeros(NX,dtype=np.float)
                    data2 = np.zeros(NX,dtype=np.float)                    
                elif pln[0] == 'ITEMS_PER_LINE:':
                    if isnoisy == True:
                        print(pln[2])                                           
                elif pln[0] == 'END:':
                    if isnoisy == True:
                        print(pln[2])
                elif pln[0] == 'TIME:':
                    if isnoisy == True:
                        print(pln[2])                     
                elif pln[0] == 'DATA:':
                    usenextdata = True
                    if isnoisy == True:
                        print(pln[2])                                        
            else: # usednextata == True
                rln = rln.rstrip(' ')
                rln = rln.rstrip(']')
                rln = rln.lstrip('[')
                if isnoisy == True:
                    print('reading data')
                    print(rln)
                pln = rln.split()
                #for ii in range(len(pln)):
                    #pln[ii] = pln[ii].lstrip("  ")
                    #pln[ii] = pln[ii].rstrip("  ")
                    #print(ii, pln[ii])
                #pln = pln.lstrip(" ")
                #print(pln[0])
                #print(pln[1])
                #print(pln[2])
                din = [float(s) for s in pln]
                #print(NX, ii, data)
                data1[ii] = din[0]
                data2[ii] = din[1]
                ii=ii+1
                #for ii in range(0,len(pln)):
                #    print(ii, pln[ii])
                #    dtemp = float(pln[ii])
                #break
                if ii == NX:
                    dontquit = False    
                           
        return [data1, data2]
       
#==============================================================================        
    def read_max_number_pairs(self, setting, widthdepth_filename):
        '''
        reads through the file to determine the maximum number of data pairs
        '''
        openfile = open(widthdepth_filename,'r')
               
        dontquit = True
        npair = 0

        while dontquit == True:
            
            rln = openfile.readline()
            if len(rln) == 0:
                dontquit = False
                if setting.inoisy == True:
                    print('end of file ' ,widthdepth_filename )
                break
            #endif len(rln) == 0
            rln = rln.lstrip()
            rln = rln.rstrip('\n')
            pln = rln.partition(" ")

            if pln[0] == 'numberPairs':
                npair = max(npair,int(pln[2]))            
            #endif pln[0] == 'numberPairs'
        #endwhile dontquit == True
        
        #print(nPairs) 
        openfile.close() 
        
        return npair
    
#==============================================================================
    def read_number_of_cells(self, setting, widthdepth_filename):
        '''
        reads through the file to count the number of cell elements
        Assumes file is keyword item pairs on each line, and each
        new cell has "begin cell"
        '''
        openfile = open(widthdepth_filename,'r')
               
        dontquit = True
        cellcount = 0

        while dontquit == True:
            
            rln = openfile.readline()
            if len(rln) == 0:
                dontquit = False
                if setting.inoisy == True:
                    print('end of file')
                break
            #endif len(rln) == 0
            
            rln = rln.lstrip()
            rln = rln.rstrip('\n')
            pln = rln.partition(" ")
            
            if pln[0] == 'begin':
                if pln[2] == 'cell':
                    cellcount = cellcount + 1
                #endif pln[2] == 'cell'
            #endif pln[0] == 'begin'
        #endwhile dontquit == True
        
        openfile.close() 
        
        return cellcount
    
#==============================================================================
    def read_restart_files(self,setting, trun):
        '''
        opens, reads, and closes restart files for volume and flowrate
        '''
        import os
        import sys
        
        # ensure we have an absolute path for the restart direcotry
        setting.restart_directory = self.absolute_path                        \
            (setting.restart_directory, setting.working_directory)

        # get the complete restart path
        setting.restart_flowratefile = os.path.join                           \
                  (setting.restart_directory, setting.restart_flowratefile) 

        setting.restart_volumefile = os.path.join                             \
                  (setting.restart_directory, setting.restart_volumefile) 
           
        # open the files    
        file_restart_flowrate = open(setting.restart_flowratefile,'r')
        file_restart_volume   = open(setting.restart_volumefile,'r')
        
        # read the line that is equal to or after the restart time
        [temp_time,volume_restart] = self.get_oneline_data_matching_time      \
            (setting.restart_time, file_restart_volume, False)

        [restart_time,flowrate_restart] = self.get_oneline_data_matching_time \
            (setting.restart_time, file_restart_flowrate, False)
         
        # close the files    
        file_restart_flowrate.close()
        file_restart_volume.close()
     
        # check that times matched on the volume and flowrate
        if restart_time != temp_time:
            print('problem with non-matching file time stamps',               \
                  restart_time, temp_time)
            sys.exit()
        #endif
        
        trun.restart_time = restart_time
        trun.flowrate_restart = flowrate_restart
        trun.volume_restart = volume_restart
        
        return trun

#==============================================================================
    def read_widthdepth_pairs(self, geo, setting, NX, widthdepth_filename):
        '''
        Reads a set of widthdepth pairs from a text file for each cell'
        Text file is organized with in the form of "keyword value" with
        the width depth pairs as two items per line.  Typical layout and 
        allowable keyword value pairs are as follow (where xxxx is a number)
            begin cell
            ID  xxxx
            Length xxxx
            xDistance xxxx
            zBottom xxxx
            Breadth xxxx
            TrapezoidAngle xxxx (in degrees)
            cellType channel_WidthDepthPairs
            numberPairs xxxx
            WidthDepthPairs follow
            xxxx xxxx
            xxxx xxxx
            end WidthDepthPairs
            end cell
            begin cell
            ....
            end cell
        '''
        import sys
        import numpy as np
        
        #----------------------------------------------------------------------                       
        # scan the file for the number of cells and max number of data pairs
        nCells = self.read_number_of_cells(setting, widthdepth_filename)
        #nPairs = self.read_max_number_pairs(setting, widthdepth_filename)
        
        #----------------------------------------------------------------------                       
        # aliases for indexes
        width  = setting.geometry_widthdepth_values['widthAtLayerTop']
        depth  = setting.geometry_widthdepth_values['depthAtLayerTop']
    
        #----------------------------------------------------------------------                       
        #error checking
        if nCells != NX:
            if (setting.geometry_add_downstream_buffer) and \
                (NX == nCells + 1):
                print('Expecting custom code for downstream buffer')
            else:
                print(nCells, NX)            
                print('error, mismatch between file and input')
                sys.exit()
            #endif
        #endif nCells != NX    
                
        #----------------------------------------------------------------------                       
        # open the files    
        openfile = open(widthdepth_filename,'r')
         
        #----------------------------------------------------------------------                       
        # loop controls
        nextdata_is_pair = False
        dontquit = True
        expecting_newcell = True
        
        # cell index        
        icell = -1
        
        #----------------------------------------------------------------------                       
        while dontquit == True:
            # read one line of data and parse.            
            rln = openfile.readline()
            rln = rln.lstrip()
            # check for end of file
            if len(rln) == 0:
                dontquit = False
                if setting.inoisy == True:
                    print('end of file ',widthdepth_filename )
                break
            # HACK skip the comment lines - later make these saved metadata
            if rln[1] == '#':  
                continue
            rln = rln.rstrip('\n')
            pln = rln.partition(" ")
            
            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            # parse headers where data is not a pair of numbers
            if nextdata_is_pair == False:  
                # re-initialize for a new cell
                if pln[0] == 'begin':
                    if pln[2] == 'cell':
                        if expecting_newcell != True:
                            print(rln)
                            print('error, likely misalignment in file')
                            sys.exit()
                        else:
                            # reset for a new cell
                            expecting_newcell = False
                            nextdata_is_pair = False
                            icell = icell+1
                    else:    
                        print ('error, not yet designed for features other than cells')
                        sys.exit()   
                    #endif pln[2] == 'cell'    
                elif pln[0] == 'ID':
                    geo['ID'][icell] = int(pln[2])
                elif pln[0] == 'ManningsN':
                    geo['manningsn'][icell] = float(pln[2])
                elif pln[0] == 'Length':
                    geo['length'][icell] = float(pln[2])
                elif pln[0] == 'zBottom':
                    geo['zbottom'][icell] = float(pln[2])
                elif pln[0] == 'xDistance':
                    geo['xvalue'][icell] = float(pln[2])
                elif pln[0] == 'Breadth':
                    geo['breadth'][icell] = float(pln[2])
                elif pln[0] == 'TrapezoidAngle':
                    geo['gA'][icell] = np.deg2rad(float(pln[2]))
                elif pln[0] == 'cellType':
                    if pln[2] == 'channel_WidthDepthPairs':
                        geo['etype'][icell] = 'widthdepth_pair'
                    else:
                        pln
                        print('error: unknown value for cellType (etype)')
                    #endif pln[0] == 'channel_WidthDepthPairs'  
                elif pln[0] == 'numberPairs':
                    geo['npair'][icell] = int(pln[2])
                elif pln[0] == 'WidthDepthPairs':
                    if pln[2] == 'follow':
                        nextdata_is_pair = True
                        ipair = -1
                    else: 
                        print('error, not designed for 2nd argument other than follow')
                elif pln[0] == 'end':
                    if pln[2] == 'WidthDepthPairs':
                        # store the data
                        #geo['widthdepth'][icell,0:geo['npair'][icell]-1,0:2] \
                        #    = this_data 
                        this_data = None
                    elif pln[2] == 'cell':
                        # end the cell
                        expecting_newcell = True
                        #print('cell ',this_ID,' ended')                        
                    else:
                        print(rln)
                        print(pln[2])
                        print('error, unknown option')
                        sys.exit()
                    #endif pln[2] == 'WidthDepthPairs'    
                #endif pln[0] == 'begin'  
            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            else:
                if ipair == -1:
                    this_data                                                 \
                        = np.zeros([2,geo['npair'][icell]],dtype=np.float64)
                    ipair = 0

                #endif ipair == -1 
                if pln[0] == 'end':
                    #print(ipair, geo['widthdepth'].shape, this_data.shape)
                    geo['widthdepth'][icell,0:ipair,width]                    \
                        = this_data[0,0:ipair]
                    geo['widthdepth'][icell,0:ipair,depth]                    \
                        = this_data[1,0:ipair]
                    nextdata_is_pair = False
                    #print('data pairs ended')
                else:
                    #print('in else',ipair, geo['npair'][icell]-1)
                    if ipair > geo['npair'][icell]-1:
                        print('error, mismatch in expected number of pairs')
                        sys.exit()
                    #endif ipair > this_numberPairs-1
                    this_data[0,ipair] = float(pln[0])
                    this_data[1,ipair] = float(pln[2])                    
                    ipair = ipair+1
                    
                #endif pln[0] == 'end'
            #endif nextdata_is_pair == False 
            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        #endwhile dontquit == True
        #----------------------------------------------------------------------                       
        
        openfile.close() 
        
        return geo
    
#==============================================================================
    def set_working_directory(self, setting):
        ''' 
        If the source code is executed from ../src, this will back up
        the directory by one level for writing.reading the output, restart etc.                         
        '''
        import os
        
        if setting.working_directory == None: 
            setting.working_directory = os.getcwd()
            temp = os.path.split(setting.working_directory)
            if temp[1] == 'src':
                setting.working_directory = temp[0]
                
        return setting
    
#==============================================================================
    def time_loop_datasave(self, el, setting, trun):
        '''
        Writes output to binaru files
        '''
        
        # HACK need to check that this is working. Suspect trun caused problem
        
        import DataStorage2
        data = DataStorage2.DataStorage()
        
        if setting.binsave == True:
            if trun.present_time >=                                           \
                (trun.time_next_save_element - 0.499*setting.dt):
                           
                trun.element_data_saved = data.saveslice                      \
                    (el, trun.element_data_saved, setting.binsave_element,    \
                     trun.index_next_save_element)
                    
                trun.element_time_saved = data.savetimes                      \
                            (trun.element_time_saved, trun.present_step,      \
                             trun.present_time, trun.index_next_save_element)  
                            
                trun.time_next_save_element = trun.present_time               \
                                            + setting.binsave_timeinterval
                                            
                trun.index_next_save_element = trun.index_next_save_element+1
                
            #endif
        #endif
        
        return trun

#==============================================================================
    def time_loop_txtout_write(self, el, fa, setting, trun, NX):  
        '''
        Writes the text output to open files.
        '''
        import numpy as np
        import math
        
#       HACK - this could be restructured for a much simpler formattin
        
        if setting.txtout_writedata == True and \
            trun.present_time >= setting.txtout_writedata_starttime:  
                
            # check to see if it is time to output text data    
            if trun.present_time >=  \
                (trun.time_next_txtout_writedata - 0.499*setting.dt):
                
                # set the time counter for the next write
                trun.time_next_txtout_writedata = trun.present_time \
                    + setting.txtout_writedata_timeinterval                
                
                # write a TIME header line to each file
                setting.txtout_file_volume.write \
                    (''.join(['TIME: ',str(trun.present_time),'\n']))
                setting.txtout_file_volume.write('DATA: NEXTLINE\n')

                setting.txtout_file_flowrate.write \
                    (''.join(['TIME: ',str(trun.present_time),'\n']))
                setting.txtout_file_flowrate.write('DATA: NEXTLINE\n')

                setting.txtout_file_flowrate_face.write \
                    (''.join(['TIME: ',str(trun.present_time),'\n']))
                setting.txtout_file_flowrate_face.write('DATA: NEXTLINE\n')

                setting.txtout_file_hyddepth.write \
                    (''.join(['TIME: ',str(trun.present_time),'\n']))
                setting.txtout_file_hyddepth.write('DATA: NEXTLINE\n')

                setting.txtout_file_elevation.write \
                    (''.join(['TIME: ',str(trun.present_time),'\n']))
                setting.txtout_file_elevation.write('DATA: NEXTLINE\n')

                setting.txtout_file_elevation_faceM.write \
                    (''.join(['TIME: ',str(trun.present_time),'\n']))
                setting.txtout_file_elevation_faceM.write('DATA: NEXTLINE\n')

                setting.txtout_file_elevation_faceP.write \
                    (''.join(['TIME: ',str(trun.present_time),'\n']))
                setting.txtout_file_elevation_faceP.write('DATA: NEXTLINE\n')

                setting.txtout_file_area.write \
                    (''.join(['TIME: ',str(trun.present_time),'\n']))
                setting.txtout_file_area.write('DATA: NEXTLINE\n')

                setting.txtout_file_timestep.write \
                    (''.join(['TIME: ',str(trun.present_time),'\n']))
                setting.txtout_file_timestep.write('DATA: NEXTLINE\n')
                
                # number of lines ot be written. This is needed because
                # a large system might need to be written over multiple lines
                nline = math.ceil(NX / setting.txtout_items_on_line_maximum)
                
                idxL = 0
                idxH = setting.txtout_items_on_line_maximum
                #print('here ',nline, NX, setting.txtout_items_on_line_maximum)
                #print(setting.dt, np.shape(setting.dt))
                ar0 = np.zeros(1,dtype=np.float64 )
                ar0[0] = setting.dt
                #print(ar0, ar0.shape)
                #ts0 = np.array_str(setting.dt,precision = 16)
                ts0 = np.array_str(ar0,precision = 16)
                setting.txtout_file_timestep.write(ts0)
                setting.txtout_file_timestep.write('\n')
                for ii in range(0,nline):
                    #ts = temp[idxL:idxH]
                    #ts = np.array_str(ts,precision = 16)
                    #setting.txtout_file_volume.write(ts)
                    
                    ts1 = el['volume'][idxL:idxH]
                    ts2 = el['flowrate'][idxL:idxH]
                    ts3 = el['hyddepth'][idxL:idxH]
                    ts4 = el['eta'][idxL:idxH]
                    ts5 = el['area'][idxL:idxH]
                    ts6 = fa['flowrate'][idxL:idxH]
                    ts7 = fa['etaM'][idxL:idxH]
                    ts8 = fa['etaP'][idxL:idxH]
                    
                    ts1 = np.array_str(ts1,precision = 16)
                    ts2 = np.array_str(ts2,precision = 16)
                    ts3 = np.array_str(ts3,precision = 16)
                    ts4 = np.array_str(ts4,precision = 16)
                    ts5 = np.array_str(ts5,precision = 16)
                    ts6 = np.array_str(ts6,precision = 16)
                    ts7 = np.array_str(ts7,precision = 16)
                    ts8 = np.array_str(ts8,precision = 16)
                    
                    setting.txtout_file_volume.write(ts1)
                    setting.txtout_file_volume.write('\n')
                    setting.txtout_file_flowrate.write(ts2)
                    setting.txtout_file_flowrate.write('\n')
                    setting.txtout_file_hyddepth.write(ts3)
                    setting.txtout_file_hyddepth.write('\n')
                    setting.txtout_file_elevation.write(ts4)
                    setting.txtout_file_elevation.write('\n')
                    setting.txtout_file_area.write(ts5)
                    setting.txtout_file_area.write('\n')
                    
                    setting.txtout_file_flowrate_face.write(ts6)
                    setting.txtout_file_flowrate_face.write('\n')

                    setting.txtout_file_elevation_faceM.write(ts7)
                    setting.txtout_file_elevation_faceM.write('\n')
                    
                    setting.txtout_file_elevation_faceP.write(ts8)
                    setting.txtout_file_elevation_faceP.write('\n')

                    idxL = idxH
                    idxH = idxL + setting.txtout_items_on_line_maximum 
                #endfor 
            #endif
        #endif                   
      
        return trun

#==============================================================================
    def write_control_setting(self, setting):
        '''
        writes a file with all the simulation settings 
        '''                    
        thisfile   = open(setting.setting_filename,'w')
        
        attlist = dir(setting)
        for ii in range(0,len(attlist)):
            if attlist[ii][0:2] != '__':
                thisfile.write(attlist[ii]+' , ')  
                titem = getattr(setting,attlist[ii])
                stritem = str(titem)
                thisfile.write(stritem)
                thisfile.write('\n')
            #endif
        #endfor
        
        thisfile.close()   
        
        return
    

#==============================================================================

    def write_header_for_text_output \
        (self, NX, thisitems_per_line, thisfile, thisdata_type,              \
         thisdata_units, this_header, timestamp ): 
        ''''
        Writes a header for the data file that can be read in this class
        '''
        import datetime
        
        if timestamp == []:
            timestamp = datetime.datetime.now().strftime('%Y%m%d_%H%M')
        #endif    
        if this_header != []:
            thisfile.write('HEADER: '+this_header+'\n')
        #endif    
        thisfile.write('FILETYPE: nTIME-1DSPACE-1DATATYPE\n')
        thisfile.write('DATE: '+timestamp+'\n')
        thisfile.write('DATATYPE: '+thisdata_type+'\n')
        thisfile.write('UNITS: '+thisdata_units+'\n')
        thisfile.write('NX: '+str(NX)+'\n') 
        thisfile.write('ITEMS_PER_LINE: '+str(thisitems_per_line)+'\n')
        thisfile.write('END: HEADER\n')

        return
#==============================================================================
