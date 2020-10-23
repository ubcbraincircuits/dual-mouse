"""
librain, the Brain Library.

Created to support data processing for the study, "Mesoscale cortical calcium imaging reveals widespread synchronized infraslow 
activity during social touch in mice" by the Murphy Lab at the University of British Columbia.

This module can generate file names, return a complete path to a file, and save results such that the naming convention
used by the source data is maintained.
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio
import os
from os.path import join, getsize, isfile, isdir
from pathlib import Path
from dateutil.parser import parse 
import sys, traceback
class Data:


    def __init__(self, directory):
        self.directory = directory
        self.dates = os.listdir(self.directory)


    def experiment(self, date, exp_num, listfiles=False):
        """
        Returns the directory of the experiment if it exists and if not,
        an exception is raised.

        :param date: the date the experiment was conducted in any format
        :type: str
        :param exp_num: the experiment number
        :type: int 
        :param listfiles: if False, only the first value below is returned; otherwise, 
        both are returned
        :type: bool

        :return: full path to the directory of the experiment
        :type: str
        :return: all files in the experiment directory
        :type: list
        """
        if type(date) == str:
            d_format = parse(date)
            d = f"/{d_format.year}{d_format:%m}{d_format:%d}"
            date_path = Path(str(self.directory) + d) 
            if isdir(date_path) is False:
                raise Exception(f'Folder {d} does not exist')
            else:
                exp = f"/Experiment_{exp_num}"
                exp_folder = Path(str(date_path) + exp)
                if isdir(exp_folder) is True:
                    if listfiles is True:
                        files = [file for file in os.listdir(exp_folder)]
                        return str(exp_folder), files
                    else:
                        return str(exp_folder)
                else:
                    raise Exception(f'Folder Experiment {exp_num} does not exist in folder {d}')
                

    def file(self, exp_folder, fname, subfolder=""):
        """
        Returns complete path to fname if it exists in exp_folder or within its subfolders, Behaviour or Derivatives. 
        If fname is not in exp_folder, an exception is raised.

        :param exp_folder: complete path to experiment folder, i.e. Experiment 1
        :type: str
        :param fname: one of 
        'timestamps', 
        'subset interpolated',
        'interpolated', 
        'h264', 
        'combined', 
        'processed',
        'RM mask'
        'LM mask', 
        'left blue', 
        'left green', 
        'right blue',
        'right green', 
        'left blue 0.01-3.0Hz', 
        'left green 0.01-3.0Hz', 
        'right blue 0.01-3.0Hz',
        'right green 0.01-3.0Hz',
        'left blue 0.01-12.0Hz', 
        'left green 0.01-12.0Hz', 
        'right blue 0.01-12.0Hz',
        'right green 0.01-12.0Hz',
        'left 0.01-12.0Hz',
        'right 0.01-12.0Hz',
        'left 0.01-3.0Hz',
        'right 0.01-3.0Hz',
        'left',
        'right',
        'left timestamps',
        'right timestamps'
        'freq split ws=3800',
        'freq split ws=5000',
        'freq split ws=1000',
        'freq split ws=2000',
        'freq split ws=2800'
        :type: str 
        :param subfolder: can be specified as 'Behaviour' or 'Derivatives'
        :type: str

        :return: full path to fname 
        :type: str
        """
        fnames = [
            'timestamps',
            'subset interpolated',
            'interpolated',
            'h264', 
            'combined', 
            'processed',
            'RM mask',
            'LM mask', 
            'left blue', 
            'left green', 
            'right blue',
            'right green', 
            'left blue 0.01-3.0Hz', 
            'left green 0.01-3.0Hz', 
            'right blue 0.01-3.0Hz',
            'right green 0.01-3.0Hz',
            'left blue 0.01-12.0Hz', 
            'left green 0.01-12.0Hz', 
            'right blue 0.01-12.0Hz',
            'right green 0.01-12.0Hz',
            'left 0.01-12.0Hz',
            'right 0.01-12.0Hz',
            'left 0.01-3.0Hz',
            'right 0.01-3.0Hz',
            'left',
            'right',
            'left timestamps',
            'right timestamps',
            'freq split ws=3800',
            'freq split ws=5000',
            'freq split ws=1000',
            'freq split ws=1500',
            'freq split ws=1750',
            'freq split ws=2500',
            'freq split ws=864',
            'freq split ws=432',
            'freq split ws=2000',
            'freq split ws=2800',
            'trunc',
            'left gsr',
            'right gsr',
            'left global',
            'right global',
            'matlab'
            ]

        if fname not in fnames:
            raise ValueError(f'{fname} is not a valid filename. Check help(<directory>.file) for a list of filenames')

        if subfolder == 'Behaviour':   
            direc = os.path.join(str(exp_folder), subfolder)        
            for root, dirs, files in os.walk(direc): 
                for f in files:
                    if 'timestamps' in f and fname == 'timestamps':
                        return str(Path(os.path.join(root, f)))
                    elif 'subset_interpolated' in f and fname == 'subset interpolated':
                        return str(Path(os.path.join(root, f)))
                    elif 'interpolated' in f and fname == 'interpolated':
                        return str(Path(os.path.join(root, f)))
                    elif 'h264' in f and fname == 'h264':
                        return str(Path(os.path.join(root, f))) 

            raise FileNotFoundError(f'File {fname} does not exist in subfolder {subfolder}') 
        
        elif subfolder == "" or subfolder == "Derivatives":
            direc = os.path.join(str(exp_folder), subfolder) 
            for root, dirs, files in os.walk(direc):
                for f in files:
                    if 'combined' in f and 'raw' in f and 'upscaled' not in f:
                        if not 'gsr' in f:
                            
                            if '0.01-3' in f:
                                if fname == 'combined':
                                    return Path(os.path.join(root, f))
                    elif 'left_mouse_gsr' in f and fname == "left gsr":
                        return str(Path(os.path.join(root, f)))
                    elif 'right_mouse_gsr' in f and fname == "right gsr":
                        return str(Path(os.path.join(root, f)))
                    elif 'l_global_signal' in f and fname == "left global":
                        return str(Path(os.path.join(root, f))) 
                    elif 'r_global_signal' in f and fname == "right global":
                        return str(Path(os.path.join(root, f))) 
#                     elif 'snips' in f and fname == "matlab":
#                         return str(Path(os.path.join(root, f)))                    
                    elif 'processed' in f and fname == 'processed':
                        return str(Path(os.path.join(root, f)))
                    elif 'RM_mask' in f and fname == 'RM mask':
                        return str(Path(os.path.join(root, f)))
                    elif 'LM_mask' in f and fname == 'LM mask':
                        return str(Path(os.path.join(root, f)))

                    elif 'frequency_split_correlation_filtered_ws' in f:
                        if '3800' in f and fname == 'freq split ws=3800':
                            return str(Path(os.path.join(root, f)))
                        if '5000' in f and fname == 'freq split ws=5000':
                            return str(Path(os.path.join(root, f)))
                        if '1000' in f and fname == 'freq split ws=1000':
                            return str(Path(os.path.join(root, f)))
                        if '1500' in f and fname == 'freq split ws=1500':
                            return str(Path(os.path.join(root, f)))
                        if '1750' in f and fname == 'freq split ws=1750':
                            return str(Path(os.path.join(root, f)))
                        if '2500' in f and fname == 'freq split ws=2500':
                            return str(Path(os.path.join(root, f)))
                        if '864' in f and fname == 'freq split ws=864':
                            return str(Path(os.path.join(root, f)))
                        if '432' in f and fname == 'freq split ws=432':
                            return str(Path(os.path.join(root, f)))
                        if '2000' in f and fname == 'freq split ws=2000':
                            return str(Path(os.path.join(root, f)))
                        if '2800' in f and fname == 'freq split ws=2800':
                            return str(Path(os.path.join(root, f)))

                    elif 'BLUE' in f:
                        if 'LEFT' in f:
                            if 'RAW' in f and fname == 'left blue':
                                return str(Path(os.path.join(root, f)))
                            if '0.01-3.0' in f and fname == 'left blue 0.01-3.0Hz':
                                return str(Path(os.path.join(root, f)))
                            if '0.01-12.0' in f and fname == 'left blue 0.01-12.0Hz':
                                return str(Path(os.path.join(root, f)))
                        if 'RIGHT' in f:
                            if 'RAW' in f and fname == 'right blue':
                                return str(Path(os.path.join(root, f)))
                            if '0.01-3.0' in f and fname == 'right blue 0.01-3.0Hz':
                                return str(Path(os.path.join(root, f)))
                            if '0.01-12.0' in f and fname == 'right blue 0.01-12.0Hz':
                                return str(Path(os.path.join(root, f)))   
                    elif 'GREEN' in f:
                        if 'LEFT' in f:
                            if 'RAW' in f and fname == 'left green':
                                return str(Path(os.path.join(root, f)))
                            if '0.01-3.0' in f and 'TRUNCATED.npy' in f and fname == 'trunc':
                                return str(Path(os.path.join(root, f)))
                            if '0.01-3.0' in f and fname == 'left green 0.01-3.0Hz':
                                return str(Path(os.path.join(root, f)))
                            if '0.01-12.0' in f and fname == 'left green 0.01-12.0Hz' and "mp4" not in f:
                                return str(Path(os.path.join(root, f)))   
                        if 'RIGHT' in f:
                            if 'RAW' in f and fname == 'right green':
                                return str(Path(os.path.join(root, f)))
                            if '0.01-3.0' in f and fname == 'right green 0.01-3.0Hz':
                                return str(Path(os.path.join(root, f)))
                            if '0.01-12.0' in f and fname == 'right green 0.01-12.0Hz' and 'mp4' not in f:
                                return str(Path(os.path.join(root, f)))   
                    
                    elif 'LEFT_corrected' in f:
                        if '0.01-3.0' in f and fname == 'left 0.01-3.0Hz':
                            return str(Path(os.path.join(root, f)))
                        if '0.01-12.0' in f and fname == 'left 0.01-12.0Hz' and "mp4" not in f:
                            return str(Path(os.path.join(root, f)))
                    elif 'RIGHT_corrected' in f:
                        if '0.01-3.0' in f and fname == 'right 0.01-3.0Hz':
                            return str(Path(os.path.join(root, f)))
                        if '0.01-12.0' in f and fname == 'right 0.01-12.0Hz':
                            return str(Path(os.path.join(root, f)))

                    elif not 'LEFT' in f and not 'RIGHT' in f:
                        if f[5] == 'L':
                            if 'timestamps' in f and fname == 'left timestamps':
                                return str(Path(os.path.join(root, f)))
                            if not 'bandpass' in f and fname == 'left':
                                return str(Path(os.path.join(root, f)))
                        if f[5] == 'R':
                            if 'timestamps' in f and fname == 'right timestamps':
                                return str(Path(os.path.join(root, f)))
                            if not 'bandpass' in f and fname == 'right':
                                return str(Path(os.path.join(root, f)))
                            
                    
            raise FileNotFoundError(f'File {fname} does not exist in {exp_folder}',)

        else:
            raise FileNotFoundError(f'Subfolder {subfolder} does not exist in folder {exp_folder}')


class Output:


    def __init__(self, directory):
        self.directory = directory


    def saveas(self, f_out, ftype, suffix=None, prefix=None, dtype=None, f_in=None, save=False, path=None, dirname="Derivatives", fig=False):
        """
        Returns file name of a result with the same naming convention as the corresponding raw data. 
        The result can also be saved.

        :param f_out: the result to be saved
        :type: any
        :param suffix: word or phrase to append to the end of the file name (i.e. <file>_PROCESSED) OR unique file name
        :type: str
        :param ftype: desired file type
        :type: str
        :param prefix: word or phrase to append to the beginning of the file name (i.e. LEFT_GREEN_<file>)
        :type: str
        :param dtype: desired numpy data type, i.e. 'float32'; ftype must be set to 'raw'
        :type: str
        :param f_in: complete path to raw file from which the result was derived; leave out to customize file name
        :type: str 
        :param save: saves the result if True
        :type: bool
        :param path: complete path to the directory in which the file will be saved; if unspecified, the file
        will be saved in "Derivatives" folder created automatically within self.directory
        :type: str
        :param dirname: name of the directory in which the file will be saved
        :type: str
        :param fig: must be set to True if f_out is a figure; the file is saved in npy format otherwise
        :type: bool

        :return: file name of result
        :type: str
        """
        if path == None:
            direc = self.directory + f'/{dirname}/'
            if isdir(direc) is False:
                direc = os.mkdir(direc)
        else:
            if isdir(path) is False:
                raise FileNotFoundError(f'{path} does not exist')
            else:
                direc = path + f'/{dirname}/'
                if isdir(direc) is False:
                    direc = os.mkdir(direc)

        if f_in == None:
            fname = suffix + f'.{ftype}'
        else:
            if isfile(f_in) is False:
                raise FileNotFoundError(f'{f_in} does not exist')
            # isolate file name
            if '/' in f_in: 
                fname = f_in.split('/')[-1].split('.')
            elif '\\' in f_in:
                fname = f_in.split('\\')[-1].split('.')
            rm_type = ""
            for i in range(len(fname)-1): 
                # Removes file type. In case there are periods in the file name, 
                # concatenate up to but not including the file type. 
                if i == len(fname)-2:
                    rm_type += fname[i]
                else:
                    rm_type += fname[i] + '.'
            
            if suffix == None:
                if prefix == None:
                    raise ValueError('Specify a value for prefix. For unique file names, use suffix')
                fname = f'{prefix}_' + rm_type + f'.{ftype}'
            else:
                fname = rm_type + f'_{suffix}' + f'.{ftype}'
            
        path = Path(str(direc)+fname)

        if save is True:
            if fig is True:
                f_out.savefig(path)
                print(f'Saved as {fname}')
                return fname
            else:
                if ftype == 'raw':
                    if dtype == None:
                        f_out.tofile(direc+fname)
                        print(f'Saved as {fname}')
                        return fname
                    else:
                        f_out.astype(np.dtype(dtype)).tofile(direc+fname)
                        print(f'Saved as {fname}')
                        return fname
                else:
                    np.save(path, f_out)
                    print(f'Saved as {fname}')
                    return fname
        else:
            return fname 




