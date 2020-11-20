'''
This script is able to record video of the mouse's brain.
It pre-allocates the memory that will be used in order to avoid dropping frames.

The script has automatic gain setters as well as intensity checkers.

Written by: Federico Bolanos, and Luis Bolanos, The University of British Columbia
Dr Timothy Murphy
'''

import picamera
from time import sleep, time
import numpy as np
from picamera.array import PiRGBAnalysis
import random
import matplotlib
matplotlib.use('QT4Agg')
from matplotlib import pyplot as plt
import RPi.GPIO as GPIO
from matplotlib.animation import FuncAnimation
from datetime import datetime
import os


#Is it the server/master pi? This can be set to false if you are recording a second mouse at the same time requiring synchronized recordings
# The other mouse will have isServer set to true which will output a rising GPIO to the trigger_output_pin pin which will trigger the behaviour camera's
# client pin and this mouse's triggered_input_pin pin
isServer = True
#GPIO Pins
non_server_input_pin = 0 #if this pi is used as a client brain imager, this will be the input pin it waits for a rising edge to begin recording
trigger_behaviour_pin = 24 #this pin is used to output a high signal as soon as recording begins to trigger behaviour pi
trigger_server_output_pin = 23 # If you have multiple setups brain imagers, this pin will trigger high as soon as recording begins
trigger_other = 18
#Video File Settings, directory and naming (mouse_id, experiment_number, conditions and the date will be used to generate a complete filename
todays_date = "20201119"
directory = "/media/pi/Data/" + todays_date + "/" #Directory to save the recording
try:
    os.mkdir(directory)
except:
    pass
mouse_id = "PG3L-s_0000R-m" #ID of the mouse used in the experiment
experiment_number = "1" #experiment number
conditions = "test" #conditions

### VIDEO SETTINGS ###
height = 256 #Recording resolution height
width = 256 #width
h_height = 1024 #High quality preview resolution for focusing
h_width = 1024
length = 20 #Length of recording in seconds

### CAMERA SETTINGS ###
iso = 800 # It is preferred to keep the iso at a high level (100,200,400,800 are accepted values) If scene is too bright, reduce shutter_speed
shutter_speed = 33333 # Shutter speed is length of time in micro seconds. At 30 frames per second, a 360 degree exposure (camera exposes during the entire length of the frame) is equal to 33,333 microseconds
framerate = 30 #Framerate to capture. Making this a high number may cause the pi to drop frames
auto_white_balance_gains = (1,1) #color balance
vertical_flip = False   #Flips the image vertically
horizontal_flip = True #Flips the image horizontally

## These settings should not be adjusted
auto_white_balance_mode = 'off' # The AWB_Mode setting automatically adjusts the ratio of AWB_Gains to create a "correctly" coloured image/video in real-time
    # The bad thing about this is that we want raw data, and having the camera adjust a scene to make it look "right" when each channel is used as a
    # unique souce of information, will cause faulty data to appear.
exposure_mode = 'off' # similar to AWB_Mode, Exposure_modes tend to adjst the shutter_speed in real-time to keep the average intensity at a middle grey
    # So the camera will automaticaly adjust the shutter if it is looking at a very bright object to prevent over-exposure or a very dark object preventing
    # pure black regions. Turning this to "off" will remove the auto-adjustment
   




class OnlineAnalysis(PiRGBAnalysis):
    def __init__(self, camera, size):
        super(PiRGBAnalysis, self).__init__(camera, size=size)
        print("Called constructor for intensity viewer")
        self.intensity_array = np.empty((width, 3), dtype=np.uint8)

    # Array comes as a (64, 64, 3) array 
    def analyze(self, array):
        #self.intensity_array_col = np.max(array,axis=1)
        self.intensity_array_col = np.median(array,axis=1)
        #self.intensity_array_row = np.max(array,axis=0)
        self.intensity_array_row = np.median(array,axis=0)
        self.intensity_array = np.concatenate(
            (self.intensity_array_col, self.intensity_array_row),
            axis=1)
        #self.intensity_array = array[:, width//2, :]

class FrameCounter(PiRGBAnalysis):
    def __init__(self, camera, size):
        super(PiRGBAnalysis, self).__init__(camera, size=size)
        print("Called constructor for frame counter")
        self.frame_count = 0

    # Array comes as a (64, 64, 3) array 
    def analyze(self, array):
        self.frame_count+=1

    def __del__(self):
        print("Total frames: ", self.frame_count)
        
class MainClass():
    def __init__(self):
        
        
        #LED SETUP
        GPIO.setmode(GPIO.BCM)
        
        GPIO.setup(non_server_input_pin, GPIO.IN, pull_up_down=GPIO.PUD_DOWN)
        GPIO.setup(trigger_behaviour_pin, GPIO.OUT, initial=GPIO.LOW)
        GPIO.setup(trigger_server_output_pin, GPIO.OUT, initial=GPIO.LOW)
        GPIO.setup(trigger_other, GPIO.OUT, initial=GPIO.LOW)

        
    def __del__(self):
        self.camera.close()
        self.output.close()
        GPIO.cleanup()
        print("camera and output objects destroyed.")


    def setup_camera(self):
        self.camera = picamera.PiCamera()
        self.camera.resolution = (width, height)
        self.camera.framerate = framerate
        self.camera.vflip = vertical_flip
        self.camera.hflip = horizontal_flip
        self.camera.start_preview()
        sleep(0.5)
        self.camera.shutter_speed = shutter_speed
        self.camera.iso = iso
        self.camera.exposure_mode = exposure_mode
        self.camera.awb_mode = auto_white_balance_mode
        self.camera.awb_gains = auto_white_balance_gains
        self.drc_strength = 'off'
        self.camera.video_denoise = False
        sleep(0.1)
        self.camera.stop_preview()
        print()
        print("FINAL Analog gains:  ", float(self.camera.analog_gain))
        print("FINAL Digital gains: ", float(self.camera.digital_gain))
        print("FINAL Shutter speed: ", self.camera.shutter_speed)
        print("FINAL ISO: ",self.camera.iso)
        print("FINAL Auto White Balance Gains: ",self.camera.awb_gains)
        print("Exposure Mode: ", self.camera.exposure_mode)
        print("Auto White Balance Mode: ",self.camera.awb_mode)


    def display_high_res_preview(self):
        self.camera.resolution = (h_width, h_height)
        self.camera.start_preview()
        print("Press CTRL-C to stop the high resolution preview.")
        try:
            while True:
                sleep(0.1)
        except KeyboardInterrupt:
            self.camera.stop_preview()
            self.camera.resolution = (width, height)

    def data_gen(self):
       yield self.output.intensity_array


    def update_plot(self, data):
        self.line_r1.set_ydata(data[:, 0])
        self.line_g1.set_ydata(data[:, 1])
        self.line_b1.set_ydata(data[:, 2])
        self.line_r2.set_ydata(data[:, 3])
        self.line_g2.set_ydata(data[:, 4])
        self.line_b2.set_ydata(data[:, 5])
        #print(self.camera.analog_gain)

        return self.line_r1, self.line_g1, self.line_b1, self.line_r2, self.line_g2, self.line_b2
    


    def setup_intensity(self):
        self.camera.start_preview(fullscreen=False, window=((10, 10, 256, 256)))
        self.output = OnlineAnalysis(self.camera, size=(width, height))
        self.camera.start_recording(self.output, format='rgb')

        for i in range(1):
            # Plot in real time.
            self.fig, (self.ax1, self.ax2) = plt.subplots(1,2)
            self.line_r1, = self.ax1.plot(self.output.intensity_array[:, 0], 'r')
            self.line_g1, = self.ax1.plot(self.output.intensity_array[:, 1], 'g')
            self.line_b1, = self.ax1.plot(self.output.intensity_array[:, 2], 'b')

            self.line_r2, = self.ax2.plot(self.output.intensity_array[:, 3], 'r')
            self.line_g2, = self.ax2.plot(self.output.intensity_array[:, 4], 'g')
            self.line_b2, = self.ax2.plot(self.output.intensity_array[:, 5], 'b')

            self.ax1.set_ylim(0, 255)
            self.ax1.set_xlim(0, width)
            self.ax1.title.set_text('Max across columns')

            self.ax2.set_ylim(0, 255)
            self.ax2.set_xlim(0, width)
            self.ax2.title.set_text('Max across rows')


            ani = FuncAnimation(self.fig, self.update_plot, self.data_gen, interval=30)
            plt.show()

        self.camera.stop_recording()
        self.camera.stop_preview()


    def record_video(self, length):
        full_filename = self.get_full_filename()

        print(full_filename)

        # allocate memory for video data
        total_frames = int(length*self.camera.framerate)
        print("Will record a total of: ", total_frames, " frames.")

        # trigger other imaging pi, behaviour and audio.
        #self.output = FrameCounter(self.camera, size=(width, height))
        #frames = []
        input("Press ENTER to begin recording and trigger pin HIGH...")
        print("Recording!")
        GPIO.output(trigger_other,GPIO.HIGH)
        
        self.camera.start_preview(fullscreen=False, window=((10, 10, 512, 512)))
        GPIO.output(trigger_server_output_pin, GPIO.HIGH)
        GPIO.output(trigger_behaviour_pin, GPIO.HIGH)
        
        #self.camera.start_recording(self.output, format='rgb', splitter_port=2)
        self.camera.start_recording(full_filename+".raw", format='rgb', total_frames=total_frames, timestamps_fn=full_filename+"_timestamps.raw")
        self.camera.wait_recording(length)
            
        GPIO.output(trigger_behaviour_pin, GPIO.LOW)
        GPIO.output(trigger_server_output_pin, GPIO.LOW)
        GPIO.output(trigger_other,GPIO.LOW)
        self.camera.stop_recording()
        #self.camera.stop_recording(splitter_port=2)
        self.camera.stop_preview()

        print("Done!")


    def save_array(self, array, filename):
        with open(filename, "w") as file:
            for entry in array:
                file.write(str(entry)+"\n")

                
    def record_triggered_video(self, length):
        full_filename = self.get_full_filename()
        
        print(full_filename)

        total_frames = int(length*self.camera.framerate)
        print("Will record a total of: ", total_frames, " frames.")

        # wait for trigger from server imaging pi
        self.camera.start_preview(fullscreen=False, window=((10, 10, 512, 512)))
        print("Waiting for RISING edge...")
        GPIO.wait_for_edge(non_server_input_pin, GPIO.RISING)
        print("Recording!")
        self.camera.start_recording(full_filename+".raw", format='rgb', total_frames=total_frames, timestamps_fn=full_filename+"_timestamps.raw")
        #start recording and output the trigger pin 
        GPIO.output(trigger_behaviour_pin, GPIO.HIGH)
        self.camera.wait_recording(length)
        self.camera.stop_recording()
        GPIO.output(trigger_behaviour_pin, GPIO.LOW)
        self.camera.stop_preview()


        print("Done!")

    def capture_bayer(self, filename):
        self.camera.capture(filename, format='jpeg', bayer=True)
        
        
    def get_full_filename(self):
        now = datetime.now()
        filename = "M" + mouse_id + "_" + now.strftime("%B-%d_%H%M")
        filename += "_experiment-" + experiment_number + "_" + conditions

        full_filename = directory + filename

        return full_filename

if __name__ == '__main__':
    mc = MainClass()
    mc.setup_camera()
    mc.display_high_res_preview()
    mc.setup_intensity()
    if isServer:
        mc.record_video(length)
    else:
        mc.record_triggered_video(length)
    print("All done.")
    mc.__del__()
    
