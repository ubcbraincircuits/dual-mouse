'''
This script is able to record video of the mouse's behaviour.
It uses encoded h264 files to avoid dropping frames.

The script has automatic gain setters as well trigerring from a server pi.

Written by: Federico Bolanos, and Luis Bolanos, The University of British Columbia
Dr Timothy Murphy
'''

import picamera
import RPi.GPIO as GPIO
from time import sleep
from datetime import datetime
import io
import numpy as np

#GPIO Pins
client_pin = 23 #Input pin from server pi that triggers the start and end of the recording to synchronize videos.

#Video  File settings that will generate a full filename
directory = "/media/pi/Data/" #Directory to save video files
mouse_id = "M001" #ID of mouse
experiment_number = "1" #experiment number
conditions = "Laval_Demo" #Condition/ Brief Description of experiment

#Camera Settings
framerate = 30 #Since we are recording at such a low resolution, the pi is able to record at 90 fps which allows he rcording of quick whisker movements. The H264 format prevents frame drop as the video is encoded
resolution = (320,180) #Low resolution to record at high fps
#Unlike the brain imager, we do not need to lock and set specific awb_gains and exposure speeds as the camera, when is first initialized, sets its shutter speed to prevent over/under exposure



def pulse_pin(pin):
    GPIO.output(pin, True)
    sleep(0.001)
    GPIO.output(pin, False)

def get_filename():
    now = datetime.now()
    filename = mouse_id + "_"
    filename += "_experiment-" + experiment_number + "_" + conditions

    return filename

GPIO.setmode(GPIO.BCM)
GPIO.setup(client_pin, GPIO.IN, pull_up_down=GPIO.PUD_DOWN)


with picamera.PiCamera() as camera:
    camera.resolution = resolution
    camera.framerate = framerate
    camera.start_preview()
    sleep(2)
    g = camera.awb_gains
    camera.awb_mode = 'off'
    camera.awb_gains = g
    camera.awb_gains = (1,1)
    camera.shutter_speed = camera.exposure_speed
    camera.exposure_mode = 'off'
    #camera.shutter_speed = 6000000
    sleep(2)
    camera.stop_preview()

    print ("Waiting...")
    filename = get_filename()
    print(filename)
    camera.start_preview(fullscreen=False, window=(10, 10, 640, 480))
    GPIO.wait_for_edge (client_pin, GPIO.RISING)
    print ("Recording")
    camera.start_recording(directory+filename+".h264", format='h264',
                           quality=10, total_frames= 67000, timestamps_fn=directory+filename+"_timestamps.raw")
    while(GPIO.input(client_pin) != GPIO.LOW):
        camera.wait_recording(0.01)
    
    camera.stop_recording()
    camera.stop_preview()
    print("Saving...")                     
print ("Done!")

GPIO.cleanup()
