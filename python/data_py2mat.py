import video_processing as vp
import librain as lb
import fnames
import numpy as np
from roipoly import RoiPoly
import matplotlib

matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from scipy.ndimage.filters import gaussian_filter
import scipy.io as sio
import scipy.signal as sig
import scipy.stats as stats
import easygui
import warnings
import h5py

warnings.filterwarnings('ignore')

path = "B:/Dual/"
direc = lb.Data(path)

fs = 28.815
fc = 0.025
together_duration = np.round(120 * fs)
translation_duration = np.round(27.5 * fs)
first_translation = np.round(119.5 * fs)
start_first_translation = first_translation
end_first_translation = first_translation + translation_duration

start_interaction = first_translation + translation_duration
end_interaction = first_translation + translation_duration + together_duration

start_second_translation = first_translation + translation_duration + together_duration
end_second_translation = first_translation + translation_duration + together_duration + translation_duration


def draw_roi(frames, region):
    plt.imshow(frames[5000], cmap='gray', vmin=0, vmax=255)
    plt.title("Draw " + region + " ROI")
    roi = RoiPoly(color='r')
    mask = roi.get_mask(frames[0])

    gradient_signal = np.abs(np.gradient(np.mean(frames[:, mask], axis=1)))
    threshold = np.mean(gradient_signal) + np.std(gradient_signal)

    return gradient_signal, threshold


# SOCIAL EXPERIMENTS
# experiments = [['20190729', '1', '2', '3', '4'],
#                ['20190808', '1', '2', '3', '4', '5', '6', '7', '8'],
#                ['20190815', '5', '6'],
#                # ['20190821', '1'],
#                ['20190822', '1', '2', '3', '4'],
#                ['20190829', '2', '4', '6', '8']]


# experiments = [['20190821', '1']]
# experiments = [['20200721', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12']]
# experiments = [['20200721', '7', '8', '9', '10', '11', '12']]


# MESH EXPERIMENTS
experiments = [#['20190823', '1', '2', '3', '4', '5', '6', '7', '8'],
               ['20190826', '5', '6'],#'1', '2', '3', '4', '5', '6'],
               ['20190829', '3', '5']]

# OPAQUE EXPERIMENTS
# experiments = [['20190827', '1', '2', '3', '4'],
#               ['20190828', '1', '2', '3', '4', '5', '6', '7', '8', '9'],
#               ['20190829', '1', '7']]

# folder_path = "B:/Social_Outputs/Matlab/social/behaviour_added_corrected/test/"
# folder_path = "B:/Social_Outputs/Matlab/social/behavior_added/"
folder_path = "B:/Social_Outputs/Matlab/mesh/"
# folder_path = "B:/Social_Outputs/Matlab/opaque/"
for i in range(len(experiments)):
    for j in range(len(experiments[i]) - 1):
        EXP = direc.experiment(experiments[i][0], int(experiments[i][1 + j]))

        subset_behaviour_file = direc.file(exp_folder=EXP, fname="interpolated", subfolder="Behaviour")
        # SOCIAL EXPERIMENTS
        # if experiments[i][0] == '20190822' or experiments[i][0] == '20190829' or experiments[i][0] == '20200721':
        #     HEIGHT = 240
        #     WIDTH = 320
        # else:
        #     HEIGHT = 480
        #     WIDTH = 640

        # MESH  & OPAQUE EXPERIMENTS
        HEIGHT = 240
        WIDTH = 320

        # load behavior frames
        # behaviour_frames = vp.extract_RAW_frames(
        #     subset_behaviour_file,
        #     width=WIDTH, height=HEIGHT,
        #     dtype='uint8', num_channels=1)
        #
        # answer = True
        # left_whisk = None
        # left_forelimb = None
        # left_groom = None
        # right_whisk = None
        # right_forelimb = None
        # right_groom = None
        # while answer:
        #     # draw ROIs
        #     left_whisk, lwt = draw_roi(behaviour_frames, "Whisking")
        #     left_forelimb, lft = draw_roi(behaviour_frames, "Forelimb")
        #
        #     right_whisk, rwt = draw_roi(behaviour_frames, "Whisking")
        #     right_forelimb, rft = draw_roi(behaviour_frames, "Forelimb")
        #
        #     # plot gradient signal
        #     plt.figure()
        #     plt.subplot(211)
        #     plt.plot(gaussian_filter(left_whisk, 20), 'b')
        #     plt.axhline(lwt, color='b')
        #     plt.plot(gaussian_filter(left_forelimb, 30), 'r')
        #     plt.axhline(lft, color='r')
        #
        #     plt.subplot(212)
        #     plt.plot(gaussian_filter(right_whisk, 20), 'b')
        #     plt.axhline(rwt, color='b')
        #     plt.plot(gaussian_filter(right_forelimb, 30), 'r')
        #     plt.axhline(rft, color='r')
        #     plt.show()
        #
        #     # check to see if we try again
        #     answer = easygui.ynbox('Retry?', 'Title', ('Yes', 'No'))

        # load brain data
        l_mouse_processed_fileg = direc.file(exp_folder=EXP, fname="left green 0.01-12.0Hz", subfolder="Derivatives")
        r_mouse_processed_fileg = direc.file(exp_folder=EXP, fname="right green 0.01-12.0Hz", subfolder="Derivatives")

        l_mouse_processed_fileb = direc.file(exp_folder=EXP, fname="left blue 0.01-12.0Hz", subfolder="Derivatives")
        r_mouse_processed_fileb = direc.file(exp_folder=EXP, fname="right blue 0.01-12.0Hz", subfolder="Derivatives")

        # l_mouse_processed_file = direc.file(exp_folder=EXP, fname="left 0.01-12.0Hz", subfolder="Derivatives")
        # r_mouse_processed_file = direc.file(exp_folder=EXP, fname="right 0.01-12.0Hz", subfolder="Derivatives")

        l_mouse_framesg = vp.extract_RAW_frames(l_mouse_processed_fileg, 256, 256, dtype='float32', num_channels=1)
        r_mouse_framesg = vp.extract_RAW_frames(r_mouse_processed_fileg, 256, 256, dtype='float32', num_channels=1)

        l_mouse_framesb = vp.extract_RAW_frames(l_mouse_processed_fileb, 256, 256, dtype='float32', num_channels=1)
        r_mouse_framesb = vp.extract_RAW_frames(r_mouse_processed_fileb, 256, 256, dtype='float32', num_channels=1)

        name = l_mouse_processed_fileg.split("\\")
        name = name[-1:]

        # load frames for registration
        l_green_frame = np.load(direc.file(exp_folder=EXP, fname="left green"))
        r_green_frame = np.load(direc.file(exp_folder=EXP, fname="right green"))

        my_dict = {"left_dFF_green": l_mouse_framesg,
                   "right_dFF_green": r_mouse_framesg,
                   "left_frame": l_green_frame,
                   "right_frame": r_green_frame,
                   "left_dFF_blue": l_mouse_framesb,
                   "right_dFF_blue": r_mouse_framesb,
                   "t1": start_first_translation,
                   "i1": start_interaction,
                   "i2": end_interaction,
                   "t2": end_second_translation
                   # "left_whisk": left_whisk,
                   # "right_whisk": right_whisk,
                   # "left_forelimb": left_forelimb,
                   # "right_forelimb": right_forelimb
                   }
        sio.savemat(folder_path + name[0][:-4] + "20200901" + ".mat", my_dict)

        print("done")
