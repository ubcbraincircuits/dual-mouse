
import numpy
from joblib import cpu_count, delayed, Parallel
from scipy import signal
from sklearn.utils import gen_even_slices
import cv2
import pandas
import matplotlib.pyplot as plt
from roipoly import RoiPoly

def extract_RAW_frames(
    filename,
    width,
    height,
    channel="all",
    dtype="uint8",
    num_channels=3,
):
    """
   Extract channels from .RAW file containing image data

   :param filename: name of .RAW file containing image data
   :type: str
   :param width: width of data
   :type: int
   :param height: height of data
   :type: int
   :param channel: name of channel to extract
   :type: optional str, one of: ('all', 'red', 'blue', 'green'). default is 'all'
   :param dtype: numpy datatype. default is 'uint8'
   :type: optional str
   :param num_channels, one of (1,3). defualt is 3
   :type: optional int

   :return: channel(s) extracted from .RAW file
   :type: numpy.ndarray
   """
    if num_channels not in (1, 3):
        raise AttributeError(
            "Keyword 'num_channels' must be one of (1,3)"
        )
    try:
        datatype = getattr(numpy, dtype)
    except AttributeError:
        raise AttributeError(
            f"dtype numpy.{dtype} does not exist"
        )
    if num_channels is 3:
        if channel not in ("all", "red", "blue", "green"):
            raise AttributeError(
                "Keyword 'channel' must be one of: ('all', 'red', 'blue', 'green')"
            )
        with open(filename, "rb") as file:
            raw_frames = numpy.fromfile(
                file, dtype=datatype
            )
        time_dim_float = raw_frames.shape[0] / (
            width * height * 3
        )
        print(time_dim_float)
        time_dim = int(time_dim_float)
        if time_dim != time_dim_float:
            raise Exception(
                "Invalid input file or arguments"
            )
        raw_frames = numpy.reshape(
            raw_frames, (time_dim, height, width, 3)
        )

        channel = {"red": 0, "green": 1, "blue": 2}.get(
            channel
        )
        if channel is None:
            # return all frames
            return raw_frames
        else:
            # return a particular channel
            return raw_frames[..., channel]
    else:
        # num_channels is 1
        with open(filename, "rb") as file:
            raw_frames = numpy.fromfile(
                file, dtype=datatype
            )
            time_dim_float = raw_frames.shape[0] / (
                width * height
            )
            time_dim = int(time_dim_float)
            if time_dim != time_dim_float:
                raise Exception(
                    "Invalid input file or arguments"
                )
            raw_frames = numpy.reshape(
                raw_frames, (time_dim, height, width)
            )
            return raw_frames


def clean_raw_timestamps(filename):
    """
    Perform cleaning routine on .RAW timestamps

    :param filename: name of .RAW file containing timestamps
    :type:str

    :return: cleaned timestamp array
    :type: numpy.ndarray
    """
    with open(filename, "rb") as file:
        raw_timestamps = numpy.fromfile(
            file, dtype=numpy.float32
        )
    raw_timestamps[0] = 1
    return raw_timestamps[numpy.where(raw_timestamps > 0)]


def get_locations_of_dropped_frames(timestamps, threshold):
    """
    Return delta(timestamps) and indices of timestamps of dropped
    frames

    :param timestamps: 1-D numpy array containing timestamps
    :type: numpy.ndarray
    :param thereshold: threshold in microseconds. suggested 50,000 for 30fps,
                       12,500 for 90fps
    :type: float

    :return differences between timestamps; dt
    :type: numpy.ndarray
    :return indices of timestamps where frames were dropped
    :type: numpy.ndarray
    """
    differences = numpy.diff(timestamps, 1)
    print(
        "Mean filtered frame difference: ",
        numpy.mean(
            differences[
                numpy.where(differences <= threshold)
            ]
        ),
    )

    return (
        differences,
        numpy.where(differences > threshold)[0],
    )


def generate_frames(
    frames, differences, locations, true_frame_rate
):
    """
    Generator function for producing args to initialise list of DroppedFrames

    :param frames: frames of channel extracted from RAW file
    :type: numpy.ndarray
    :param differences: differences between timestamps
    :type: numpy.ndarray
    :param locations: indices of timestamps where frames were dropped
    :type: numpy.ndarray
    :param true_frame_rate: true frame rate of footage
    :type: float

    :return tuple:
        :first_frame: last frame before missing frame(s)
        :type: numpy.ndarray
        :last_frame: first frame after missing frame(s)
        :type: numpy.ndarray
        :number_of_dropped_frames: number of dropped frames
        :type: int
        :location: indices where frames are missing
        :type: numpy.ndarray
    """
    for location in locations:
        number_of_dropped_frames = (
            int(
                numpy.round(
                    differences[location]
                    / (1.0e6 / true_frame_rate)
                )
            )
            - 1
        )
        first_frame, last_frame = frames[
            location : location + 2
        ]
        yield (
            first_frame,
            last_frame,
            number_of_dropped_frames,
            location,
        )


class DroppedFrames:
    """
    Used to fill in dropped frames by interpolating between closest
    available data pairs
    """

    def __init__(
        self,
        first_frame,
        last_frame,
        num_dropped_frames,
        location,
    ):
        """
        Create the DroppedFrames object

        :param first_frame: last frame before missing frame(s)
        :type: numpy.ndarray
        :param last_frame: first frame after missing frame(s)
        :type: numpy.ndarray
        :param num_dropped_frames: number of dropped frames
        :type: int
        :param location: indices where frames are missing
        :type: numpy.ndarray
        """
        self.height, self.width = first_frame.shape
        self.num_dropped_frames = num_dropped_frames
        # the frames need to be added
        self.location = location
        # Must eventually make it to use all channels.
        self.first_frame = first_frame
        self.last_frame = last_frame
        self.interpolated_frames = False

    def interpolate(self):
        """
        Produce missing frames by interpolating between first frame and last frame
        """
        diff_per_frame = (
            self.last_frame - self.first_frame
        ) / (self.num_dropped_frames + 1)

        interpolated_frames = numpy.empty(
            (
                self.num_dropped_frames,
                self.height,
                self.width,
            ),
            dtype=numpy.float32,
        )
        for frame_index in range(self.num_dropped_frames):
            interpolated_frames[frame_index] = (
                frame_index + 1
            ) * diff_per_frame + self.first_frame
        del self.first_frame
        del self.last_frame
        self.interpolated_frames = interpolated_frames
        return self


def insert_interpolated_frames(
    frames, list_of_interpolated_frames
):
    """
    Insert interpolated frames into input array as a substitute for
    dropped frames

    :param frames: array of frames from RAW channel with dropped frames
    :type: numpy.ndarray
    :param list_of_interpolated_frames: list of interpolated DroppedFrames objects
    :type: list
    :return: array of frames with dropped frames filled
    :type: numpy.ndarray
    """
    shifting_index = 1
    for interpolated_frames in list_of_interpolated_frames:
        frames = numpy.insert(
            frames,
            interpolated_frames.location
            + shifting_index,
            interpolated_frames.interpolated_frames,
            0,
        )
        shifting_index += (
            interpolated_frames.num_dropped_frames
        )
    return frames


class DarkFramesSlice:
    @staticmethod
    def threshold_method(frames, threshold=4):
        """
        Remove the dark frames at the start and end of the footage
        (Assume there are no dark frames in between the remaining frames)

        :param frames: frames of channel extracted from RAW file
        :type: numpy.ndarray
        :param threshold: threshold for average value per pixel
        :type: float

        :return: slice object
        :type: slice
        """
        temporal_means = abs(numpy.mean(frames, axis=(1, 2)))
        start, end = 0, temporal_means.shape[0]

        for i, mean in enumerate(temporal_means):
            if mean > threshold:
                start = i+1
                break
        reversed_temporal_means = numpy.flip(
            temporal_means, axis=0
        )
        del temporal_means

        for i, mean in enumerate(reversed_temporal_means):
            if mean > threshold:
                end = end - (i+1)
                break
        return slice(start, end)

    @staticmethod
    def gradient_method(
        behaviour_frames, sigma=15, spacetime=False
    ):
        """

        :param behaviour_frames:
        :param sigma:
        :param spacetime:
        :return:
        """
        if spacetime:
            means = numpy.mean(behaviour_frames, axis=1)
        else:
            means = numpy.mean(
                behaviour_frames, axis=(1, 2)
            )
        grads = numpy.gradient(means)
        mean = numpy.mean(grads)
        std = numpy.std(grads)
        start, end = 0, means.shape[0]
        threshold = mean + std * sigma
        for i, grad in enumerate(grads):
            if abs(grad) > threshold:
                start = i
                break
        reversed_grads = numpy.flip(grads, axis=0)
        del grads
        for i, grad in enumerate(reversed_grads):
            if abs(grad) < threshold:
                end = end - i
                break
        return slice(start, end)


def calculate_df_f0(frames):
    """
    Calculate df/f0, the fractional change in intensity for each pixel
    and the variance of df/f0


    :param frames: 3D array of image frames
    :type: numpy.ndarray

    :return: df/f0 with nans masked to -1
    :type: numpy.ndarray
    :return: variance in df/d0
    :type: numpy.ndarray
    """
    frames = frames.astype(numpy.float32)
    baseline = numpy.mean(frames, axis=0)
    df_f0 = numpy.divide(
        numpy.subtract(frames, baseline), baseline
    )
    del frames, baseline
    df_f0[
        numpy.where(numpy.isnan(df_f0))
    ] = -1  # Make the nans black.

    return df_f0, numpy.var(df_f0, axis=0)


def calculate_df_f0_moving(frames, n=144, axis=0):
    """
    Calculate df/f0, the fractional change in intensity for each pixel
    and the variance of df/f0 with a moving baseline (f0)


    :param n: size of moving window
    :param axis: axis to compute mean along (default 0)
    :param frames: 3D array of image frames
    :type: numpy.ndarray

    :return: df/f0 with nans masked to -1
    :type: numpy.ndarray
    """
    #
    frames = frames.astype(numpy.float32)
    baseline = numpy.cumsum(frames, axis, dtype=numpy.float32)
    baseline[n:] = baseline[n:] - baseline[:-n]
    baseline[n-1:] = baseline[n-1:] / n

    prelim = numpy.arange(n)

    # df_f0[:n-1] = df_f0[:n-1] / prelim[1:]
    baseline[:n - 1] = baseline[:n - 1] / prelim[1:, None, None]

    df_f0 = frames - baseline

    # df_f0 = numpy.divide(
    #     numpy.subtract(frames, baseline), baseline
    # )
    # del frames, baseline
    # df_f0[
    #     numpy.where(numpy.isnan(df_f0))
    # ] = -1  # Make the nans black.

    return df_f0

    # ret = np.cumsum(frames, axis=0, dtype=float)
    # ret[n:] = ret[n:] - ret[:-n]
    # ret[n - 1:] = ret[n - 1:] / n
    # return ret




class Filter:
    def __init__(
        self,
        low_freq_cutoff,
        high_freq_cutoff,
        frame_rate,
        order=4,
        rp=0.1,
    ):
        """
        Create Bandpass Filter object with frequency cutoff attributes

        :param low_freq_cutoff: first critical frequency
        :type: float>0
        :param high_freq_cutoff: second critical frequency
        :type: float>0
        :param frame_rate: frame rate of the footage
        :type: float>0
        :param order: the order of the filter
        :type: int
        :param rp: maximum ripple allowed below unity gain in the passband. Specified in decibels, as a positive number
        :type: float>0
        """
        nyq = frame_rate * 0.5
        passband = low_freq_cutoff / nyq
        stopband = high_freq_cutoff / nyq

        numerator, denominator = signal.cheby1(
            order,
            rp,
            Wn=[passband, stopband],
            btype="bandpass",
            analog=False,
        )
        self.numerator, self.denominator = (
            numerator,
            denominator,
        )

    @staticmethod
    def lfilter(numerator, denominator, data, axis=0):
        """
        Wrapper on signal.lfilter to constrain axis to time-axis
        """
        return signal.lfilter(
            numerator, denominator, data, axis
        )

    def filter(self, frames, n_jobs=None):
        """
        Use joblib to apply scipy.filter.lfilter to the data in parallel

        :param frames: frames to apply filter to
        :type: numpy.ndarray
        :param n_jobs: number of workers to utilise for parallel jobs
        :type: None or int>0

        :return: frames with applied lfilter
        :type: numpy.ndarray
        """
        n_frames, height, width = frames.shape
        frames = frames.reshape(n_frames, height * width)

        if n_jobs is None:
            n_jobs = cpu_count()
        bandpass_filter = delayed(Filter.lfilter)
        result = Parallel(n_jobs=n_jobs, verbose=0)(
            bandpass_filter(
                self.numerator,
                self.denominator,
                frames[:, s]
            )
            for s in gen_even_slices(
                frames.shape[1], n_jobs
            )
        )

        return numpy.hstack(result).reshape(
            n_frames, height, width
        )


def correct_channel_a_by_b(a, b):
    """
    Frames of channel a corrected by frames of channel b
    a/(1+b)

    :param a: Frames of channel a
    :type: numpy.ndarray
    :param b: Frames of channel b
    :type: numpy.ndarray

    :return: a/(1+b)
    """
    #return a / (1 + b)
    return a-b


def load_frames(filename, color):
    """
    Load frames of .h264/5 as color channel(s) or B&W frames as numpy array

    :param filename: path to video file
    :type: str
    :param color: one of ('red', 'green', 'blue', 'all', False), with False for B&W
    :type: str or bool

    :return: video frames
    :type: numpy.ndarray
    """
    channel = {"red": 0, "green": 1, "blue": 2, False:'', "all":''}.get(color)
    if channel is None:
        raise AttributeError(
            "Argument 'color' must be one of ('red', 'green', 'blue', 'all', False)"
        )
    cap = cv2.VideoCapture(filename)
    num_frames = 0
    first_frame = True
    while cap.isOpened:
        ret, frame = cap.read()
        if not ret:
            break
        if first_frame:
            height, width, _ = frame.shape
            first_frames = False
        num_frames +=1
    cap.release()
    if num_frames == 0:
        raise Exception(f"No frames found in 'filename'")
    if color == 'all':
        frames = numpy.zeros((num_frames, height, width, 3), dtype=numpy.uint8)
    else: 
        frames = numpy.zeros((num_frames, height, width), dtype=numpy.uint8)
    index = 0
    cap = cv2.VideoCapture(filename)
    while cap.isOpened():
        ret, frame = cap.read()
        if not ret:
            break
        if not color:
            frames[index] = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)
        elif color == 'all':
            frames[index] = cv2.cvtColor(frame, cv2.COLOR_BGR2RGB)
        else:
            frames[index] = cv2.cvtColor(frame, cv2.COLOR_BGR2RGB)[
                ..., channel
            ]
        index+=1
    cap.release()
    return frames


def video_synchronisation_indices(
    A_period, A_frames, B_period, B_frames
):
    if B_period < A_period:
        raise ValueError(
            "First video has lower frame rate than second video"
        )
    # Total length of the videos. The first frame comes in at t=0
    # so the adjustment has to be made when doing the calculation
    A_total_time, B_total_time = (
        A_period * (A_frames - 1),
        B_period * (B_frames - 1),
    )
    t_max = min(A_total_time, B_total_time)
    A = numpy.arange(0, t_max, A_period)
    B = numpy.arange(0, t_max, B_period)

    A_indices, B_indices = [], []
    for B_index, timestamp in enumerate(B):
        abs_difference = numpy.abs(A - timestamp)
        minimum = abs_difference.min()
        # If the closest match to the timestamp in B is within one
        # period of B, we can use it
        # We also keep the corresponding index of B if B is shorter than A
        if minimum <= B_period:
            A_indices.append(
                numpy.where(abs_difference == minimum)[0][0]
            )
            B_indices.append(B_index)
    return A_indices, B_indices


def downsample(array, new_shape):
    """Rebin last two dimensions of 2D or 3D array arr to shape new_shape by averaging.
       :param array: 2D or 3D, array-like
       :type: numpy.ndrray
       :param new_shape: tuple-like describing new shape of last two dimensions
       :type: tuple-like
    """
    dims = array.shape
    n_dims = len(dims)
    new_n_dims = len(new_shape)
    if new_n_dims != 2:
        raise ValueError(
            f"Argument `new_shape` should have len 2, len {new_n_dims} given instead"
        )
    if n_dims == 2:
        shape = (
            new_shape[0],
            dims[0] // new_shape[0],
            new_shape[1],
            dims[1] // new_shape[1],
        )
        return numpy.mean(
            array.reshape(shape), axis=(1, -1)
        )
    if n_dims == 3:
        shape = (
            dims[0],
            new_shape[0],
            dims[1] // new_shape[0],
            new_shape[1],
            dims[2] // new_shape[1],
        )
        return numpy.mean(
            array.reshape(shape), axis=(2, -1)
        )
    else:
        raise ValueError(
            f"Argument `array` was expected to have 2 or 3 dimensions, {n_dims} given instead"
        )


def global_signal(frames):
    """
    Calculate: 1. 'mean' of `frames`
               2. 'beta', the matrix product of the Moore-Penrose pseudo-inverse of 'mean', and `frames`
               3. 'globalsignal', the matrix product of `mean` and `beta`

    :param frames: array-like
    :type: numpy.ndarray
    :return: globalsignal, mean, beta
    """
    mean_g = numpy.mean(frames, axis=1)
    g_plus = numpy.squeeze(numpy.linalg.pinv([mean_g]))
    mean_g = numpy.expand_dims(mean_g, axis=1)
    g_plus = numpy.expand_dims(g_plus, axis=0)
    beta_g = numpy.matmul(g_plus, frames)
    globalsignal = numpy.matmul(mean_g, beta_g)

    return globalsignal, mean_g, beta_g


def draw_mask(frame):
    plt.imshow(frame, cmap='gray', vmin=0, vmax=255)
    roi = RoiPoly(color='r')
    mask_left = roi.get_mask(frame)
    mask_right = roi.get_mask(frame)

    mask = numpy.logical_or(
        mask_left.get_mask(frame),
        mask_right.get_mask(frame)
    )
    return mask
