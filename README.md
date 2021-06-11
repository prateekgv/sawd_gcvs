## Gait Cycle Validation and Segmentation using Inertial Sensors

### Background
In the last decade, gait analysis has moved away from equipment-intensive, laboratory-based analyses toward the use of wearable sensors. Several gait segmentation methods for ambulatory gait analysis using wearable technology have been developed, with varying success. In the case of inertial sensors, the existing methods are based on thresholding sensor measurements, template matching via dynamic time warping, or machine learning based methods such as hidden Markov models. Across all these methods, gyroscope measurements in the sagittal plane are the best choice for gait segmentation because the measurements contain typical time-series patterns such as "valleys," "peaks," and "plateaus." These patterns are respectively, toe-off and heel-strike, mid-swing, and midstance. While thresholding methods work well in practice, these methods do not have a mechanism to verify and validate the observable patterns of the gyroscope signal in the sagittal plane, which leads to low precision. As an alternative, dynamic time warping is a pattern matching method that computes a similarity measure between a template representing a valid gait cycle and an input sequence. However, the threshold required to validate an input sequence as a gait cycle is not fixed and exhibits a large variance. Yet another alternative, machine learning methods, requires large samples of training data across many participants and manual segmentation of gait phases, depending on the granularity of the task, to label the training data. Moreover, it is computationally expensive to learn the parameters of the machine learning task.

<p align="center">
  <img width="800" height="200" src="https://github.com/prateekgv/sawd_gcvs/blob/master/matlab/figures/left_seg_gait_cycle_1.png">
</p>


### Our research
* We propose a modular approach that employs pattern matching and thresholding, to validate and detect three gait events in the gait cycle, namely midstance, toe-off, and heel-strike (see Fig. 1).
* To identify the data as stationary or moving, in the detection module, we first use physical models that describe zero-velocity events or stationary events of the sensor data obtained from a foot-mounted inertial navigation system.
* Next, in the gait cycle validation and segmentation module, we simultaneously combine linear time-invariant filters, wavelets, and sparsity-based methods to extract a discrete wavelet transform (DWT) coefficient vector as a sparse representation of the moving segment of the gyroscope measurements in the sagittal plane.
* We generate a template of the DWT coefficient vector by taking the average of the DWT coefficient vectors obtained using the SAWD algorithm for all valid gait cycles of a given trial.
* Thereafter, to validate any moving segment as a gait cycle, we compute the root-mean-square error between the generated template and the sparse representation of the moving segment of the gyroscope data in the sagittal plane, obtained using the SAWD algorithm.
* Our proposed method demonstrates an average F1 score of 87.78% across all groups for a fixed sampling rate, and an average F1 score of 92.44% across all Parkinson disease participants for a variable sampling rate.

## License
The source code is released under the [MIT](LICENSE.md) license.

## Citation
If you find this work interesting and useful, please cite the following publication:
> G. V. Prateek, P. Mazzoni, G. M. Earhart and A. Nehorai, "Gait Cycle Validation and Segmentation using Inertial Sensors," in IEEE Transactions on Biomedical Engineering. doi: 10.1109/TBME.2019.2955423
