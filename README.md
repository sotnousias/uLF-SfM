# uLF-SfM

Code for our CVPR 2019 paper: Large-scale, metric structure from motion for unordered light fields.

Incremental 3D reconstruction for light field systems (e.g., multi-camera arrays or microlens based cameras).

An example of the reconstuction process is on main.m.

Key steps are the following.

## Features

- Light field feature extraction.
- Pairwise feature matching.
- Multi-image matching with visibility matrix calculation.
- Initial frame selection based on geometric verification.

## Structure from Motion

- Estimate relative pose of the initial camera pair.
- Traingulate 3D points from the initial camera pair.
- Select next best view, estimate the pose and incrementally triangulate points.
- Occasioanlly run light field bundle adjustment to refine camera poses and 3D points.

## Requirements

- MATLAB
- Please install vlfeat prior to running the code: https://www.vlfeat.org/.

## Getting Started

Before running the script, ensure that you have all the necessary input files placed in the `./txtFiles/` directory.

1. **Input Files**:
   - `K.txt`: Camera intrinsic parameters for the central sub-aperture images.
   - `subApertureRelPoses.txt`: Relative poses of sub-aperture images. Each row has 7 parameters, the first 4 parameters are the quaternion representation of the rotation matrix and the last 3 the translation vector.



## Running the Script

1. Open MATLAB.
2. Navigate to the directory containing `main.m`.
3. Run the script by typing `main` in the MATLAB command window.

## Output

The script will output the following:
- 3D reconstructed points of the scene.
- Camera poses.
- Feature matches and visibility matrices.


## Notes

- The script includes warnings suppression. Modify this behavior if needed.
  

For detailed information about each function and parameter, refer to the inline comments in the `main.m` file.


