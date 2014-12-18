This module is to initialize tracker. I wrote this for a part of my MS thesis "ROI based Parametric Video Coding Scheme".
This code uses first two frames of a video, does motion segmentation using phase-correlation based motion estmation and color segmentation. Then it combines the two segmentation and other heuristics to identify potential ROI candidate regions. Finally using greedy method, the ROI candidate regions are clustered to form final (single/or multiple) ROIs. 

usage:
tracker_init.m is the start-up routine

Please cite my MS thesis if you use this code.

Subarna Tripathi, MS(R) thesis, IIT Delhi, 2011, "ROI based Parametric Video Coding Scheme"