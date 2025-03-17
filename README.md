# X-Seg: Scalable Bone Segmentation Algorithm for Hand X-rays

This MATLAB code was developed as part of the final project for Stanford EE 368: Digital Image Processing, Winter 2025, taught by Prof. David Stork.

## Overview

X-Seg is an automated, scalable image processing algorithm designed to segment bones in pediatric hand X-rays and extract key physiological parameters. The algorithm was used to construct a preliminary bone age growth chart based on over 10,000 X-rays, providing a foundation for future population-representative growth charts.

Unlike traditional bone age assessment methods reliant on outdated reference atlases, X-Seg offers an objective, transparent, and efficient alternative that requires no manual labeling. The algorithm was validated using the RSNA 2017 Bone Age dataset (available: https://www.rsna.org/rsnai/ai-image-challenge/rsna-pediatric-bone-age-challenge-2017) and achieved a Dice similarity coefficient of 0.909, indicating strong correspondence with manual segmentations.

## Running the code

The `src` directory contains all code files. `pipeline.m` can be used to perform segmentations.

Before running code, the following settings should be configured:
- **sample:** true if you want to segment a random X-ray, otherwise false
- **make_plot:** true if you want to plot the resulting segmentation, otherwise false
- **p:** the number of images you want to segment; must be less than 5 if `make_plot` is true
- note: all images must be in the same directory as `pipeline.m`, saved under the subdirectory `boneage-training-dataset/`. Metadata should be in the same directory in a file named `boneage-training-dataset.csv`.

## Acknowledgments

- Prof. David Stork, Stanford University
- Rhea Prem, EE 368 Course Assistant
- RSNA 2017 Bone Age Dataset contributors

