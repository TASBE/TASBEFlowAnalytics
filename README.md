# [TASBE Flow Analytics](https://tasbe.github.io/TASBE/)

TASBE Flow Analytics is a flow cytometry analysis package. 
Documentation, templates, and the Excel user interface are stored in the [TASBE Tutorial project](https://github.com/TASBE/TASBEFlowAnalytics-tutorial) and in the [TASBE project website](https://tasbe.github.io/).

## Build status:
[![Build Status](https://travis-ci.org/TASBE/TASBEFlowAnalytics.svg?branch=master)](https://travis-ci.org/TASBE/TASBEFlowAnalytics) for this project

[![Build Status](https://travis-ci.org/TASBE/TASBEFlowAnalytics-Tutorial.svg?branch=master)](https://travis-ci.org/TASBE/TASBEFlowAnalytics-Tutorial) for tutorial and templates

## Features

- Runs on both Matlab and Octave, and from Python/Jupyter via oct2py
- User-friendly Excel interface
- Flow cytometry analysis
- Plotting and comparison templates for many experiments, including comparative analysis and parametric analysis
- Output of "point clouds", histograms, and key statistics
- Unit conversion to ERF from multiple channels
- Compensation for autofluorescence and spectral overlap
- Automation-assisted gating
- Automatic "sanity check" warnings and detailed diagnostic reports
- Scriptable for high-throughput analysis
- Distributed under a permissive free and open license

**If you make use of the TASBE Flow Analytics package, please cite
the following two publications:**

* Jacob Beal, "Bridging the Gap: A Roadmap to Breaking the Biological
  Design Barrier," Frontiers in Bioengineering and Biotechnology,
  2:87. doi:10.3389/fbioe.2014.00087, January 2015.

* Jacob Beal, Ron Weiss, Fusun Yaman, Noah Davidsohn, and Aaron Adler,
  "A Method for Fast, High- Precision Characterization of Synthetic
  Biology Devices," MIT CSAIL Tech Report 2012-008, April 2012. 
  http://hdl.handle.net/1721.1/69973

## Installation

- Dependencies:

  - If you are using Matlab, TASBE prefers to have the `stats` package installed (though it can use alternate code) 
  - If you are using Octave, TASBE depends on the `io` package (recommended version 2.4.10)
  
     - This can typically be installed with `octave --no-gui --quiet --eval "pkg install -forge io"`

- Installation on Apple OSX or GNU/Linux using the shell:

    ```bash
    git clone https://github.com/TASBE/TASBEFlowAnalytics.git
    cd TASBEFlowAnalytics
    make install
    ```
    This will add the TASBEFlowAnalyics directory to the Matlab and/or GNU Octave searchpath. If both Matlab and GNU Octave are available on your machine, it will install TASBEFlowAnalyics for both.

- Installation on Windows (or manual installation on Mac/Linux):
  - Download the package from [GitHub](https://github.com/TASBE/TASBEFlowAnalyics)
  - Start Matlab or Octave
  - Go to the ``TASBEFlowAnalytics`` directory
  - Run the installation command by hand:
  
      ```
    tasbe_set_path(); savepath();
    ```
- **Optional:** If you want to run the test files, install [MOxUnit](https://github.com/MOxUnit/MOxUnit), and install the tutorial/example package from [GitHub](https://github.com/TASBE/TASBEFlowAnalyics-Tutorial) in a sibling directory to TASBEFlowAnalytics named `TASBEFlowAnalyics-Tutorial`.

## Usage

In use of this package, you will typically want to split your
processing into three stages:

- Creation of a ColorModel that translates raw FCS to comparable unit data
- Using a ColorModel for batch processing of experimental data
- Comparison and plotting of the results of batch processing

Example files are provided in the [TASBE Tutorial](https://github.com/TASBE/TASBEFlowAnalytics-tutorial) that show how these stages typically work.
