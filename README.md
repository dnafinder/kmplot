[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=dnafinder/kmplot&file=kmplot.m)

# kmplot

## 📘 Overview
kmplot is a MATLAB function that computes and plots the Kaplan–Meier (product-limit) estimator of the survival function from right-censored time-to-event data.

Kaplan–Meier analysis is a cornerstone in survival analysis and medical statistics. It allows estimation of the survival function S(t) without assuming any specific parametric distribution for event times and can properly handle censored observations such as:

- patients lost to follow-up
- patients still event-free at the end of the study
- incomplete observations due to study termination

kmplot:
- constructs the Kaplan–Meier step function S(t)
- accounts for censored observations
- computes Greenwood-type standard errors and confidence intervals
- estimates the median survival time (if it exists)
- provides a simple estimate of the average hazard rate
- optionally returns all underlying data structures for further analysis

## ✨ Features
- Handles right-censored survival data in a simple [time, censorFlag] format
- Accepts:
  - a vector of times (all uncensored), or
  - an N-by-2 matrix [time, censorFlag]
- Produces Kaplan–Meier survival curves with:
  - main survival function S(t)
  - approximate (1 – α)·100% confidence bands
  - optional marking of censored observations
  - optional display of median survival time
- Computes:
  - table of event counts and numbers at risk
  - median survival time (via interpolation if needed)
  - average hazard rate estimate
- Can be used:
  - interactively to visualize survival curves
  - programmatically to return data for further inferential procedures (e.g. log-rank tests)

## 📥 Installation
1. Download or clone the repository:
   https://github.com/dnafinder/kmplot

2. Add the folder to your MATLAB path, for example:
      addpath('path_to_kmplot_folder')

3. Verify that MATLAB can see the function:
      which kmplot

If MATLAB returns the path to kmplot.m, the installation is successful.

## ⚙️ Requirements
- MATLAB (any recent version)
- Statistics and Machine Learning Toolbox (required for:
  - tabulate
  - basic statistical operations)

No additional toolboxes are strictly required beyond the above.

## 📈 Usage

Basic usage with an N-by-2 matrix (time, censorFlag):

    % x(:,1) = survival times
    % x(:,2) = 0 (uncensored) or 1 (censored)
    kmplot(x);

Using a vector of event times (all uncensored):

    t = [2; 3; 5; 7; 9; 10];
    kmplot(t);

Specify a different significance level alpha (for confidence intervals):

    kmplot(x, 0.01);   % 99% confidence bands

Control plotting of censored observations:

    % cflag = 0: spread censored marks along the horizontal segment
    % cflag = 1: place censored marks at the censoring time
    kmplot(x, 0.05, 0);
    kmplot(x, 0.05, 1);

Get outputs without plotting (for internal use, e.g. by LOGRANK):

    % flag = 0 turns off plotting and some display elements
    [table1, table12, t1, T1, xcg, ycg, lambda] = kmplot(x, 0.05, 0, 0);

## 🔢 Inputs

kmplot(x)  
kmplot(x, alpha)  
kmplot(x, alpha, cflag)  
kmplot(x, alpha, cflag, flag)

- x  
  - Type:
    - numeric vector (N×1), or
    - numeric matrix (N×2)
  - Description:
    - If x is a vector:
      - interpreted as survival times
      - all observations are considered uncensored
      - internally converted to [x(:) zeros(N,1)]
    - If x is N-by-2:
      - x(:,1) = survival times
      - x(:,2) = censoring indicator (0 = uncensored, 1 = censored)

- alpha (optional)  
  - Type: scalar numeric  
  - Default: 0.05  
  - Description:
    - significance level for confidence intervals
    - the function constructs (1 – alpha)·100% confidence bands around S(t)

- cflag (optional)  
  - Type: scalar numeric (0 or 1)  
  - Default: 0  
  - Description:
    - controls how censored observations are plotted:
      - 0 → censored marks are spread along the horizontal segment to avoid overlap
      - 1 → censored marks are placed at the actual censoring time on the curve

- flag (optional)  
  - Type: scalar numeric (0 or 1)  
  - Default: 1  
  - Description:
    - controls full plotting behaviour:
      - 1 → standard usage: produces plots, confidence intervals, median time, etc.
      - 0 → internal mode: used by other functions (e.g. LOGRANK) to obtain survival
        estimates and tables without generating figures

## 📤 Outputs

When called with output arguments, kmplot returns:

    [table1, table12, t1, T1, xcg, ycg, lambda] = kmplot(...);

- table1  
  - N1-by-3 numeric matrix containing:
    - column 1: distinct event times
    - column 2: number of events (deaths) at each time (after censoring is accounted for)
    - column 3: number at risk at the beginning of each interval

- table12  
  - Nc-by-3 numeric matrix derived from tabulate on censored times:
    - column 1: distinct censoring times
    - column 2: number of censored observations at each time
    - column 3: percentage (ignored in computations but useful for reference)

- t1  
  - Time points used for the step function, including the initial time 0 and the last extended time (if needed when censored times exceed the last event time).

- T1  
  - Kaplan–Meier survival estimates S(t1).

- xcg, ycg  
  - Coordinates of censored marks for plotting:
    - xcg: x-positions of censored observations
    - ycg: y-positions (survival level at which censored marks are drawn)

- lambda  
  - Average hazard rate estimate computed from the log-survival decrements.

If you request fewer than seven outputs, kmplot simply returns the first k elements in the order above.

![](https://github.com/dnafinder/kmplot/blob/master/kmplot.jpg)

## 🧠 Interpretation

- The step function T1 vs t1 represents the estimated survival function S(t).
- Confidence bands are computed as:
  - S(t) ± cv · SE(t)
  - where SE(t) is obtained via a Greenwood-type formula and cv is based on the standard normal approximation.
- The median survival time is:
  - either the time at which S(t) equals 0.5,
  - or, more often, obtained via linear interpolation between the two steps that bracket 0.5.
- The hazard rate lambda is a rough average hazard over the observed time span:
  - based on finite differences of log-survival
  - values should be interpreted with caution for highly non-constant hazards.

Good survival curves typically show:
- a non-increasing step function S(t)
- confidence bands that behave consistently (never below 0 or above 1, and trim where appropriate)
- a stable median time (if S(t) falls below 0.5 within the observed range)

## 📌 Notes

- Censoring is assumed to be independent and non-informative with respect to the event process.
- When x is provided as a vector, all data are treated as uncensored; if censoring is present, you should use the N-by-2 matrix format with explicit censoring flags.
- The optional flag argument is mainly intended for internal usage by log-rank or other survival-related functions, so ordinary users typically do not need to modify it.
- Censored marks can be visually crowded when many observations are censored at the same time; cflag = 0 mitigates this by spreading them horizontally.

## 🧾 Citation

If you use kmplot in research, analysis, or publications, please cite:

Cardillo G. (2008). KMPLOT: Kaplan–Meier estimation of the survival function.  
Available at: https://github.com/dnafinder/kmplot

## 👤 Author

Giuseppe Cardillo  
Email: giuseppe.cardillo.75@gmail.com  
GitHub: https://github.com/dnafinder

## 📄 License

The code is provided as-is, without any explicit warranty.  
Please refer to the repository for licensing details if a LICENSE file is present.  
kmplot is distributed under the terms specified in the LICENSE file:
https://github.com/dnafinder/kmplot
