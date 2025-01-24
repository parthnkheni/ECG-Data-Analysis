# ECG Data Analysis

This repository contains the MATLAB implementation of an **ECG Data Analysis** project, which processes and analyzes electrocardiogram (ECG) signals to extract meaningful insights and features. The project focuses on signal preprocessing, peak detection, and deriving metrics such as heart rate variability (HRV) and other cardiac indicators.

---

## Overview

Electrocardiography (ECG) is a vital tool for monitoring heart health. This project leverages signal processing techniques to analyze raw ECG data, aiming to identify key features like R-peaks, compute heart rate, and provide a foundation for more advanced cardiac diagnostics.

---

## Features

- **Signal Preprocessing**: Noise reduction using filtering techniques (e.g., low-pass and high-pass filters).
- **R-Peak Detection**: Accurate identification of R-peaks in the ECG signal.
- **Heart Rate Calculation**: Derives heart rate based on detected R-R intervals.
- **Heart Rate Variability (HRV) Analysis**: Calculates HRV metrics to assess cardiac health.
- **Visualizations**: Provides plots of the raw and processed ECG signals, along with detected peaks.

---

## How It Works

1. **Data Input**: Import raw ECG data from a file or data source.
2. **Preprocessing**:
   - Remove noise using bandpass filtering.
   - Normalize the ECG signal.
3. **R-Peak Detection**:
   - Apply an algorithm (e.g., Pan-Tompkins or similar) to detect R-peaks in the ECG signal.
4. **Feature Extraction**:
   - Compute heart rate and other metrics based on R-R intervals.
5. **Visualization**:
   - Plot raw and filtered signals, marking detected R-peaks.
   - Generate HRV plots if applicable.

---

## Requirements

- **MATLAB** (R2021a or newer recommended)
- Signal Processing Toolbox

---

## Setup & Usage

1. Clone this repository:
   ```bash
   git clone https://github.com/parthnkheni/ecg-data-analysis.git
   ```
2. Open the `projectfinalv3.mlx` file in MATLAB.
3. Follow the instructions in the script to load your ECG data.
4. Run the script to preprocess, analyze, and visualize the data.

---

## Example Results

The following outputs can be expected:
- Processed ECG signal with noise removed.
- Marked R-peaks on the ECG waveform.
- Calculated heart rate and HRV metrics.
- Visualizations for better understanding of the data.

---

## Future Improvements

1. Extend the algorithm to classify arrhythmias using machine learning.
2. Improve preprocessing to handle noisy or low-quality signals.
3. Add support for real-time ECG analysis with wearable devices.
