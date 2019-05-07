SignalBox : Read Me
========================
_A set of tools for signal analysis and design in the Time and Frequency domains._

Extensions to the __Signal__ class include:


### Cosine basis

| Method | Description |
| ----------- | ----------- |
| `*cosineFill` | Fill a Signal of the given size with a sum of cosines at the given amplitudes and phases. |
| `cosineFill` | A sum of cosines with the given amplitudes and phases, returned in place. |
| `cosineFill2` | A sum of cosines with the given harmonics, amplitudes and phases, returned in place. |


### Wave packet

| Method | Description |
| ----------- | ----------- |
| `*readPacket` | Read a Hann enveloped wave packet from a soundfile. |
| `*packetFill` | Fill a Signal of the given size with a sum of Hann enveloped cosines at the given amplitudes and phases. |
| `packetFill` | A sum of Hann enveloped cosines at the given amplitudes and phases, returned in place. |
| `packetFill2` | A sum of Hann enveloped cosines with the given harmonics, amplitudes and phases, returned in place. |


### Reading & writing

| Method | Description |
| ----------- | ----------- |
| `*read` | Read a soundfile. |
| `*readWave` | Read a periodic waveform from a soundfile. |
| `write` | Writes the entire Signal to a soundfile. |
| `writeFile` | Writes the entire Signal to a file, where every chunk of four bytes is a 32-bit floating point sample. |


### Fourier Transform

| Method | Description |
| ----------- | ----------- |
| `*rfftCosTable` | Fill a Signal with the cosine table needed by the Real-FFT methods. |
| `*rfftTwoCosTable` | Fill a Signal with the cosine table needed by the Real-FFT-Two methods. |
| `fftToRfft` | Return a complex Real-FFT spectrum from a complex FFT spectrum. |
| `rfftToFft` | Return a complex FFT spectrum from a complex Real-FFT spectrum. |
| `rfft` | Perform a Real-FFT on a real signal. |
| `irfft` | Perform an inverse Real-FFT on a real and imaginary signal. |
| `rfftTwo` | Perform a Real-FFT on two real signals. |
| `irfftTwo` | Perform an inverse Real-FFT on two real and imaginary signals. |
| `dft` | Perform a DFT on a real and imaginary signal. |
| `idft` | Perform an inverse DFT on a real and imaginary signal. |
| `dftZoom` | Perform a Zoom DFT on a real and imaginary signal. |
| `rdftZoom` | Perform a Zoom DFT on a real signal. |
| `goertzel` | Return individual terms of a DFT on a real and imaginary signal. |
| `rgoertzel` | Return individual terms of a DFT on a real signal. |


### Chirp z-Transform

| Method | Description |
| ----------- | ----------- |
| `czt` | Perform a Chirp z-Transform on a real and imaginary signal. |
| `rczt` | Perform a Chirp z-Transform on a real signal. |


### Cepstrum

| Method | Description |
| ----------- | ----------- |
| `rceps` | Return the real part of the cepstrum of a real signal. |
| `irceps` | Return the real part of the inverse of the real part of the cepstrum. |


### Analytic & Hilbert

| Method | Description |
| ----------- | ----------- |
| `*analyticFill` | Return a complex analytic Signal of the given size with a sum of cosines and a sum of sines at the given amplitudes and phases. |
| `*hilbert` | Hilbert Transform: Return complex Hilbert Transform coefficients. |
| `analytic` | Hilbert Transform: Return a complex analytic signal from a real signal. |


### Rotation & Phase

| Method | Description |
| ----------- | ----------- |
| `rotateWave` | Rotate the Signal by a value in radians, in place. |
| `rotatePhase` | Rotate the phase of Signal by a value in radians, in place. |
| `linearPhase` | Return a linear phase kernel, preserving the magnitude, in place. |
| `minimumPhase` | Return a minimum phase kernel, preserving the magnitude, in place. |
| `even` | Return the even part of a signal. |
| `odd` | Return the odd part of a signal. |

### RMS & Magnitude

| Method | Description |
| ----------- | ----------- |
| `rms` | Return the RMS of a signal. |
| `peakMagnitude` | Return the peak of the frequency domain magnitude of a signal. |
| `normalizeMagnitude` | Normalize the Signal in place such that the maximum magnitude in the frequency domain is 1. |



&nbsp;

&nbsp;

Installing
==========

Distributed via
[SC3 SignalBox](https://github.com/ambisonictoolkit/SignalBox).

Start by reviewing the Quark installation instructions
[found here](https://github.com/supercollider-quarks/quarks#installing). See
also [Using Quarks](http://doc.sccode.org/Guides/UsingQuarks.html).

With [git](https://git-scm.com/) installed, you can easily install the
[SC3 SignalBox](https://github.com/ambisonictoolkit/SignalBox)
directly by running the following line of code in SuperCollider:

    Quarks.install("https://github.com/ambisonictoolkit/SignalBox.git");



Feedback and Bug Reports
========================

Known issues are logged at
[GitHub](https://github.com/ambisonictoolkit/SignalBox/issues).

&nbsp;


List of Changes
---------------

Version 0.1.0

* First Public Release.

&nbsp;

Credits
=======

Joseph Anderson 2019.

* J Anderson : [[e-mail]](mailto:joanders[at]uw.edu)

&nbsp;

The development of SignalBox for SuperCollider3 is
supported by
[DXARTS, Center for Digital Arts and Experimental Media](https://dxarts.washington.edu/).

&nbsp;


Contributors
------------

Version 0.1.0
*  Joseph Anderson (@joslloand)
