/*
	Copyright Joseph Anderson, 2019
		J Anderson	joanders@uw.edu

	This file is part of the SignalBox quark for SuperCollider 3 and is free software:
	you can redistribute it and/or modify it under the terms of the GNU General
	Public License as published by the Free Software Foundation, either version 3
	of the License, or (at your option) any later version.

	This code is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
	without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
	See the GNU General Public License for more details.
	<http://www.gnu.org/licenses/>.
*/


//---------------------------------------------------------------------
//
// 	Extension: Signal
//
//---------------------------------------------------------------------


+ Signal {

	/* misc */

	// maxdb { arg aNumber; ^this.as(Array).maxdb(aNumber).as(Signal) }
	// mindb { arg aNumber; ^this.as(Array).mindb(aNumber).as(Signal) }
	// clipdb2 { arg aNumber; ^this.as(Array).clipdb2(aNumber).as(Signal) }
	// threshdb { arg aNumber; ^this.as(Array).threshdb(aNumber).as(Signal) }
	// clipdb { arg lo, hi; ^this.as(Array).clipdb(lo, hi).as(Signal) }

	flip {
		^this.deepCopy.reverse.rotate(1)
	}

	complex { arg imag;
		^Complex.new(this, imag)
	}

	/* zeros */

	*zeroFill { arg size;
		^Signal.newClear(size)
	}

	/* wrapping */
	wrapExtend { arg size;
		^this.as(Array).wrapExtend(size).as(Signal)
	}

	/* sine */

	/* cosine */
	*cosineFill { arg size, amplitudes, phases;
		^Signal.newClear(size).cosineFill(amplitudes, phases).normalize
	}

	addCosine { arg harmonicNumber = 1, amplitude = 1.0, phase = 0.0;
		^this.addSine(harmonicNumber, amplitude, phase + 0.5pi)
	}

	cosineFill { arg amplitudes, phases;
		this.fill(0.0);
		if (phases.isNil, { phases = #[0]; });
		amplitudes.do({ arg amp, i; this.addCosine(i+1, amp, phases @@ i) });
	}

	cosineFill2 { arg list;
		this.fill(0.0);
		list.do({ arg item, i;
			var harm, amp, phase;
			# harm, amp, phase = item;
			this.addCosine(harm, amp ? 1.0, phase ? 0.0);
		});
	}

	/* read waveform */

	*readWave { arg size, path, frame = 0, channel = 0, freq = 440.0, alpha = 3, winScale = 1;
		var soundFile, soundFrames;
		var startFrame, endFrame;
		var k, winSize, perSize, readSize, overDubIndex;
		var soundData, window, foldedWin, signal, foldedSig;

		if(File.exists(path), {

			// open...
			soundFile = SoundFile.new(path);
			soundFile.openRead;

			// calcs...
			perSize = (soundFile.sampleRate / freq).asInteger;  // period size
			k = (1 + alpha.squared).sqrt;  // folding frequency
			winSize = (2.pow(winScale) * k * perSize).asInteger;  // window size
			startFrame = (frame - (winSize/2)).asInteger;  // window is centered
			endFrame = startFrame + winSize;  // window is centered
			soundFrames = (soundFile.numFrames / soundFile.numChannels).asInteger;

			case
			// within bounds
			{ (startFrame >= 0) && (endFrame <= soundFrames) } {
				readSize = winSize;
				overDubIndex = 0;
			}
			// before soundFile start
			{ (startFrame < 0) && (endFrame <= soundFrames) } {
				readSize = startFrame + winSize;  // negagive start frame
				overDubIndex = startFrame.neg;
			}
			// after soundFile end
			{ (startFrame >= 0) && (endFrame > soundFrames) } {
				readSize = startFrame + winSize - soundFrames;
				overDubIndex = 0;
			}
			// before AND after soundFile end
			{ (startFrame < 0) && (endFrame > soundFrames) } {
				readSize = soundFrames;  // negagive start frame
				overDubIndex = startFrame.neg;
			};

			// read
			soundData = FloatArray.newClear(readSize * soundFile.numChannels);  // interleaved
			(startFrame >= 0).if({
				soundFile.seek(startFrame).readData(soundData);
				}, {
				soundFile.seek(0).readData(soundData);
			});

			// multi-channel??
			(soundFile.numChannels > 1).if({
				soundData = soundData.clump(soundFile.numChannels).flop.at(channel)
			});
			// ... close
			soundFile.close;

			// window & folded window
			window = Signal.kaiserWindow(winSize, a: alpha);
			foldedWin = Signal.newClear(perSize);
			window.clump(perSize).do({ arg item;
				foldedWin.overDub(item.as(Signal))
			});

			// signal & folded signal
			signal = window * Signal.newClear(winSize).overDub(soundData.as(Signal), overDubIndex);
			foldedSig = Signal.newClear(perSize);
			signal.clump(perSize).do({ arg item;
				foldedSig.overDub(item.as(Signal))
			});
			foldedSig = foldedSig * foldedWin.reciprocal;  // normalize the window!

			// normalize phase
			foldedSig = foldedSig.rotateWave(foldedSig.rgoertzel(1).phase.neg);

			// resample to size...
			(size == nil).not.if({
				foldedSig = foldedSig.resize(size);
			});

			// remove dc... & normalize scale
			foldedSig.discardDC;
			foldedSig.normalize;

			^foldedSig;
		}, {
			Error("No sound file at:" + path).throw;
		})
	}

	/* read & write */

	writeFile { arg path;
		var file;
		file = File.new(path, "wb");
		if (file.notNil, {
			file.write(this);
			file.close
		})
	}

	writeSoundFile { arg path, headerFormat = "WAV", sampleFormat = "float", sampleRate = 44100;
		var soundFile = SoundFile.new
		.headerFormat_(headerFormat)
		.sampleFormat_(sampleFormat)
		.sampleRate_(sampleRate);

		soundFile.openWrite(path);
		soundFile.writeData(this);
		soundFile.close
	}



	/* magnitude & envelope */

	rms { ^this.squared.mean.sqrt }

	peakMagnitude { arg oversample = 1;
		var osSize, osThis;

		osSize = (oversample.isInteger || this.size.isPowerOfTwo).if({
			(oversample * this.size).nextPowerOfTwo
		}, {
			(oversample * this.size).floor
		}).asInteger;
		osThis = this ++ Signal.newClear(osSize - this.size);

		^(osSize.isPowerOfTwo).if({  // rfft
			var cosTable = Signal.rfftCosTable(osSize/2 + 1);
			osThis.rfft(cosTable).magnitude.maxItem
		}, {  // czt via dft
			osThis.dft(Signal.newClear(osSize)).magnitude.maxItem
		})
	}

	// Normalize the Signal in place such that the maximum absolute peak value is 1.
	normalizeMagnitude { arg oversample = 1;
		^this.scale(this.peakMagnitude(oversample).reciprocal)
	}

	/* DC */
	discardDC { ^this.offset(this.mean.neg) }

	/* rotation & phase */

	// NOTE: match Hilbert phase rotation
	rotateWave { arg phase;
		var precision = 0.1;  // corresponds to integer
		var n = phase / 2pi.neg * this.size;
		var intN = n.asInteger;

		^this.overWrite(
			n.equalWithPrecision(intN, precision).if({  // direct
				this.rotate(n)
			}, {
				(this.size.isPowerOfTwo).if({  // fft
					var cosTable = Signal.rfftCosTable(this.size / 2 + 1);
					var complex = this.rfft(cosTable);
					var rotatedComplex;
					var start = 0, step, phaseOffset;

					step = phase;
					phaseOffset = Array.series(this.size / 2 + 1, start, step);  // positive freqs only

					rotatedComplex = complex.rotate(phaseOffset);  // converts to Array
					rotatedComplex.real.as(Signal).irfft(rotatedComplex.imag.as(Signal), cosTable)
				}, {  // dft
					var complex = this.dft(Signal.newClear(this.size));
					var rotatedComplex;
					var start, step, phaseOffset;

					step = phase;

					this.size.even.if({
						start = step.neg * this.size / 2;  // start with negative freqs
						phaseOffset = Array.series(this.size, start, step);
						phaseOffset = phaseOffset.rotate((this.size / 2).asInteger)  // rotate
					}, {
						start = step.neg * (this.size-1) / 2;  // start with negative freqs
						phaseOffset = Array.series(this.size, start, step);
						phaseOffset = phaseOffset.rotate(((this.size+1) / 2).asInteger)  // rotate
					});
					rotatedComplex = complex.rotate(phaseOffset);    // converts to Array
					rotatedComplex.real.as(Signal).idft(rotatedComplex.imag.as(Signal)).real
				})
			})
		)
	}

	rotatePhase { arg phase;
		var complex;

		^this.overWrite(
			(this.size.isPowerOfTwo).if({  // rfft
				var cosTable = Signal.rfftCosTable(this.size/2 + 1);

				complex = this.rfft(cosTable);  // real fft
				complex = Complex.new(phase.cos, phase.sin) * complex;  // rotate

				complex.real.irfft(complex.imag, cosTable)  // irfft
			}, {  // czt via analytic
				complex = this.analytic;

				// equivalent to: (Complex.new(phase.cos, phase.sin) * complex).real
				(phase.cos * complex.real) + (phase.sin.neg * complex.imag)
			})
		)
	}

	linearPhase { arg sym = false;
		var magnitude, phase, complex;

		^this.overWrite(
			(this.size.isPowerOfTwo).if({  // rfft
				var cosTable = Signal.rfftCosTable(this.size/2 + 1);
				var rcomplex, rcomplexMin;

				rcomplex = this.rfft(cosTable);  // real fft
				magnitude = rcomplex.real.rfftToFft(rcomplex.imag).magnitude;  // mirror spectrum & magnitude
				phase = magnitude.linearPhase(sym);  // linear phase

				complex = Polar.new(magnitude.as(Signal), phase.as(Signal)).asComplex;  // linear phase spectrum
				rcomplexMin = complex.real.fftToRfft(complex.imag);  // discard negative freqs

				rcomplexMin.real.irfft(rcomplexMin.imag, cosTable)  // irfft
			}, {  // czt via dft
				magnitude = this.dft(Signal.newClear(this.size)).magnitude;  // magnitude
				phase = magnitude.linearPhase(sym);  // linear phase
				complex = Polar.new(magnitude.as(Signal), phase.as(Signal)).asComplex;  // linear phase spectrum

				complex.real.idft(complex.imag).real  // idft
			})
		)
	}

	// Hilbert minimum phase
	minimumPhase { arg mindb = -120.0, oversample = 1;
		var magnitude, phase, complex;
		var osSize, osThis;

		osSize = (oversample.isInteger || this.size.isPowerOfTwo).if({
			(oversample * this.size).nextPowerOfTwo
		}, {
			(oversample * this.size).floor
		}).asInteger;
		osThis = this ++ Signal.newClear(osSize - this.size);

		^this.overWrite(
			(osSize.isPowerOfTwo).if({  // rfft
				var cosTable = Signal.rfftCosTable(osSize/2 + 1);
				var rcomplex, rcomplexMin;

				rcomplex = osThis.rfft(cosTable);  // real fft
				magnitude = rcomplex.real.rfftToFft(rcomplex.imag).magnitude;  // mirror spectrum & magnitude
				phase = magnitude.minimumPhase(mindb);  // minimum phase

				complex = Polar.new(magnitude.as(Signal), phase.as(Signal)).asComplex;  // minimum phase spectrum
				rcomplexMin = complex.real.fftToRfft(complex.imag);  // discard negative freqs

				rcomplexMin.real.irfft(rcomplexMin.imag, cosTable).keep(this.size)  // irfft
			}, {  // czt via dft
				magnitude = osThis.dft(Signal.newClear(osSize)).magnitude;  // magnitude
				phase = magnitude.minimumPhase(mindb);  // minimum phase
				complex = Polar.new(magnitude.as(Signal), phase.as(Signal)).asComplex;  // minimum phase spectrum

				complex.real.idft(complex.imag).real.keep(this.size)  // idft
			})
		)
	}

	/* real even and odd */
	even {
		^(0.5 * (this + this.flip))
	}

	odd {
		^(0.5 * (this - this.flip))
	}

	/* frequency domain resampling */

	prInterpolate { arg newSize;
		var complex, real, imag;
		var scale = newSize / this.size;

		// source - real spectrum
		case
		{ this.size.isPowerOfTwo } {  // rfft
			complex = this.rfft(Signal.rfftCosTable(this.size/2 + 1));
			complex = Complex.new(  // rescale Nyquist
				complex.real.put(this.size/2, 0.5 * complex.real.at(this.size/2)),
				complex.imag.put(this.size/2, 0.5 * complex.imag.at(this.size/2))
			)
		}
		{ this.size.even } {  // zoomDft
			complex = this.rdftZoom(this.size/2 + 1, 0, this.size/2);
			complex = Complex.new(  // rescale Nyquist
				complex.real.put(this.size/2, 0.5 * complex.real.at(this.size/2)),
				complex.imag.put(this.size/2, 0.5 * complex.imag.at(this.size/2))
			)
		}
		{ this.size.odd } {  // zoomDft
			complex = this.rdftZoom((this.size + 1)/2, 0, (this.size - 1)/2)
		};

		// target
		^newSize.isPowerOfTwo.if({  // irfft
			real = Signal.newClear(newSize/2 + 1).overWrite(complex.real);
			imag = Signal.newClear(newSize/2 + 1).overWrite(complex.imag);
			scale * real.irfft(imag, Signal.rfftCosTable(newSize/2 + 1))
		}, {  // idft
			// full spectrum
			real = Signal.newClear(newSize).overWrite((complex.real.flip.drop(1) ++ complex.real).as(Signal));
			real = real.rotate(1 - complex.real.size);
			imag = Signal.newClear(newSize).overWrite((complex.imag.neg.flip.drop(1) ++ complex.imag).as(Signal));
			imag = imag.rotate(1 - complex.imag.size);
			scale * real.idft(imag).real
		})
	}

	prDecimate { arg newSize;
		var complex, real, imag;
		var scale = newSize / this.size;
		var halfsize;

		// source - real spectrum
		case
		{ this.size.isPowerOfTwo } {  // rfft
			complex = this.rfft(Signal.rfftCosTable(this.size/2 + 1));
		}
		{ this.size.even } {  // zoomDft
			complex = this.rdftZoom(this.size/2 + 1, 0, this.size/2);
		}
		{ this.size.odd } {  // zoomDft
			complex = this.rdftZoom((this.size + 1)/2, 0, (this.size - 1)/2)
		};

		// target
		^newSize.isPowerOfTwo.if({  // irfft
			halfsize = (newSize/2).asInteger;

			real = Signal.newClear(halfsize + 1).overWrite(complex.real);
			imag = Signal.newClear(halfsize + 1).overWrite(complex.imag);

			// rescale Nyquist
			real.put(halfsize, 2 * real.at(halfsize));
			imag.put(halfsize, 2 * imag.at(halfsize));

			scale * real.irfft(imag, Signal.rfftCosTable(halfsize + 1))
		}, {  // idft
			// real spectrum
			real = complex.real;
			imag = complex.imag;

			// full spectrum
			newSize.even.if({
				halfsize = (newSize/2).asInteger;

				// rescale Nyquist
				real.put(halfsize, 2 * real.at(halfsize));
				imag.put(halfsize, 2 * imag.at(halfsize));

				real = (real.copySeries(0, last: halfsize) ++ real.copySeries(halfsize - 1, last: 1)).as(Signal);
				imag = (imag.copySeries(0, last: halfsize) ++ imag.copySeries(halfsize - 1, last: 1).neg).as(Signal);
			}, {
				halfsize = ((newSize - 1)/2).asInteger;

				real = (real.copySeries(0, last: halfsize) ++ real.copySeries(halfsize, last: 1)).as(Signal);
				imag = (imag.copySeries(0, last: halfsize) ++ imag.copySeries(halfsize, last: 1).neg).as(Signal);
			});
			scale * real.idft(imag).real
		})
	}

	resize { arg newSize;
		^case
		{ this.size == newSize } { this }
		{ this.size < newSize } { this.prInterpolate(newSize) }
		{ this.size > newSize } { this.prDecimate(newSize) }
	}

	/* real fft */

	*rfftCosTable { arg rfftsize;
		^this.newClear(((rfftsize - 1) / 4 + 1)).fftCosTable
	}

	*rfftTwoCosTable { arg rfftsize;
		^this.newClear(((rfftsize - 1) / 2) + 1).fftCosTable
	}

	// method 3
	// [1] H. Sorensen; D. Jones ; M. Heideman ; C. Burrus, Real-valued fast Fourier transform algorithms, IEEE TRANSACTIONS ON ACOUSTICS, SPEECH, AND SIGNAL PROCESSING, VOL. 'ASSP-35, NO. 6, JUNE 1987
	// [2] http://www.robinscheibler.org/2013/02/13/real-fft.html
	// [3] https://dsp.stackexchange.com/questions/20115/real-valued-fft-in-matlab

	/*
	NOTE: If this.size is odd, then append zero to make even... OR, throw an error?
	      How would this compare to the results from method 1, above?
	*/

	rfft { arg cosTable;
		// argCosTable must contain 1/4 cycle of a cosine (use rfftCosTable)
		// fftsize is the next greater power of two than 1/2 the receiver's length

		var halfsize = (this.size/2).asInteger;  // time domain 1/2 size
		var rfftsize = (4 * (cosTable.size - 1)) + 1;  // number of returned coeffs
		var complexZ, complexZconj, complexE, complexO;

		// complex fft: decimate in time
		complexZ = fft(
			(this[Array.series(halfsize, 0, 2)]).as(Signal),  // even
			(this[Array.series(halfsize, 1, 2)]).as(Signal),  // odd
			cosTable
		);
		complexZconj = complexZ.conjugate;

		// reflect in frequency
		complexZconj = complexZconj.real.flip.complex(
			complexZconj.imag.flip
		);

		complexE = 0.5 * (complexZ + complexZconj);
		complexO = 0.0.complex(0.5.neg) * (complexZ - complexZconj);

		^(
			Complex.new(
				complexE.real.as(Array).wrapExtend(rfftsize).as(Signal),
				complexE.imag.as(Array).wrapExtend(rfftsize).as(Signal)
			) + (
				Complex.new(
					Signal.newClear(rfftsize).addCosine(0.5 * (rfftsize) / (rfftsize - 1)),
					Signal.newClear(rfftsize).addSine(0.5 * (rfftsize) / (rfftsize - 1), -1)
				) * Complex.new(
					complexO.real.as(Array).wrapExtend(rfftsize).as(Signal),
					complexO.imag.as(Array).wrapExtend(rfftsize).as(Signal)
				)
			)
		)
	}

	irfft { arg imag, cosTable;
		// argCosTable must contain 1/4 cycle of a cosine (use rfftCosTable)
		// fftsize is the next greater power of two than the receiver's (length - 1) / 2 (correct?)
		var size = (this.size - 1);  // fftsize == 4 * (cosTable.size - 1)
		var complex, complexConj, complexE, complexO, complexZ;  // frequency domain
		var timeZ;  // time domain - complex

		complex = this.complex(imag);
		complexConj = complex.conjugate;

		// reverse complexConj - use Array-reverse, as Signal-reverse changes receiver!
		complexConj = Complex.new(
				complexConj.real.as(Array).reverse.as(Signal),
				complexConj.imag.as(Array).reverse.as(Signal)
		);

		// discard copy of DC
		complex = Complex.new(complex.real.drop(-1), complex.imag.drop(-1));
		complexConj = Complex.new(complexConj.real.drop(-1), complexConj.imag.drop(-1));

		complexE = 0.5 * (complex + complexConj);
		complexO = 0.5 * (complex - complexConj) * Complex.new(
			Signal.newClear(size).addCosine(0.5),
			Signal.newClear(size).addSine(0.5)
		);

		complexZ = complexE + (0.0.complex(1.0) * complexO);

		// ifft - returns complex result
		timeZ = complexZ.real.ifft(complexZ.imag, cosTable);

		^[ timeZ.real.as(Array), timeZ.imag.as(Array) ].lace.as(Signal)
	}

	rfftTwo { arg real2, cosTable;
		// argCosTable must contain 1/4 cycle of a cosine (use rfftTwoCosTable)
		// fftsize is the next greater power of two than the receiver's length

		var fftsize = 4 * (cosTable.size - 1);
		var rfftsize = (fftsize / 2).asInteger + 1;  // number of returned coeffs
		var complexZ, complexZconj, complex1, complex2;

		// complex fft
		complexZ = fft(
			this,  // real1
			real2,  // real2
			cosTable
		);
		complexZconj = complexZ.conjugate;

		// reflect in frequency
		complexZconj = complexZconj.real.flip.complex(
			complexZconj.imag.flip
		);

		complex1 = 0.5 * (complexZ + complexZconj);
		complex2 = 0.0.complex(0.5.neg) * (complexZ - complexZconj);

		^Dictionary.with(*[
			\rfft1->Complex.new(
				complex1.real.keep(rfftsize),
				complex1.imag.keep(rfftsize),
			),
			\rfft2->Complex.new(
				complex2.real.keep(rfftsize),
				complex2.imag.keep(rfftsize),
			)
		])
	}

	irfftTwo { arg imag1, real2, imag2, cosTable;
		// argCosTable must contain 1/4 cycle of a cosine (use rfftTwoCosTable)
		// fftsize is the next greater power of two than the receiver's (length - 1) / 2

		var complexZ, timeZ;

		// mirror complex conjugate coefficients & add
		complexZ = Complex.new(
			// -rfftToFft
			this ++ (this.deepCopy.drop(1).drop(-1).reverse),
			imag1 ++ (imag1.deepCopy.drop(1).drop(-1).reverse.neg)
		) + Complex.new(
			(imag2 ++ (imag2.deepCopy.drop(1).drop(-1).reverse.neg)).neg,
			real2 ++ (real2.deepCopy.drop(1).drop(-1).reverse)
		);

		// ifft - returns complex result
		timeZ = complexZ.real.ifft(complexZ.imag, cosTable);

		^Dictionary.with(*[
			\irfft1->timeZ.real,
			\irfft2->timeZ.imag
		])
	}

	/* fft & rfft conversion */

	rfftToFft { arg imag;
		(this.size - 1).isPowerOfTwo.if({
			^Complex.new(
				this ++ (this.deepCopy.drop(1).drop(-1).reverse),
				imag ++ (imag.deepCopy.drop(1).drop(-1).reverse.neg)
			)
		}, {
			Error("rfftsize % is not a power of two + 1!".format(this)).throw;
		})
	}

	fftToRfft { arg imag;
		this.size.isPowerOfTwo.if({
			var rfftsize = (this.size / 2).asInteger + 1;  // number of returned coeffs
			var complexZ, complexZconj, complex1;

			// complex fft
			complexZ = this.complex(imag);
			complexZconj = complexZ.conjugate;

			// reflect in frequency
			complexZconj = complexZconj.real.flip.complex(
				complexZconj.imag.flip
			);

			// the complex spectrum of the original real signal
			complex1 = 0.5 * (complexZ + complexZconj);

			^Complex.new(
				complex1.real.keep(rfftsize),
				complex1.imag.keep(rfftsize),
			)
		}, {
			Error("fftsize % is not a power of two!".format(this)).throw;
		})
	}

	/* dft */

	dft { arg imag, method = 'czt';
		^method.switch(
			'czt', {
				this.czt(imag)
			},
			'dir', {
				var complex = this.size.collect({ arg i;
					[\real, \imag].collect({ arg msg;
						(
							Complex.new(
								Signal.newClear(this.size).addCosine(i, 1.0, 0.0),
								Signal.newClear(this.size).addSine(i, -1.0, 0.0)
							) * this.complex(imag)
						).perform(msg).sum
					});
				}).flop;

				Complex.new(complex.at(0).as(Signal), complex.at(1).as(Signal))
			}, {
				Error("DFT method % is not available!".format(method)).throw;
			}
		)
	}

	// Inverse FFT Method #3
	// [1] Duhamel P., el al, "On Computing the Inverse DFT", IEEE Trans. on Acoustics, Speech, and Signal Processing, Vol. 36, No. 2, Feb. 1988.
	// [2] https://www.dsprelated.com/showarticle/800.php
	idft { arg imag, method = 'czt';
		var complex = imag.dft(this, method) / this.size;

		^Complex.new(
			complex.imag,
			complex.real
		)
	}

	/* real dft */
	/*
	NOTE: There are ambiguities regarding even and odd sizes in the time domain.
	*/


	// default - rdft
	dftZoom { arg imag, zoomsize, k0, k1;
		var kc0, kc1;
		var cwtsize, step, start;

		// assign algorithm vars
		kc0 = (k0 != nil).if({
			k0
		}, {  // DC
			0
		});
		kc1 = (k1 != nil).if({
			k1
		}, {
			this.size.even.if({
				this.size / 2  // Nyquist
			}, {
				(this.size - 1) / 2  // highest +freq
			})
		});

		cwtsize = (zoomsize != nil).if({
			zoomsize
		}, {
			(kc1 - kc0).asInteger + 1
		});

		start = (0.0.complex(2pi * kc0 / this.size)).exp;
		step = (0.0.complex(-2pi * (kc1-kc0) / (this.size * (cwtsize-1)))).exp;

		^this.czt(imag, cwtsize, step, start)
	}

	// just a convenience - could optimize with decimation for size.even
	rdftZoom { arg zoomsize, k0, k1;
		^this.dftZoom(Signal.newClear(this.size), zoomsize, k0, k1)
	}


	/* chirp z-transform */

	/*
	Consider:
	- compute complex chirps via .sin, .cos
	- Rabiner part III optimisations
	- optimised rczt
	*/

	// complex
	// The Chirp z-Transform Algorithm, Rabiner, et al.
	// real.cwt(imag, m, w, a)
	czt { arg imag, cwtsize, step, start;
		var xn, n, m, l, w, a;
		var yn, ynChirp;
		var vn;
		var yr, vr, gr;
		var gk, xk, xkChirp;
		var cosTable;
		var reshapeComplexArray;

		// function to return Complex of two Signals, from Array of complex values
		reshapeComplexArray = { arg cArr;
			var arr = cArr.collect({ arg item; [ item.real, item.imag ] }).flop;
			Complex.new(
				arr.at(0).as(Signal),  // real
				arr.at(1).as(Signal)  // imag
			)
		};

		// assign algorithm vars
		xn = this.complex(imag);

		n = this.size;
		m = cwtsize.isNil.if({ n }, { cwtsize.asInteger });
		l = (n + m - 1).nextPowerOfTwo;  // fftsize

		w = step.isNil.if({ (0.0.complex(-2pi/m)).exp }, { step });  // matlab; Rabiner uses -2pi/n
		a = start.isNil.if({ 1.0 }, { start });

		// generate yn (xn * chirp)
		ynChirp = reshapeComplexArray.value(
			n.collect({ arg i;
				a.pow(i.neg) * w.pow(i.squared / 2)
			})
		);

		yn = ynChirp * xn;

		// generate vn (chirp)
		vn = reshapeComplexArray.value(
			m.collect({ arg i;  // 0 <= n <= M - 1
				w.pow(i.squared.neg / 2)
			}) ++ (l-n-m+1).collect({   // M <= n <= L - N
				0.0.complex(0.0)  // "arbitrary"
			}) ++ Array.series(n-1, l-n+1).collect({ arg i;  // L - N + 1 <= n <= L - 1
				w.pow((l-i).squared.neg / 2)
			})
		);

		// find gk via circular convolution
		cosTable = Signal.fftCosTable(l);

		yr = yn.real.zeroPad(l).fft(yn.imag.zeroPad(l), cosTable);
		vr = vn.real.fft(vn.imag, cosTable);
		gr = vr * yr;
		gk = gr.real.ifft(gr.imag, cosTable);

		// generate xk (gk * chirp)
		xkChirp = reshapeComplexArray.value(
			m.collect({ arg i;
				w.pow(i.squared / 2)
			})
		);

		xk = xkChirp * Complex.new(
			gk.real.keep(m),
			gk.imag.keep(m)
		);

		^xk
	}

	// just a convenience - could optimize with decimation for size.even
	rczt { arg cwtsize, step, start;
		^this.czt(Signal.newClear(this.size), cwtsize, step, start);
	}


	/* Goertzel */

	// Sysel and Rajmic: Goertzel algorithm generalized to non-integer multiples
	// of fundamental frequency. EURASIP Journal on Advances in Signal Processing 2012 2012:56.
	//
	// Vitali; The Goertzel algorithm to compute individual terms of the discrete
	// Fouier transform (DFT)

	// real "generalized" Goertzel ; k = integer or float, or array thereof
	rgoertzel { arg k, method = 'iir';
		(k.class == Array).if({
			// an arry of ks

			// function to return Complex of two Signals, from Array of complex values
			var reshapeComplexArray = { arg cArr;
				var arr = cArr.collect({ arg item; [ item.real, item.imag ] }).flop;
				Complex.new(
					arr.at(0).as(Signal),  // real
					arr.at(1).as(Signal)  // imag
				)
			};

			^reshapeComplexArray.value(
				k.collect({ arg item;
					this.rgoertzel(item, method)
				})
			)
		}, {
			// single k
			^method.switch(
				'iir', {
					var w, cosW, sinW, c;
					var w2, cosW2, sinW2;
					var z0, z1, z2;
					var realT, imagT;

					w = 2pi * k/this.size;
					cosW = w.cos;
					sinW = w.sin;
					c = 2 * cosW;

					w2 = 2pi * k;
					cosW2 = w2.cos;
					sinW2 = w2.sin;

					// initialize states
					z0 = 0.0; z1 = 0.0; z2 = 0.0;

					// recursion
					this.do({ arg item;
						z0 = item + (c * z1) - z2;
						z2 = z1;
						z1 = z0;
					});

					realT = (cosW * z1) - z2;
					imagT = (sinW * z1);

					Complex.new(
						(realT * cosW2) + (imagT * sinW2),
						(realT * sinW2.neg) + (imagT * cosW2)
					)
				},
				'czt', {
					var w, start;
					var complexSignal;

					w = 2pi * k/this.size;
					start = (0.0.complex(w)).exp;
					complexSignal = this.rczt(1, start: start);

					Complex.new(complexSignal.real.first, complexSignal.imag.first)
				},
				'dir', {
					Complex.new(
						(this * Signal.newClear(this.size).addCosine(k, 1.0)).sum,
						(this * Signal.newClear(this.size).addSine(k, -1.0)).sum
					)
				}, {
					Error("Goertzel method % is not available!".format(method)).throw;
				}
			)
		})
	}

	// complex "generalized" Goertzel ; k = integer or float, or array thereof
	goertzel { arg imag, k, method = 'iir';
		(k.class == Array).if({
			// an arry of ks

			// function to return Complex of two Signals, from Array of complex values
			var reshapeComplexArray = { arg cArr;
				var arr = cArr.collect({ arg item; [ item.real, item.imag ] }).flop;
				Complex.new(
					arr.at(0).as(Signal),  // real
					arr.at(1).as(Signal)  // imag
				)
			};

			^reshapeComplexArray.value(
				k.collect({ arg item;
					this.goertzel(imag, item, method)
				})
			)
		}, {
			// single k
			^method.switch(
				'iir', {
					var thisG, imagG;

					thisG = this.rgoertzel(k, 'iir');
					imagG = 0.0.complex(1.0) * imag.rgoertzel(k, 'iir');

					thisG + imagG
				},
				'czt', {
					var w, start;
					var complexSignal;

					w = 2pi * k/this.size;
					start = (0.0.complex(w)).exp;
					complexSignal = this.czt(imag, 1, start: start);

					Complex.new(complexSignal.real.first, complexSignal.imag.first)
				},
				'dir', {
					var complex = [\real, \imag].collect({ arg msg;
						(
							Complex.new(
								Signal.newClear(this.size).addCosine(k, 1.0),
								Signal.newClear(this.size).addSine(k, -1.0)
							) * this.complex(imag)
						).perform(msg).sum
					});

					Complex.new(complex.at(0), complex.at(1))
				}, {
					Error("Goertzel method % is not available!".format(method)).throw;
				}
			)
		})
	}

	/* cepstrum */

	rceps {
		^this.size.isPowerOfTwo.if({  // rfft
			var rfftSize = (this.size/2+1).asInteger;
			var cosTable = Signal.rfftCosTable(rfftSize);
			var imag = Signal.newClear(rfftSize);
			this.rfft(cosTable).magnitude.log.as(Signal).irfft(imag, cosTable)
		}, {  // dft
			var imag = Signal.newClear(this.size);
			this.dft(imag).magnitude.log.as(Signal).idft(imag).real
		})
	}

	irceps {
		^this.size.isPowerOfTwo.if({  // rfft
			var rfftSize = (this.size/2+1).asInteger;
			var cosTable = Signal.rfftCosTable(rfftSize);
			var complex = this.rfft(cosTable).exp;
			complex.real.irfft(complex.imag, cosTable)
		}, {  // dft
			var imag = Signal.newClear(this.size);
			var complex = this.dft(imag).exp;
			complex.real.idft(complex.imag).real
		})
	}


	/* Hilbert & Analytic */

	/*
	Include -hilbertEnvelope, -hilbertPhase ?
	*/

	*hilbert { arg size, pad = 0, sym = false;
		var hbSize = size - pad;
		var rad;
		var real, imag;

		case
		{ (hbSize.odd && sym) || (hbSize.even && sym.not) } {  // integer
			hbSize.odd.if({  // odd
				rad = pi * Array.series(hbSize, ((hbSize - 1)/2).neg);
				real = Signal.newClear(hbSize).put((hbSize - 1)/2, 1.0)  // sinc
			}, {  // even
				rad = pi * Array.series(hbSize, (hbSize/2).neg);
				real = Signal.newClear(hbSize).put(hbSize/2, 1.0)  // sinc
			});
		}
		{ (hbSize.odd && sym.not) || (hbSize.even && sym) } {  // fractional
			hbSize.odd.if({  // odd
				rad = pi * Array.series(hbSize, (hbSize/2).neg);
			}, {  // even
				rad = pi * Array.series(hbSize, ((hbSize - 1)/2).neg);
			});
			real = rad.sincPi.as(Signal);
		};
		imag = rad.collect({|i| (i == 0).if({ 0 }, { 1 - cos(i) / (i) }) }).as(Signal);

		^Complex.new(
			real,
			imag
		)
	}

	// https://www.mathworks.com/help/signal/ref/hilbert.html
	// Marple, S. L. “Computing the Discrete-Time Analytic Signal via FFT.” IEEE® Transactions on Signal Processing. Vol. 47, 1999, pp. 2600–2603.

	analytic {
		^(this.size.isPowerOfTwo).if({  // fft
			var rfft, h, real, imag;

			rfft = this.rfft(Signal.rfftCosTable(this.size/2 + 1));  // real fft
			h = ([1] ++ Array.fill(rfft.real.size - 2, { 2 }) ++ [1]).as(Signal);  // norm for +freqs

			real = (rfft.real * h).zeroPad(this.size);
			imag = (rfft.imag * h).zeroPad(this.size);

			real.ifft(imag, Signal.fftCosTable(this.size))  // complex fft
		}, {  // czt
			var czt, cztsize, step;
			var h, real, imag;
			var complex;

			step = (0.0.complex(-2pi/(this.size))).exp;  // step, same as dft

			this.size.even.if({  // even
				cztsize = (this.size/2).asInteger + 1;  // up to Nyquist
				h = ([1] ++ Array.fill(cztsize - 2, { 2 }) ++ [1]).as(Signal);  // norm for +freqs
			}, {  // odd
				cztsize = ((this.size + 1) / 2).asInteger;  // up to highest +freq
				h = ([1] ++ Array.fill(cztsize - 1, { 2 })).as(Signal);  // norm for +freqs
			});

			czt = this.czt(Signal.newClear(this.size), cztsize, step);
			real = (czt.real * h);
			imag = (czt.imag * h);

			complex = imag.czt(real, this.size, step) / this.size;  // inverse...
			Complex.new(  // ... "real" czt
				complex.imag,
				complex.real
			)
		})
	}

	*analyticFill { arg size, amplitudes, phases;
		var analytic = Signal.cosineFill(size, amplitudes, phases).analytic;
		^analytic * analytic.magnitude.maxItem.reciprocal  // normalize
	}

	/* convolution */

	convolve { arg aSignal, method = 'fft';
		var size = this.size + aSignal.size - 1;
		var complex;

		^method.switch(
			'fft', {
				var fftsize = size.nextPowerOfTwo.asInteger;
				var rfftsize = fftsize/2 + 1;
				var rfftDict = this.zeroPad(fftsize).rfftTwo(
					aSignal.zeroPad(fftsize),
					Signal.rfftTwoCosTable(rfftsize)
				);
				complex = rfftDict[\rfft1] * rfftDict[\rfft2];
				complex.real.irfft(complex.imag, Signal.rfftCosTable(rfftsize)).keep(size)
			},
			'dft', {
				complex = (this ++ Signal.newClear(aSignal.size - 1)).dft(
					Signal.newClear(size)
				) * (aSignal ++ Signal.newClear(this.size - 1)).dft(
					Signal.newClear(size)
				);
				complex.real.idft(complex.imag).real;
			},
			'dir', {
				var fac1, fac2;
				var sum = Array.zeroFill(size);
				(aSignal.size < this.size).if({
					fac1 = aSignal.as(Array);
					fac2 = this.as(Array) ++ Array.zeroFill(aSignal.size - 1);
				}, {
					fac1 = this.as(Array);
					fac2 = aSignal.as(Array) ++ Array.zeroFill(this.size - 1);
				});
				fac1.do({ arg item, i;
					sum = sum + (item * fac2.shift(i))
				});
				sum.as(Signal)
			}, {
				Error("Convolve method % is not available!".format(method)).throw;
			}
		)
	}

}
