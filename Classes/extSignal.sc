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

	/* rotation */

	// NOTE: match Hilbert phase rotation
	// -rotatePeriod ??
	rotateWave { arg phase;
		var precision = 0.1;  // corresponds to integer
		var n = phase / 2pi.neg * this.size;
		var intN = n.asInteger;

		^n.equalWithPrecision(intN, precision).if({  // direct
			this.deepCopy.rotate(n)
		}, {
			(this.size.isPowerOfTwo).if({  // fft
				var cosTable = Signal.rfftCosTable(this.size / 2 + 1);
				var complex = this.rfft(cosTable);
				var rotatedComplex =  complex.magnitude.polarComplexSignal(
					complex.phase + phase
				);
				rotatedComplex.real.irfft(rotatedComplex.imag, cosTable)
			}, {  // dft
				var complex = this.dft(Signal.newClear(this.size));
				var rotatedComplex;
				var halfsize = (this.size/2).floor;

				rotatedComplex = this.size.even.if({
					complex.magnitude.polarComplexSignal(
						complex.phase + (
							Array.fill(halfsize, { phase }) ++ Array.fill(halfsize, { phase.neg })
						)
					)
				}, {
					complex.magnitude.polarComplexSignal(
						complex.phase + (
							Array.fill(halfsize + 1, { phase }) ++ Array.fill(halfsize, { phase.neg })
						)
					)
				});
				rotatedComplex.real.idft(rotatedComplex.imag).real
			})
		})
	}

	rotatePhase { arg phase;
		var complex;

		^(this.size.isPowerOfTwo).if({  // rfft
			var cosTable = Signal.rfftCosTable(this.size/2 + 1);

			complex = this.rfft(cosTable);  // real fft
			complex = Complex.new(phase.cos, phase.sin) * complex;  // rotate

			complex.real.irfft(complex.imag, cosTable)  // irfft
		}, {  // czt via analytic
			complex = this.analytic;

			// equivalent to: (Complex.new(phase.cos, phase.sin) * complex).real
			(phase.cos * complex.real) + (phase.sin.neg * complex.imag)
		})
	}

	/* real even and odd */
	even {
		^(0.5 * (this + this.flip))
	}

	odd {
		^(0.5 * (this - this.flip))
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


	// /* real dft */
	//
	// rdft {
	// 	var complex;
	// 	var rdftSize;
	//
	// 	rdftSize = this.size.even.if({
	// 		(this.size / 2).asInteger + 1
	// 	}, {
	// 		((this.size + 1)/ 2).asInteger
	// 	});
	//
	// 	complex = rdftSize.collect({ arg i;
	// 		[\real, \imag].collect({ arg msg;
	// 			(
	// 				Complex.new(
	// 					Signal.newClear(this.size).addSine(i, 1.0, 0.5pi) * this,
	// 					Signal.newClear(this.size).addSine(i, 1.0, pi) * this
	// 				)
	// 			).perform(msg).sum
	// 		});
	// 	}).flop;
	//
	// 	^Complex.new(complex.at(0).as(Signal), complex.at(1).as(Signal))
	// }
	//
	//
	// /*
	// Revise to use design pattern(s) found in -idft & -rdft
	// */
	// // match irfft return
	// // doesn't account for this.size == odd
	// // irdft { arg imag;
	// // 	var size;
	// // 	var magnitudePhase;
	// // 	var sineFill2List;
	// //
	// // 	size = 2 * (this.size - 1);
	// //
	// // 	magnitudePhase = [\magnitude, \phase].collect({ arg msg;
	// // 		Complex.new(this, imag).perform(msg)
	// // 	});
	// //
	// // 	sineFill2List = this.size.collect({ arg i;
	// // 		[ i, magnitudePhase.at(0).at(i), 0.5pi + magnitudePhase.at(1).at(i) ]
	// // 	});
	// //
	// // 	^Signal.newClear(size).sineFill2(sineFill2List).scale(2 * size.reciprocal)
	// // }
	// irdft { arg imag;
	// 	var size;
	// 	var magnitudePhase;
	// 	var sineFill2List;
	//
	// 	size = this.size.even.if({
	// 		(2 * this.size) - 1
	// 	}, {
	// 		2 * (this.size - 1)
	// 	});
	//
	// 	magnitudePhase = [\magnitude, \phase].collect({ arg msg;
	// 		Complex.new(this, imag).perform(msg)
	// 	});
	//
	// 	/*
	// 	Take account of scaling at DC and Nyquist.
	// 	Different for size.even and size.odd
	//
	// 	this.size == power or two - 1;  // irdft works!
	// 	this.size == power or two;  // irdft works!
	// 	this.size == power or two + 1;  // irdft fails!
	//
	// 	*/
	// 	magnitudePhase.put(
	// 		0,
	// 		magnitudePhase.at(0).performInPlace(
	// 			'*',
	// 			1,
	// 			this.size.even.if({
	// 				this.size - 1
	// 			}, {
	// 				this.size - 2
	// 			}),
	// 			2.0
	// 		)
	// 	);
	//
	// 	sineFill2List = this.size.collect({ arg i;
	// 		[ i, magnitudePhase.at(0).at(i), 0.5pi + magnitudePhase.at(1).at(i) ]  // real
	// 	});
	//
	// 	^Signal.newClear(size).sineFill2(sineFill2List).scale(size.reciprocal)
	// }

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
		m = cwtsize.isNil.if({ n }, { cwtsize });
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
					this.goertzel(item, method)
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

	/* Hilbert & Analytic */

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
				cztsize = ((this.size + 1) /2).asInteger;  // up to highest +freq
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

}
