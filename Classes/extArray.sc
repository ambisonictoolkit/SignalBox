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
// 	Extension: Array
//
//---------------------------------------------------------------------


+ Array {

	/* zeros */

	*zeroFill { arg size;
		^Array.fill(size, {0.0} )
	}

	/* magnitude */
	*logShelf { arg size, freq0, freq1, gainDC, gainNy, sampleRate;
		var f0, f1;
		var delta, sigma;
		var kdc, kny;
		var freqs;

		(freq0.abs < freq1.abs).if({
			f0 = freq0.abs;
			f1 = freq1.abs;
		}, {
			f0 = freq1.abs;
			f1 = freq0.abs;
		});

		delta = (f1 / f0).log2.reciprocal;
		sigma = (f0 * f1).log2;
		kdc = gainDC.dbamp;
		kny = gainNy.dbamp;

		freqs = size.isPowerOfTwo.if({
			size.fftFreqs(sampleRate)
		}, {
			size.dftFreqs(sampleRate)
		});

		^freqs.collect({ arg freq;
			var beta, sinBeta2, cosBeta2;
			var freqAbs = freq.abs;

			case
			{ freqAbs <= f0 } {
				kdc
			}
			{ (freqAbs > f0) && (freqAbs < f1) } {
				beta = (pi/4) * (1 + (delta * (sigma - (2 * freqAbs.log2))));  // direct beta
				sinBeta2 = beta.sin.squared;
				cosBeta2 = 1 - sinBeta2;

				kdc.pow(sinBeta2) * kny.pow(cosBeta2)  // return as scale
			}
			{ freqAbs >= f1 } {
				kny
			}
		})
	}


	/* phase helpers */

	// input is magnitude of +/-frequencies
	linearPhase { arg sym = false;
		var start, step;
		var phase;

		sym.if({
			step = pi.neg * (this.size-1) / this.size;
		}, {
			step = pi.neg;
		});

		this.size.even.if({
			start = step.neg * this.size / 2;  // start with negative freqs
			phase = Array.series(this.size, start, step);
			phase = phase.rotate((this.size / 2).asInteger)  // rotate
		}, {
			start = step.neg * (this.size-1) / 2;  // start with negative freqs
			phase = Array.series(this.size, start, step);
			phase = phase.rotate(((this.size+1) / 2).asInteger)  // rotate
		});
		^phase
	}

	// input is magnitude of +/-frequencies
	minimumPhase { arg mindb = -120.0;
		var logMag = this.max(this.maxItem * mindb.dbamp).log;
		^logMag.as(Signal).analytic.imag.neg.as(Array)  // -1 * Hilbert
	}

}
