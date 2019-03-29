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
// 	Extension: Integer
//
//---------------------------------------------------------------------


+ Integer {

	/* fft helper */

	// the receiver is fft.real.size or fft.imag.size
	fftFreqs { arg sampleRate;
		var fftHalfSize, nyquist, step;

		this.isPowerOfTwo.if({
			fftHalfSize = this / 2;
			nyquist = (sampleRate == nil).if({
				1  // use normalized freq
			}, {
				sampleRate / 2  // or actual sampleRate
			});
			step = nyquist / fftHalfSize;

			^Array.series(fftHalfSize, 0, step) ++ Array.series(fftHalfSize, nyquist.neg, step)
		}, {
			Error("fftsize % is not a power of two!".format(this)).throw;
		})
	}


	/* rfft helper */

	// the receiver is rfft.real.size or rfft.imag.size
	rfftFreqs { arg sampleRate;
		var nyquist, step;

		(this - 1).isPowerOfTwo.if({
			nyquist = (sampleRate == nil).if({
				1  // use normalized freq
			}, {
				sampleRate / 2  // or actual sampleRate
			});
			step = nyquist / (this - 1);

			^Array.series(this, step: step)
		}, {
			Error("rfftsize % is not a power of two + 1!".format(this)).throw;
		})
	}


	/* dft helper */

	// the receiver is signalSize
	dftFreqs { arg sampleRate;
		var nyquist, step;

		nyquist = (sampleRate == nil).if({
			1  // use normalized freq
		}, {
			sampleRate / 2  // or actual sampleRate
		});
		step = 2 / this;  // normalized frequency

		^(nyquist * ((Array.series(this, step: step) + 1).mod(2) - 1))
	}

}
