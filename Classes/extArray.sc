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

	/* Goertzel, DFT helper */

	// the receiver is an array of frequencies
	freqKs { arg size, sampleRate;
		^this.collect({ arg item;
			item.freqK(size, sampleRate)
		})
	}

	// the receiver is an array of ks
	kFreqs { arg size, sampleRate;
		^this.collect({ arg item;
			item.kFreq(size, sampleRate)
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
	minimumPhase { var mindb = -120.0;
		var logMag = this.max(this.maxItem * mindb.dbamp).log;
		^logMag.as(Signal).analytic.imag.neg.as(Array)  // -1 * Hilbert
	}

}
