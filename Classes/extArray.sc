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

	/* complex helper */

	// return Complex from magnitude, phase - useful for fft
	polarComplexSignal { arg phase;
		^(this.as(Signal) * phase.cos.as(Signal).complex(phase.sin.as(Signal)))
	}
}
