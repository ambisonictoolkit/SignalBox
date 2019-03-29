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
// 	Extension: SimpleNumber
//
//---------------------------------------------------------------------


+ SimpleNumber {

	/* Goertzel, DFT helper */

	// the receiver is frequency
	freqK { arg size, sampleRate;
		var samplePeriod;

		samplePeriod = (sampleRate != nil).if({
			sampleRate.reciprocal
		}, {
			0.5
		});

		^((samplePeriod * this).mod(1) * size)
	}

	// the receiver is k
	kFreq { arg size, sampleRate;
		^(sampleRate != nil).if({
			((this / size + 0.5).mod(1) - 0.5) * sampleRate
		}, {
			((this / size + 0.5).mod(1) - 0.5) * 2
		})
	}
}
