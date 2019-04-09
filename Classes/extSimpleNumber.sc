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

	// dB utilities - add also to SequenceableCollection
	maxdb { arg aNumber, adverb;
		var amp = aNumber.dbamp;
		^this.isPositive.if({
			this.max(amp, adverb)
		}, {
			this.neg.max(amp, adverb).neg
		})
	}

	mindb { arg aNumber, adverb;
		var amp = aNumber.dbamp;
		^this.isPositive.if({
			this.min(amp, adverb)
		}, {
			this.neg.min(amp, adverb).neg
		})
	}

	clipdb { arg lo, hi;
		var loAmp = lo.dbamp;
		var hiAmp = hi.dbamp;
		^this.isPositive.if({
			this.max(loAmp).min(hiAmp)
		}, {
			this.neg.max(loAmp).min(hiAmp).neg
		})
	}

	clipdb2 { arg aNumber, adverb;
		var loAmp = aNumber.neg.dbamp;
		var hiAmp = aNumber.dbamp;
		^this.isPositive.if({
			this.max(loAmp).min(hiAmp)
		}, {
			this.neg.max(loAmp).min(hiAmp).neg
		})
	}

	threshdb { arg aNumber, adverb;
		var amp = aNumber.dbamp;
		^this.isPositive.if({
			this.thresh(amp, adverb)
		}, {
			this.neg.thresh(amp, adverb).neg
		})
	}

}
