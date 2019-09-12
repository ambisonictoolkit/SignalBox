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
// 	Extension: Complex
//
//---------------------------------------------------------------------


+ Complex {

	// FreqSpectrum class
	*newFreqSpectrum { arg spectrum;
		^spectrum.asComplex
	}

	asFreqSpectrum { ^FreqSpectrum.newComplex(this) }

	// could add to MathLib Quark...
	rotate { arg angle;
		^this.asPolar.rotate(angle).asComplex
	}

}
