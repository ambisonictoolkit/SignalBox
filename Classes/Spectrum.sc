Spectrum : Number {
	var <>magnitude, <>phase;

	*new { arg magnitude, phase;
		^(phase == nil).if({
			super.newCopyArgs(magnitude, Array.zeroFill(magnitude.size))
		}, {
			super.newCopyArgs(magnitude, phase)
		})
	}

	*newComplex { arg complex;
		var polar = complex.asPolar;
		^Spectrum.new(polar.magnitude, polar.phase)
	}

	*newLinear { arg magnitude, sym = false;
		^Spectrum.new(magnitude, Array.zeroFill(magnitude.size)).linearPhase(sym)
	}

	*newMinimum { arg magnitude, mindb = -120.0;
		^Spectrum.new(magnitude, Array.zeroFill(magnitude.size)).minimumPhase(mindb)
	}

	rho { ^magnitude }

	angle { ^phase }
	theta { ^phase }

	// Signal classs / fft
	real { ^(magnitude * cos(phase)).as(Signal) }
	imag { ^(magnitude * sin(phase)).as(Signal) }

	asSpectrum { ^this }
	asPolar { ^Polar.new(this.magnitude, this.phase) }
	asComplex { ^Complex.new(this.real, this.imag) }

	size { ^this.magnitude.size }

	// phase
	linearPhase { arg sym = false;
		var start, step;
		var phase;

		sym.if({
			step = pi.neg * (this.size-1) / this.size;
		}, {
			step = pi.neg;
		});

		phase = this.size.even.if({
			start = step.neg * this.size / 2;  // start with negative freqs
			Array.series(this.size, start, step).rotate((this.size / 2).asInteger)  // then rotate
		}, {
			start = step.neg * (this.size-1) / 2;  // start with negative freqs
			Array.series(this.size, start, step).rotate(((this.size+1) / 2).asInteger)  // then rotate
		});

		^Spectrum.new(this.magnitude, phase)
	}

	minimumPhase { arg mindb = -120.0;
		var logMag = this.magnitude.max(this.magnitude.maxItem * mindb.dbamp).log;
		var phase = logMag.as(Signal).analytic.imag.neg.as(Array);  // -1 * Hilbert
		^Spectrum.new(this.magnitude, phase)
	}

	// rotatePhase {}


	/*
	Implement ....
	*/

	// scale { arg scale;
	// 	^Polar.new(rho * scale, theta)
	// }
	// rotate { arg angle; // in radians
	// 	^Polar.new(rho, theta + angle)
	// }
	//
	// // do math as Complex
	// + { arg aNumber;  ^this.asComplex + aNumber  }
	// - { arg aNumber;  ^this.asComplex - aNumber  }
	// * { arg aNumber;  ^this.asComplex * aNumber  }
	// / { arg aNumber;  ^this.asComplex / aNumber  }
	//
	// == { arg aPolar;
	// 	^aPolar respondsTo: #[\rho, \theta] and: {
	// 		rho == aPolar.rho and: { theta == aPolar.theta }
	// 	}
	// }

	hash {
		^magnitude.hash bitXor: phase.hash
	}

	// neg { ^Polar.new(rho, theta + pi) }
	//
	// performBinaryOpOnSomething { |aSelector, thing, adverb|
	// 	^thing.asComplex.perform(aSelector, this, adverb)
	// }
	//
	// performBinaryOpOnUGen { arg aSelector, aUGen;
	// 	^Complex.new(
	// 		BinaryOpUGen.new(aSelector, aUGen, this.real),
	// 		BinaryOpUGen.new(aSelector, aUGen, this.imag)
	// 	);
	// }


	printOn { arg stream;
		stream << "Spectrum( " << magnitude << ", " << phase << " )";
	}

	storeArgs { ^[magnitude, phase] }
}
