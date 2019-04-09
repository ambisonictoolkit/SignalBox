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
// 	Extension: SequenceableCollection
//
//---------------------------------------------------------------------


+ SequenceableCollection {

	// max { arg aNumber=0, adverb; ^this.performBinaryOp('max', aNumber, adverb) }

	maxdb { arg aNumber, adverb; ^this.performBinaryOp('maxdb', aNumber, adverb) }
	mindb { arg aNumber, adverb; ^this.performBinaryOp('mindb', aNumber, adverb) }
	clipdb2 { arg aNumber, adverb; ^this.performBinaryOp('clipdb2', aNumber, adverb) }
	threshdb { arg aNumber, adverb; ^this.performBinaryOp('threshdb', aNumber, adverb) }

	clipdb { arg ... args; ^this.multiChannelPerform('clipdb', *args) }
}
