package colouredKmer

import (
	. "github.com/smartystreets/goconvey/convey"
	"testing"
)

func TestKmerBasics(t *testing.T) {
	Convey("Given an empty kmer", t, func() {
		kmer := NewKmer()
		Convey("The kmer should not be nil", func() {
			So(kmer, ShouldNotBeNil)
		})
	})

	Convey("Given a Kmer encoding 'ACGT'", t, func() {
		//kmer := NewKmerFromInt(228)
		Convey("The kmer should be able to return the sequence string", func() {
		})
	})

	Convey("Given a short string of nucleotides", t, func() {
		seq := "acgt"
		Convey("The package should be able to calculate the big.Int that encodes it", func() {
			a, err := stringToBigInt(seq)
			So(err, ShouldBeNil)
			So(a.Int64(), ShouldEqual, 27)
		})

		Convey("You can create a enw Kmer from the string", func() {
			a, err := NewKmerFromString(seq)
			So(err, ShouldBeNil)
			So(a.Int64(), ShouldEqual, 27)

			Convey("Which should recalculate the original sequence", func() {
				So(a.Letters(len(seq)).String(), ShouldEqual, seq)
			})
		})
	})

	Convey("Given a string with non-[acgt] nucleotides", t, func() {
		seq := "acgtx"
		Convey("Creating a Kmer from the string should create an error", func() {
			_, err := NewKmerFromString(seq)
			So(err, ShouldNotBeNil)
		})
	})
}
