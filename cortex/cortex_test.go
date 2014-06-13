package cortex

import (
	. "github.com/smartystreets/goconvey/convey"
	"os"
	"testing"
)

func TestReadingBinaries(t *testing.T) {
	Convey("Subject: Opening cortex_var binaries", t, func() {

		Convey("When opening a non-existant cortex binary file", func() {
			filename := "../test/data/bad_binaries/non-existant.ctx"
			_, err := Open(filename)
			Convey("It should return an error", func() {
				So(err, ShouldNotBeNil)
				So(os.IsNotExist(err), ShouldBeTrue)
			})
		})

		Convey("When opening a cortex binary with corrupted magic string", func() {
			filename := "../test/data/bad_binaries/bad_magic_string.ctx"
			_, err := Open(filename)
			Convey("It should return an err", func() {
				So(err, ShouldNotBeNil)
			})
		})

		Convey("When opening a cortex binary with a missing byte in the header", func() {
			filename := "../test/data/bad_binaries/missing_byte.ctx"
			_, err := Open(filename)
			Convey("It should return an err", func() {
				So(err, ShouldNotBeNil)
			})
		})

		Convey("When opening a correct cortex binary file", func() {
			filename := "../test/data/three_colours.ctx"
			cortexBinary, err := Open(filename)

			Convey("It should not return an error", func() {
				So(err, ShouldBeNil)
				So(cortexBinary, ShouldNotBeNil)
			})

			Convey("We should be able to read the header information", func() {
				So(cortexBinary.Version, ShouldEqual, 6)
				So(cortexBinary.KmerSize, ShouldEqual, 21)
				So(cortexBinary.ColourCount, ShouldEqual, 3)
				So(cortexBinary.Colours[0].MeanReadLength, ShouldEqual, 70)
				So(cortexBinary.Colours[1].MeanReadLength, ShouldEqual, 70)
				So(cortexBinary.Colours[2].MeanReadLength, ShouldEqual, 70)
				So(cortexBinary.Colours[0].TotalSequenceLength, ShouldEqual, 14000)
				So(cortexBinary.Colours[1].TotalSequenceLength, ShouldEqual, 14000)
				So(cortexBinary.Colours[2].TotalSequenceLength, ShouldEqual, 14000)
				So(cortexBinary.Colours[0].Name, ShouldEqual, "org1")
				So(cortexBinary.Colours[1].Name, ShouldEqual, "org2")
				So(cortexBinary.Colours[2].Name, ShouldEqual, "org3")
			})

			Convey("We should be able to read all the kmers", func() {
				kmers := cortexBinary.Kmers()
				kmerCount := 0
				for _ = range kmers {
					kmerCount += 1
				}
				cortexBinary.Close()
				So(kmerCount, ShouldEqual, cortexBinary.KmerCount)
			})

		})

	})
}
