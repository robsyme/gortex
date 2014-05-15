package cortex

import (
	. "github.com/smartystreets/goconvey/convey"
	"os"
	"testing"
)

func TestReadingBinaries(t *testing.T) {

	Convey("Subject: Opening the binary", t, func() {

		Convey("Given a non-existant cortex binary file", func() {
			filename := "../test/data/bad_binaries/non-existant.ctx"

			Convey("When opened", func() {
				_, err := Open(filename)

				Convey("It should return an error", func() {
					So(err, ShouldNotBeNil)
					So(os.IsNotExist(err), ShouldBeTrue)
				})

			})

		})

		Convey("Given a cortex binary with corrupted magic string", func() {
			filename := "../test/data/bad_binaries/bad_magic_string.ctx"

			Convey("When opened", func() {
				_, err := Open(filename)

				Convey("It should return an err", func() {
					So(err, ShouldNotBeNil)
				})
			})
		})

		Convey("Given a cortex binary with a missing byte in the header", func() {
			filename := "../test/data/bad_binaries/missing_byte.ctx"

			Convey("When opened", func() {
				_, err := Open(filename)

				Convey("It should return an err", func() {
					So(err, ShouldNotBeNil)
				})
			})
		})

		Convey("Given a correct cortex binary file", func() {
			filename := "../test/data/three_colours.ctx"

			Convey("When opened", func() {
				binary, err := Open(filename)

				Convey("It should not return an error", func() {
					So(err, ShouldBeNil)
					So(binary, ShouldNotBeNil)
				})

				Convey("We should be able to read the header information", func() {
					So(binary.Version, ShouldEqual, 6)

					So(binary.KmerSize, ShouldEqual, 21)

					So(binary.ColourCount, ShouldEqual, 3)

					So(binary.Colours[0].MeanReadLength, ShouldEqual, 70)
					So(binary.Colours[1].MeanReadLength, ShouldEqual, 70)
					So(binary.Colours[2].MeanReadLength, ShouldEqual, 70)

					So(binary.Colours[0].TotalSequenceLength, ShouldEqual, 14000)
					So(binary.Colours[1].TotalSequenceLength, ShouldEqual, 14000)
					So(binary.Colours[2].TotalSequenceLength, ShouldEqual, 14000)

					So(binary.Colours[0].Name, ShouldEqual, "org1")
					So(binary.Colours[1].Name, ShouldEqual, "org2")
					So(binary.Colours[2].Name, ShouldEqual, "org3")
				})

				Convey("Given the first kmer", func() {
					kmer, err := binary.NextKmer()

					Convey("There should be no error", func() {
						So(err, ShouldBeNil)
						So(kmer, ShouldNotBeNil)
					})

					Convey("We should be able to retrive the nucleotide sequence", func() {
						nucs := binary.KmerNucleotides(kmer)
						So(nucs, ShouldEqual, "CGCACTTTGCATCGCTGGCCA")
					})
				})

			})

		})

	})

}
