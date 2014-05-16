package cortex

import (
	. "github.com/smartystreets/goconvey/convey"
	"os"
	"testing"
)

func TestReadingBinaries(t *testing.T) {

	Convey("Subject: Operating on kmers", t, func() {
		Convey("Given two kmers that represent the same node", func() {
			kmer1 := Kmer{
				BinaryKmer: []uint64{857198321847},
				Coverages:  []uint32{10, 10, 10},
			}
			kmer2 := Kmer{
				BinaryKmer: []uint64{857198321847},
				Coverages:  []uint32{12, 8, 20},
			}

			Convey("They should pass a test for equality", func() {
				So(kmer1.Equals(kmer2), ShouldBeTrue)
			})
		})

		Convey("We should be able to generate a kmer from utf-8 sequence", func() {
			seqs := "TATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
			kmer := NewKmerFromSequence(seqs)
			So(kmer, ShouldNotBeNil)
			So(kmer.BinaryKmer[0], ShouldEqual, uint64(0x3fffffffffffffff))
			So(kmer.BinaryKmer[1], ShouldEqual, uint64(0x3))
		})

		Convey("Given a kmer with k < 33", func() {
			seq := "TATTTTTTTGTTTTTTTTTTTTTTTTTGT"
			k := uint32(len(seq))
			kmer := NewKmerFromSequence(seq)

			Convey("There should be no error", func() {
				So(kmer, ShouldNotBeNil)
			})

			Convey("We should be able to retrive the nucleotide sequence", func() {
				nucs := kmer.nucleotides(k)
				So(string(nucs), ShouldEqual, seq)
				rev_nucs := kmer.reverse_nucleotides(k)
				So(string(rev_nucs), ShouldEqual, reverse(seq))
			})

			Convey("We should be able to calculate the bitsting of the incoming nodes", func() {
				kmer.ColouredEdges = []edges{edges(4), edges(5), edges(6)}
			})
		})

		Convey("Given a kmer wih large values for k", func() {
			seq := "TCCGTTTTTTTTATGCATGCATGCTTGATCGTATGCGTTTTTTTTGACGTATGCATGCTGACTGATCGATGCTGACTT"
			k := uint32(len(seq))

			kmer := NewKmerFromSequence(seq)
			Convey("We should be able to retrive the nucleotide sequence", func() {
				nucs := kmer.nucleotides(k)
				So(string(nucs), ShouldEqual, seq)
				rev_nucs := kmer.reverse_nucleotides(k)
				So(string(rev_nucs), ShouldEqual, reverse(seq))
			})

		})

		Convey("Given a kmer with incoming and outgoing nodes", func() {
			seq := "TATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTGT"
			k := uint32(len(seq))
			kmer := NewKmerFromSequence(seq)
			kmer.ColouredEdges = []edges{
				edges(0xA4), // tc -> seq -> c
				edges(0xC2)} // tg -> seq -> g

			Convey("It should be able to generate the kmers for prev and next nodes", func() {
				leftKmers := kmer.LeftKmers(k)
				So(len(leftKmers), ShouldEqual, 3)
				leftKmerC := NewKmerFromSequence("C" + seq[:len(seq)-1])
				leftKmerG := NewKmerFromSequence("G" + seq[:len(seq)-1])
				leftKmerT := NewKmerFromSequence("T" + seq[:len(seq)-1])
				So(leftKmers[0].Equals(leftKmerC), ShouldBeTrue)
				So(leftKmers[1].Equals(leftKmerG), ShouldBeTrue)
				So(leftKmers[2].Equals(leftKmerT), ShouldBeTrue)

				rightKmers := kmer.RightKmers(k)
				So(len(rightKmers), ShouldEqual, 2)
				So(string(rightKmers[0].nucleotides(k)[:]), ShouldEqual, seq[1:]+"G")
				So(string(rightKmers[1].nucleotides(k)[:]), ShouldEqual, seq[1:]+"C")
				rightKmerG := NewKmerFromSequence(seq[1:] + "G")
				rightKmerC := NewKmerFromSequence(seq[1:] + "C")
				So(rightKmers[0].Equals(rightKmerG), ShouldBeTrue)
				So(rightKmers[1].Equals(rightKmerC), ShouldBeTrue)
			})
		})

	})

	Convey("Subject: Opening cortex_var binaries", t, func() {

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

			})

		})

	})

}

func reverse(input string) string {
	// Get Unicode code points.
	n := 0
	rune := make([]rune, len(input))

	for _, r := range input {
		rune[n] = r
		n++
	}

	rune = rune[0:n]
	// Reverse
	for i := 0; i < n/2; i++ {
		rune[i], rune[n-1-i] = rune[n-1-i], rune[i]
	}

	return string(rune)
}
