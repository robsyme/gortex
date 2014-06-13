package colouredKmer

import (
	"code.google.com/p/biogo/alphabet"
	"fmt"
	"math/big"
)

type Kmer struct {
	*big.Int
	Coverages []uint32
	Edges     []uint8
}

func NewKmer() *Kmer {
	kmer := Kmer{Int: big.NewInt(0)}
	return &kmer
}

func NewKmerFromString(seq string) (*Kmer, error) {
	kmer := new(Kmer)
	bits, err := stringToBigInt(seq)
	if err != nil {
		return kmer, err
	}
	kmer.Int = bits
	return kmer, nil
}

func (kmer *Kmer) Letters(k int) *alphabet.Letters {
	out := make(alphabet.Letters, k)
	bits := kmer.Bits()

	for i := range out {
		wordIndex := i / 32
		shiftDist := uint(i-wordIndex*32) * 2
		bitPair := int(bits[wordIndex] >> shiftDist & 3)
		out[len(out)-i-1] = alphabet.DNA.Letter(bitPair)
	}
	return &out
}

// Takes a slice of alphabet.Letter and returns the big.Int binary
// representation of that sequence.
func lettersToBigInt(seq alphabet.Letters) (*big.Int, error) {
	out := big.NewInt(0)
	words := make([]big.Word, len(seq)/33+1)
	for i := range seq {
		index := alphabet.DNA.IndexOf(seq[len(seq)-i-1])
		if index < 0 {
			return out, fmt.Errorf("Sequence is not a valid DNA sequence at position %d\n", i+1)
		} else {
			wordIndex := i / 32
			shiftDist := uint(i-wordIndex*32) * 2
			words[wordIndex] |= big.Word(index << shiftDist)
		}
	}
	return out.SetBits(words), nil
}

func stringToBigInt(seq string) (*big.Int, error) {
	return lettersToBigInt(alphabet.Letters(seq))
}
