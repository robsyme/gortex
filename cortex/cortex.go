package cortex

import (
	"bufio"
	"bytes"
	"encoding/binary"
	"errors"
	"io"
	"os"
)

type Colour struct {
	Name                      string
	MeanReadLength            uint32
	TotalSequenceLength       uint64
	ErrorRate                 [16]byte
	TopClippingPerformed      bool
	LowCovSupernodesRemoved   bool
	LowCovKmersRevoved        bool
	WasCleanedAgainstGraph    bool
	LowCovSupernodesThreshold uint32
	LowCovKmersThreshold      uint32
	NameOfCleaningGraph       string
}

type CortexVarBinary struct {
	reader              io.Reader
	Version             uint32
	KmerSize            uint32
	wordsPerKmer        uint32
	ColourCount         uint32
	Colours             []Colour
	kmerStartFileOffset int64
}

func (bin *CortexVarBinary) Read(data interface{}) error {
	err := binary.Read(bin.reader, binary.LittleEndian, data)
	return err
}

func (bin *CortexVarBinary) NextKmer() (Kmer, error) {
	kmer := NewKmer(bin.wordsPerKmer, bin.ColourCount)

	for i := uint32(0); i < bin.wordsPerKmer; i++ {
		err := bin.Read(&kmer.BinaryKmer[i])
		if err != nil {
			return kmer, err
		}
	}

	for i := range bin.Colours {
		err := bin.Read(&kmer.Coverages[i])
		if err != nil {
			return kmer, err
		}
	}
	for i := range bin.Colours {
		err := bin.Read(&kmer.ColouredEdges[i])
		if err != nil {
			return kmer, err
		}
	}

	return kmer, nil
}

func (bin *CortexVarBinary) readHeader() error {
	if !bin.hasMagicString() {
		return errors.New("Cortex file does not have correct format.")
	}

	bin.Read(&bin.Version)
	bin.Read(&bin.KmerSize)
	bin.Read(&bin.wordsPerKmer)
	bin.Read(&bin.ColourCount)
	bin.Colours = make([]Colour, bin.ColourCount)

	for i := range bin.Colours {
		bin.Read(&bin.Colours[i].MeanReadLength)
	}

	for i := range bin.Colours {
		bin.Read(&bin.Colours[i].TotalSequenceLength)
	}

	for i := range bin.Colours {
		// How many bytes should we read before the next name?
		var nameLength uint32
		bin.Read(&nameLength)
		// If we see a really long name, it probably means the
		// file is corrupt.
		if nameLength > 10000 {
			return errors.New("Cortex file does not have the correct format.")
		}

		// Create a buffer to hold the name bytes as a whole
		var b bytes.Buffer
		//   and a char to read each byte into
		var c byte
		for j := uint32(0); j < nameLength; j++ {
			bin.Read(&c)
			// We're not interested in adding zero bytes
			if c > 0 {
				b.WriteByte(c)
			}
		}
		// Store the resulting string in our bin object.
		bin.Colours[i].Name = b.String()
	}

	// For the moment, the error values are just stored. Perhaps
	// they should be stored as a custom type with methods?
	for i := range bin.Colours {
		bin.Read(&bin.Colours[i].ErrorRate)
	}

	var tmpByte byte
	var nameLength uint32
	for i := range bin.Colours {
		colour := bin.Colours[i]
		bin.Read(&tmpByte)
		colour.TopClippingPerformed = tmpByte == 0x1
		bin.Read(&tmpByte)
		colour.LowCovSupernodesRemoved = tmpByte == 0x1
		bin.Read(&tmpByte)
		colour.LowCovKmersRevoved = tmpByte == 0x1
		bin.Read(&tmpByte)
		colour.WasCleanedAgainstGraph = tmpByte == 0x1
		bin.Read(&colour.LowCovSupernodesThreshold)
		bin.Read(&colour.LowCovKmersThreshold)

		bin.Read(&nameLength)
		// If we see a really long name, it probably means the
		// file is corrupt.
		if nameLength > 10000 {
			return errors.New("Cortex file does not have the correct format.")
		}
		var b bytes.Buffer
		var c byte
		for j := uint32(0); j < nameLength; j++ {
			bin.Read(&c)
			if c > 0 {
				b.WriteByte(c)
			}
		}
		colour.NameOfCleaningGraph = b.String()
		bin.Colours[i] = colour
	}

	if !bin.hasMagicString() {
		return errors.New("Cortex file does not have correct format.")
	}

	return nil
}

func (bin *CortexVarBinary) hasMagicString() bool {
	var magicString [6]byte
	bin.Read(&magicString)

	expected := [6]byte{'C', 'O', 'R', 'T', 'E', 'X'}
	if magicString != expected {
		return false
	}

	return true
}

func (bin *CortexVarBinary) KmerNucleotides(kmer Kmer) string {
	nucs := kmer.nucleotides(bin.KmerSize)
	return string(nucs[:])
}

func (bin *CortexVarBinary) KmerNucleotidesReverse(kmer Kmer) string {
	nucs := kmer.reverse_nucleotides(bin.KmerSize)
	return string(nucs[:])
}

type edges byte

type Kmer struct {
	BinaryKmer    []uint64
	Coverages     []uint32
	ColouredEdges []edges
}

func NewKmer(wordsPerKmer, colourCount uint32) (kmer Kmer) {
	kmer.BinaryKmer = make([]uint64, wordsPerKmer)
	kmer.Coverages = make([]uint32, colourCount)
	kmer.ColouredEdges = make([]edges, colourCount)
	return kmer
}

func NewKmerFromSequence(seq string) Kmer {
	var kmer Kmer
	k := len(seq)
	base := make([]uint64, (len(seq)-1)/32+1)

	for i, c := range seq {
		wordIndex := (k - i - 1) / 32
		switch c {
		case 'A':
			base[wordIndex] = base[wordIndex]<<2 | 0
		case 'C':
			base[wordIndex] = base[wordIndex]<<2 | 1
		case 'G':
			base[wordIndex] = base[wordIndex]<<2 | 2
		case 'T':
			base[wordIndex] = base[wordIndex]<<2 | 3
		}
	}

	kmer.BinaryKmer = base
	return kmer
}

func (kmer *Kmer) Equals(other Kmer) bool {
	for i, v := range kmer.BinaryKmer {
		if v != other.BinaryKmer[i] {
			return false
		}
	}
	return true
}

func (kmer *Kmer) nucleotides(k uint32) []byte {
	//fmt.Printf("GENERATING NEW NUCS, k=%d\n", k)
	nucs := make([]byte, k)
	for i := k; i > 0; i-- {
		j := (i - 1) % 32
		wordIndex := (i - 1) / 32
		mask := uint64(3 << (2 * j))
		// fmt.Printf("i=%d, j=%d k-1=%d\n", i, j, k-i)
		// fmt.Printf("mask:   %064b\n", mask)
		// fmt.Printf("binary: %064b\n\n", kmer.BinaryKmer[wordIndex])
		switch mask & kmer.BinaryKmer[wordIndex] >> (j * 2) {
		case 0:
			nucs[k-i] = 'A'
		case 1:
			nucs[k-i] = 'C'
		case 2:
			nucs[k-i] = 'G'
		case 3:
			nucs[k-i] = 'T'
		}
	}
	return nucs
}

func (kmer *Kmer) reverse_nucleotides(k uint32) []byte {
	nucs := make([]byte, k)
	for i := k; i > 0; i-- {
		j := (i - 1) % 32
		wordIndex := (i - 1) / 32
		mask := uint64(3 << (2 * j))
		switch mask & kmer.BinaryKmer[wordIndex] >> (j * 2) {
		case 0:
			nucs[i-1] = 'A'
		case 1:
			nucs[i-1] = 'C'
		case 2:
			nucs[i-1] = 'G'
		case 3:
			nucs[i-1] = 'T'
		}
	}
	return nucs
}

func (kmer *Kmer) reverse_nucleotides2(k uint32) []byte {
	nucs := make([]byte, k)
	for i := uint32(0); i < k; i++ {
		j := i % 32
		mask := uint64(3 << (j * 2))
		wordIndex := (k - i - 1) / 32
		switch mask & kmer.BinaryKmer[wordIndex] >> (2 * j) {
		case 0:
			nucs[k-i-1] = 'A'
		case 1:
			nucs[k-i-1] = 'C'
		case 2:
			nucs[k-i-1] = 'G'
		case 3:
			nucs[k-i-1] = 'T'
		}
	}
	return nucs
}

func Open(filename string) (binary CortexVarBinary, err error) {
	var bin CortexVarBinary

	// Create the io.Reader
	fi, err := os.Open(filename)
	defer fi.Close()

	// Check if the file exists and that it is readable
	if err != nil {
		return bin, err
	}

	// Create a buffered io.Reader object.
	bin.reader = bufio.NewReader(fi)

	headerErr := bin.readHeader()
	if headerErr != nil {
		return bin, headerErr
	}

	bin.kmerStartFileOffset, _ = fi.Seek(0, os.SEEK_CUR)

	return bin, nil
}
