package cortex

import (
	"bufio"
	"bytes"
	"encoding/binary"
	"errors"
	"fmt"
	"io"
	"math/big"
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
	return binary.Read(bin.reader, binary.LittleEndian, data)
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

	if magicString != [6]byte{'C', 'O', 'R', 'T', 'E', 'X'} {
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
	bits          *KmerBits
	BinaryKmer    []uint64
	Coverages     []uint32
	ColouredEdges []edges
}

func NewKmer() (kmer Kmer) {
	return Kmer{bits: NewKmerBits()}
}

func NewKmerFromSequence(seq string) Kmer {
	kmer := NewKmer()
	var tmpByte byte
	var b bytes.Buffer
	for i, v := range seq {
		switch v {
		case 'C':
			newBits := byte(1) << (6 - (uint(i%4) * 2))
			tmpByte = tmpByte | newBits
		case 'G':
			newBits := byte(2) << (6 - (uint(i%4) * 2))
			tmpByte = tmpByte | newBits
		case 'T':
			newBits := byte(3) << (6 - (uint(i%4) * 2))
			tmpByte = tmpByte | newBits
		}
		if i%4 == 3 {
			b.WriteByte(tmpByte)
			tmpByte = 0
		}
	}
	b.WriteByte(tmpByte)
	kmer.bits.SetBytes(b.Bytes())
	return kmer
}

func (kmer *Kmer) Cmp(other *Kmer) int {
	return other.bits.Int.Cmp(kmer.bits.Int)
}

func (kmer *Kmer) Nucleotides() string {
	return "fail"
}

func (kmer *Kmer) nucleotides(k uint32) []byte {
	nucs := make([]byte, k)
	for i := k; i > 0; i-- {
		j := (i - 1) % 32
		wordIndex := (i - 1) / 32
		mask := uint64(3 << (2 * j))
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

func (kmer *Kmer) allEdges() (all edges) {
	all = edges(0)
	for _, edges := range kmer.ColouredEdges {
		all = all | edges
	}
	return all
}

func (kmer *Kmer) LeftKmers(k uint32) []Kmer {
	wordCount := len(kmer.BinaryKmer)

	// Shift the existing bitsting to the right
	shiftedBinaryKmer := make([]uint64, wordCount)
	mem := uint64(0)
	for i := range kmer.BinaryKmer {
		bitString := kmer.BinaryKmer[wordCount-i-1]
		shiftedBinaryKmer[wordCount-i-1] = (bitString >> 2) | (mem << 62)
		mem = bitString & 3
	}

	allEdges := kmer.allEdges()
	baseBytes := []uint64{
		uint64(0), // A
		uint64(1), // C
		uint64(2), // G
		uint64(3), // T
	}

	incomingKmers := make([]Kmer, 0)
	for i, bits := range baseBytes {
		if allEdges>>(4+uint(i))%2 == 1 {
			newKmer := Kmer{BinaryKmer: make([]uint64, wordCount)}
			copy(newKmer.BinaryKmer, shiftedBinaryKmer)
			newKmer.BinaryKmer[wordCount-1] = shiftedBinaryKmer[wordCount-1] | bits<<((k%32-1)*2)
			incomingKmers = append(incomingKmers, newKmer)
		}
	}
	return incomingKmers
}

func (kmer *Kmer) RightKmers(k uint32) []Kmer {
	wordCount := len(kmer.BinaryKmer)

	// Shift the existing bits to the left
	shiftedBinaryKmer := make([]uint64, wordCount)
	mem := uint64(0)
	for i := range kmer.BinaryKmer {
		bitString := kmer.BinaryKmer[i]
		shiftedBinaryKmer[i] = (bitString << 2) | (mem >> 62)
		mem = (bitString & 0xC000000000000000)
	}

	// Clean up the two bits now hanging on the left-hand edge of
	// the bitstring.
	shiftedBinaryKmer[wordCount-1] = shiftedBinaryKmer[wordCount-1] & (1<<((k%32)*2) - 1)

	allEdges := kmer.allEdges()
	baseBytes := []uint64{
		uint64(3), // T
		uint64(2), // G
		uint64(1), // C
		uint64(0), // A
	}

	outGoingKmers := make([]Kmer, 0)

	for i, bits := range baseBytes {
		if allEdges>>uint(i)%2 == 1 {
			newKmer := Kmer{BinaryKmer: make([]uint64, wordCount)}
			copy(newKmer.BinaryKmer, shiftedBinaryKmer)
			newKmer.BinaryKmer[0] = shiftedBinaryKmer[0] | bits
			outGoingKmers = append(outGoingKmers, newKmer)
		}
	}

	return outGoingKmers
}

type KmerBits struct {
	*big.Int
}

func NewKmerBits() *KmerBits {
	return &KmerBits{big.NewInt(0)}
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

	// Let's remember where the kmers start so that we can quickly
	// rewind if we have to.

	bin.kmerStartFileOffset, _ = fi.Seek(0, os.SEEK_CUR)

	return bin, nil
}
