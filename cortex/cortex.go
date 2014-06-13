package cortex

import (
	"bufio"
	"bytes"
	"encoding/binary"
	"errors"
	"github.com/robsyme/multifinder/colouredKmer"
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
	extraBytes          uint32
	KmerCount           int64
	fi                  *os.File
}

func (bin *CortexVarBinary) NextKmer() (*colouredKmer.Kmer, error) {
	kmer := colouredKmer.Kmer{
		Int:       big.NewInt(0),
		Coverages: make([]uint32, bin.ColourCount),
		Edges:     make([]uint8, bin.ColourCount),
	}

	bits := make([]uint64, bin.wordsPerKmer)
	err := bin.Read(&bits)
	if err != nil {
		return &kmer, err
	}
	words := make([]big.Word, bin.wordsPerKmer)
	for i, v := range bits {
		words[i] = big.Word(v)
	}
	kmer.SetBits(words)
	bin.Read(kmer.Coverages)
	bin.Read(kmer.Edges)
	return &kmer, nil
}

func (bin *CortexVarBinary) Kmers() chan *colouredKmer.Kmer {
	bin.fi.Seek(bin.kmerStartFileOffset, os.SEEK_SET)
	kc := make(chan *colouredKmer.Kmer)
	go func() {
		for k, err := bin.NextKmer(); err == nil; k, err = bin.NextKmer() {
			kc <- k
		}
		close(kc)
	}()
	return kc
}

func (bin *CortexVarBinary) Read(data interface{}) error {
	return binary.Read(bin.reader, binary.LittleEndian, data)
}

func (bin *CortexVarBinary) readHeader() error {
	bin.fi.Seek(0, os.SEEK_SET)
	// Create a buffered io.Reader object.
	bin.reader = bufio.NewReader(bin.fi)

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
		bin.extraBytes += nameLength
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
		bin.extraBytes += nameLength
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

func (bin *CortexVarBinary) Close() {
	bin.fi.Close()
}

func Open(filename string) (binary CortexVarBinary, err error) {
	var bin CortexVarBinary

	// Create the io.Reader
	bin.fi, err = os.Open(filename)
	if err != nil {
		return bin, err
	}

	headerErr := bin.readHeader()
	if headerErr != nil {
		return bin, headerErr
	}

	fileInfo, _ := bin.fi.Stat()
	headerSize := bin.extraBytes + 28 + 48*bin.ColourCount
	bin.KmerCount = (fileInfo.Size() - int64(headerSize)) / 23

	// Let's remember where the kmers start so that we can quickly
	// rewind if we have to.
	bin.kmerStartFileOffset, _ = bin.fi.Seek(0, os.SEEK_CUR)
	return bin, nil
}
