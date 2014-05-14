package main

import (
	"bufio"
	"bytes"
	"encoding/binary"
	"fmt"
	"io"
	"os"
)

type CortexHeader struct {
	top     cortexHeaderTop
	colours []colourMetaData
}

type colourMetaData struct {
	meanReadLength                  uint32
	totalSequenceLength             uint64
	name                            string
	errorRate                       float64
	topClipping                     bool
	removeLowCovSupernodes          bool
	removeLowCovKmers               bool
	cleanedAgainstGraph             bool
	removeLowCovSupernodesThreshold uint32
	removeLowCovKmersThreshold      uint32
	cleaningGraphName               string
}

type cortexHeaderTop struct {
	MagicChars    [6]byte
	VersionNumber uint32
	KmerSize      uint32
	WordsPerKmer  uint32
	ColourCount   uint32
}

type simpleKmer struct {
	BinaryKmer uint64
	Coverages  [3]uint32
	Edges      [3]byte
}

func (kmer *simpleKmer) nucs(k uint) []byte {
	nucs := make([]byte, k)
	mask := uint64(3 << (k*2 - 2))
	for i := uint(0); i < k; i++ {
		switch mask & kmer.BinaryKmer >> ((k - 1 - i) * 2) {
		case 0:
			nucs[i] = 'A'
		case 1:
			nucs[i] = 'C'
		case 2:
			nucs[i] = 'G'
		case 3:
			nucs[i] = 'T'
		}
		mask = mask >> 2
	}
	return nucs
}

func ReadHeader(r io.Reader) CortexHeader {
	var header CortexHeader
	binary.Read(r, binary.LittleEndian, &header.top)

	header.colours = make([]colourMetaData, header.top.ColourCount)

	for i := range header.colours {
		binary.Read(r, binary.LittleEndian, &header.colours[i].meanReadLength)
	}

	for i := range header.colours {
		binary.Read(r, binary.LittleEndian, &header.colours[i].totalSequenceLength)
	}

	var nameLength uint32
	var singleChar byte
	for i := range header.colours {
		binary.Read(r, binary.LittleEndian, &nameLength)
		b := make([]byte, nameLength)
		for j := range b {
			binary.Read(r, binary.LittleEndian, &singleChar)
			b[j] = singleChar
		}
		header.colours[i].name = string(b[:])
	}

	// Skip over the error rates encoded as unsigned long double.
	var errorRate float64
	for i := 0; i < len(header.colours); i++ {
		binary.Read(r, binary.LittleEndian, &errorRate)
		binary.Read(r, binary.LittleEndian, &errorRate)
	}

	var topClipping byte
	var removeLowCovSupernodes byte
	var removeLowCovKmers byte
	var cleanedAgainstGraph byte
	for i := range header.colours {
		colour := header.colours[i]
		binary.Read(r, binary.LittleEndian, &topClipping)
		switch topClipping {
		case 0x1:
			colour.topClipping = true
		case 0x0:
			colour.topClipping = false
		default:
			panic("topClipping neither 0x1 nor 0x0")
		}

		binary.Read(r, binary.LittleEndian, &removeLowCovSupernodes)
		switch removeLowCovSupernodes {
		case 0x1:
			colour.removeLowCovSupernodes = true
		case 0x0:
			colour.removeLowCovSupernodes = false
		default:
			panic("removeLowCovSupernodes neither 0x1 nor 0x0")
		}

		binary.Read(r, binary.LittleEndian, &removeLowCovKmers)
		switch removeLowCovKmers {
		case 0x1:
			colour.removeLowCovKmers = true
		case 0x0:
			colour.removeLowCovKmers = false
		default:
			panic("removeLowCovKmers neither 0x1 nor 0x0")
		}

		binary.Read(r, binary.LittleEndian, &cleanedAgainstGraph)
		switch cleanedAgainstGraph {
		case 0x1:
			colour.cleanedAgainstGraph = true
		case 0x0:
			colour.cleanedAgainstGraph = false
		default:
			panic("cleanedAgainstGraph neither 0x1 nor 0x0")
		}

		binary.Read(r, binary.LittleEndian, &colour.removeLowCovSupernodesThreshold)
		binary.Read(r, binary.LittleEndian, &colour.removeLowCovKmersThreshold)
		var cleaningGraphNameLength uint32
		binary.Read(r, binary.LittleEndian, &cleaningGraphNameLength)
		b := make([]byte, nameLength)
		for i := uint32(0); i < cleaningGraphNameLength; i++ {
			binary.Read(r, binary.LittleEndian, &b[i])
		}
		colour.cleaningGraphName = string(b)
		header.colours[i] = colour
	}

	// The last "CORTEX" characters
	var c [6]byte
	binary.Read(r, binary.LittleEndian, &c)
	return header
}

func (header CortexHeader) String() string {
	var b bytes.Buffer
	fmt.Fprintf(&b, "Magic chars:       %s\n", header.top.MagicChars)
	fmt.Fprintf(&b, "Version:           %d\n", header.top.VersionNumber)
	fmt.Fprintf(&b, "KmerSize:          %d\n", header.top.KmerSize)
	fmt.Fprintf(&b, "Words per kmer:    %d\n", header.top.WordsPerKmer)
	fmt.Fprintf(&b, "Number of colours: %d\n", header.top.ColourCount)
	return b.String()
}

func main() {
	fi, _ := os.Open("test/data/all_colours.ctx")
	defer fi.Close()
	r := bufio.NewReader(fi)
	header := ReadHeader(r)
	k := header.top.KmerSize
	var kmer simpleKmer
	for err := binary.Read(r, binary.LittleEndian, &kmer); err == nil; {
		fmt.Println(string(kmer.nucs(uint(k))))
		err = binary.Read(r, binary.LittleEndian, &kmer)
	}
}
