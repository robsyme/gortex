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
	Edges      [3]edge
}

type edge byte

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

func (kmer *simpleKmer) reverse_nucs(k uint) []byte {
	nucs := make([]byte, k)
	mask := uint64(3)
	for i := uint(0); i < k; i++ {
		switch mask & kmer.BinaryKmer >> (2 * i) {
		case 0:
			nucs[i] = 'A'
		case 1:
			nucs[i] = 'C'
		case 2:
			nucs[i] = 'G'
		case 3:
			nucs[i] = 'T'
		}
		mask = mask << 2
	}
	return nucs
}

func AppendIfMissing(slice []byte, c byte) []byte {
	for _, ele := range slice {
		if ele == c {
			return slice
		}
	}
	return append(slice, c)
}

func (kmer simpleKmer) String() string {
	var b bytes.Buffer
	outgoing_edges := make([]byte, 0)
	for _, edges := range kmer.Outgoing() {
		for _, edge := range edges {
			outgoing_edges = AppendIfMissing(outgoing_edges, edge)
		}
	}

	for _, edge := range outgoing_edges {
		fmt.Fprintf(&b, "%s\t%c\n", kmer.nucs(21), edge)
	}
	return b.String()
}

func (kmer *simpleKmer) Incoming() [][]byte {
	incoming := make([][]byte, len(kmer.Edges))
	for i, v := range kmer.Edges {
		incoming[i] = v.incoming()
	}
	return incoming
}

func (kmer *simpleKmer) Outgoing() [][]byte {
	outgoing := make([][]byte, len(kmer.Edges))
	for i, v := range kmer.Edges {
		outgoing[i] = v.outgoing()
	}
	return outgoing
}

func (kmer *simpleKmer) PrintForColex() {
	var b bytes.Buffer
	outgoing_edges := make([]byte, 0)
	for _, edges := range kmer.Outgoing() {
		for _, edge := range edges {
			outgoing_edges = AppendIfMissing(outgoing_edges, edge)
		}
	}

	for _, edge := range outgoing_edges {
		fmt.Fprintf(&b, "%s\t%c\n", kmer.reverse_nucs(21), edge)
	}

	if kmer.hasNoIncoming() {
		root := make([]byte, 21)
		for i := range root {
			root[i] = '$'
		}
		nucs := kmer.reverse_nucs(21)
		for i := 1; i <= 21; i++ {
			fmt.Printf("%s%s\t%c\n", nucs[i:], root[0:i], nucs[i-1])
		}
	}

	fmt.Print(b.String())
}

func (kmer simpleKmer) hasNoIncoming() bool {
	for _, edge := range kmer.Edges {
		masked := (edge & 0xF0) >> 4
		if masked > 0 {
			return false
		}
	}
	return true
}

func (e edge) incoming() []byte {
	masked := (e & 0xF0) >> 4
	count := 0
	for i := uint(0); i < 4; i++ {
		if masked>>i%2 == 1 {
			count += 1
		}
	}

	edges := make([]byte, count)
	j := 0
	for i := uint8(0); i < 4; i++ {
		if masked>>i%2 == 1 {
			switch i {
			case 0:
				edges[j] = 'T'
			case 1:
				edges[j] = 'G'
			case 2:
				edges[j] = 'C'
			case 3:
				edges[j] = 'A'
			}
			j++
		}
	}
	return edges
}

func (e edge) outgoing() []byte {
	masked := e & 0x0F
	count := 0
	for i := uint(0); i < 4; i++ {
		if masked>>i%2 == 1 {
			count += 1
		}
	}

	edges := make([]byte, count)
	j := 0
	for i := uint8(0); i < 4; i++ {
		if masked>>i%2 == 1 {
			switch i {
			case 0:
				edges[j] = 'T'
			case 1:
				edges[j] = 'G'
			case 2:
				edges[j] = 'C'
			case 3:
				edges[j] = 'A'
			}
			j++
		}
	}
	return edges
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
	ReadHeader(r)
	var kmer simpleKmer
	for err := binary.Read(r, binary.LittleEndian, &kmer); err == nil; {
		kmer.PrintForColex()
		err = binary.Read(r, binary.LittleEndian, &kmer)
	}
}
