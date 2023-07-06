package main

import (
	"encoding/csv"
	"flag"
	"fmt"
	"io"
	"log"
	"os"
	"strconv"
	"strings"
	"sync"
	"bufio"
)

var mu sync.Mutex

type AlignData struct {
	refLen     int
	pos        int
	meanCounts float64
}

// readCSV reads a CSV file and processes it into a map where each key is a header profile and each value is a slice of alignment data.
// The function takes two parameters: the path to the CSV file and a boolean flag noZeros. If noZeros is true, any replicate with 0 reads aligned is ignored.
// It returns the length of the small RNA (srnaLen) and the map of header profiles to alignment data.
// readCSV handles the CSV file using the encoding/csv package. It reads and processes each record one by one, skipping the header.
func readCSV(inCSV string, noZeros bool) (int, map[string][]AlignData) {
	file, _ := os.Open(inCSV)
	defer file.Close()

	r := csv.NewReader(file)
	headerProfiles := make(map[string][]AlignData)

	// Get the length of the small RNA from the file name
	parts := strings.Split(inCSV, "_")
	last := parts[len(parts)-1]
	sirnaLenStr := strings.Split(last, ".")[0]
	sirnaLen, _ := strconv.Atoi(sirnaLenStr)

	// Read and process records one by one
	_, _ = r.Read() // Skip header
	for {
		record, err := r.Read()
		if err == io.EOF {
			break
		}

		counts := make([]float64, len(record[6:]))

		for i, countStr := range record[6:] {
			count, _ := strconv.ParseFloat(countStr, 64)
			counts[i] = count
		}

		if noZeros {
			hasZero := false
			for _, count := range counts {
				if count == 0.0 {
					hasZero = true
					break
				}
			}

			if hasZero {
				continue
			}
		}

		refLen, _ := strconv.Atoi(record[1])
		row3, _ := strconv.Atoi(record[3])
		meanCounts := mean(counts)

		rowData := AlignData{
			refLen:     refLen,
			pos:        row3,
			meanCounts: meanCounts,
		}

		key := strings.Split(record[0], " ")[0]
		headerProfiles[key] = append(headerProfiles[key], rowData)
	}

	return sirnaLen, headerProfiles
}

// mean calculates the mean of a slice of float64 values.
// It takes a slice of float64 values as a parameter and returns the mean as a float64.
// The mean is calculated as the sum of the values divided by the number of values.
func mean(floats []float64) float64 {
	sum := 0.0
	for _, number := range floats {
		sum += number
	}
	return sum / float64(len(floats))
}

type Result struct {
	Data []float64
}

// prepareDataForWriting prepares the data for writing to a WIG file.
// It takes a header (chromosome name) and the corresponding result data as parameters.
// It formats the data into a slice of strings, each string representing a line to be written to the WIG file.
// It returns the slice of strings.
func prepareDataForWriting(header string, result Result) []string {
	var lines []string
	lines = append(lines, fmt.Sprintf("variableStep chrom=%s", header))
	pos := 1

	for _, value := range result.Data {
		if value > 0.0 {
			lines = append(lines, fmt.Sprintf("%d %f", pos, value))
		}
		pos++
	}
	return lines
}


// writeWig writes a slice of lines to a WIG file.
// It takes the output file path and the slice of lines as parameters.
// It opens the file in append mode and writes the lines to it.
func writeWig(outWig string, lines []string) {
	// Obtain the lock before writing to the file
	mu.Lock()
	defer mu.Unlock()

	file, _ := os.OpenFile(outWig, os.O_APPEND|os.O_CREATE|os.O_WRONLY, 0644)
	defer file.Close()

	// Create a new Writer
	w := bufio.NewWriter(file)

	for _, line := range lines {
		_, _ = w.WriteString(line + "\n")
	}

	// Don't forget to flush!
	w.Flush()
}

// handleFlags handles the command line flags and returns the flag values.
// It checks whether the input and output file paths are provided and removes the output file if it already exists.
func handleFlags() (string, string, int, bool) {
	inCSV := flag.String("input", "", "Input CSV file (required)")
	outWig := flag.String("output", "", "Output WIG file (required)")
	batchSize := flag.Int("batch", 16, "Batch size (optional, default is 16)")
	noZeros := flag.Bool("noZeros", false, "No Zeros flag -- if used, all replicates must have > 0 reads aligned for the position to be added (optional, default is false)")

	flag.Parse()

	if *inCSV == "" || *outWig == "" {
		flag.PrintDefaults()
		os.Exit(1)
	}

	if err := os.Remove(*outWig); err != nil && !os.IsNotExist(err) {
		log.Fatalf("Failed to delete existing file: %v", err)
	}

	return *inCSV, *outWig, *batchSize, *noZeros
}

// processData reads the input CSV file, processes the data in batches and writes the output to a WIG file.
// It divides the keys (header) into batches and calls writeData function for each batch.
func processData(inCSV string, noZeros bool, batchSize int, outWig string) {
	srnaLen, headerProfiles := readCSV(inCSV, noZeros)
	fmt.Println("Read length =", srnaLen)

	keys := make([]string, 0, len(headerProfiles))
	for k := range headerProfiles {
		keys = append(keys, k)
	}

	for i := 0; i < len(keys); i += batchSize {
		j := i + batchSize
		if j > len(keys) {
			j = len(keys)
		}
		fmt.Println("Processing and writing batch ", i, " to ", j)
		writeData(keys[i:j], headerProfiles, srnaLen, outWig)
	}
}

// writeData takes a slice of keys, the headerProfiles map, srnaLen and output file path as parameters.
// It creates a goroutine for each key in the keys slice to prepare the data for writing and write it to the WIG file.
// It waits for all the goroutines to finish before returning.
func writeData(keys []string, headerProfiles map[string][]AlignData, srnaLen int, outWig string) {
	var wg sync.WaitGroup
	for _, header := range keys {
		alignments := headerProfiles[header]
		arr := make([]float64, alignments[0].refLen)

		for _, alignment := range alignments[0:] {
			for i := 0; i < srnaLen; i++ {
				if alignment.pos-1+i < len(arr) {
					arr[alignment.pos-1+i] += alignment.meanCounts
				}
			}
		}

		wg.Add(1)
		go func(header string, result Result) {
			defer wg.Done()
			lines := prepareDataForWriting(header, result)
			writeWig(outWig, lines)
		}(header, Result{Data: arr})
	}
	wg.Wait()
}

// main function is the entry point of the program.
// It handles the command-line flags and calls the processData function.
func main() {
	inCSV, outWig, batchSize, noZeros := handleFlags()
	fmt.Println("Reading input CSV")
	processData(inCSV, noZeros, batchSize, outWig)
}
