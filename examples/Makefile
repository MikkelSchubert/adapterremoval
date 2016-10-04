
TARGETS := single_end paired_end phred_64_encoded phred_33_encoded \
	collapsed_q50 multiple_adapters identify_adapters \
	interleaved_input interleaved_output interleaved_input_and_output \
	demultiplexing demultiplexing_only

.PHONY: all ${TARGETS}

all: ${TARGETS}


single_end:
	AdapterRemoval --file1 reads_1.fq --basename output_single --trimns --trimqualities --gzip

paired_end:
	AdapterRemoval --file1 reads_1.fq --file2 reads_2.fq --basename output_paired --trimns --trimqualities --collapse

interleaved_input:
	AdapterRemoval --interleaved-input --file1 interleaved.fq --basename output_interleaved_input

interleaved_output:
	AdapterRemoval --interleaved-output --file1 reads_1.fq --file2 reads_2.fq --basename output_interleaved_output

interleaved_input_and_output:
	AdapterRemoval --interleaved --file1 interleaved.fq --basename output_interleaved_both

phred_64_encoded:
	AdapterRemoval --qualitybase 64 --file1 reads_q64.fq --basename output_phred_64

phred_33_encoded:
	AdapterRemoval --qualitybase 64 --qualitybase-output 33 --file1 reads_q64.fq --basename output_phred_33

collapsed_q50:
	AdapterRemoval --qualitymax 50 --file1 reads_1.fq --file2 reads_2.fq --collapse --basename output_collapsed_q50

multiple_adapters:
	AdapterRemoval --file1 reads_1.fq --file2 reads_2.fq --basename output_multi --trimns --trimqualities --collapse --adapter-list adapters.txt

identify_adapters:
	AdapterRemoval --identify-adapters --file1 reads_1.fq --file2 reads_2.fq

demultiplexing:
	AdapterRemoval --file1 demux_1.fq --file2 demux_2.fq --basename output_demux --barcode-list barcodes.txt

demultiplexing_only:
	AdapterRemoval --file1 demux_1.fq --file2 demux_2.fq --basename output_only_demux --barcode-list barcodes.txt --demultiplex-only

