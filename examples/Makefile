EXE := adapterremoval3

TARGETS := single_end paired_end phred_64_encoded \
	merged_q50 multiple_adapters identify_adapters \
	interleaved_input interleaved_output interleaved_input_and_output \
	demultiplexing demultiplexing_only

.PHONY: all ${TARGETS} clean

all: ${TARGETS}

clean:
	$(RM) output_*

single_end:
	$(EXE) --file1 reads_1.fastq --basename output_single --gzip

paired_end:
	$(EXE) --file1 reads_1.fastq --file2 reads_2.fastq --basename output_paired

interleaved_input:
	$(EXE) --interleaved-input --file1 interleaved.fastq --basename output_interleaved_input

interleaved_output:
	$(EXE) --interleaved-output --file1 reads_1.fastq --file2 reads_2.fastq --basename output_interleaved_output

interleaved_input_and_output:
	$(EXE) --interleaved --file1 interleaved.fastq --basename output_interleaved_both

phred_64_encoded:
	$(EXE) --quality-format 64 --file1 reads_q64.fastq --basename output_phred_64

merged:
	$(EXE) --qualitymax 50 --file1 reads_1.fastq --file2 reads_2.fastq --merge --basename output_merged

multiple_adapters:
	$(EXE) --file1 reads_1.fastq --file2 reads_2.fastq --basename output_multi --adapter-list adapters.txt

identify_adapters:
	$(EXE) --report-only --file1 reads_1.fastq --file2 reads_2.fastq --basename output_adapters

demultiplexing:
	$(EXE) --file1 demux_1.fastq --file2 demux_2.fastq --basename output_demux --barcode-list barcodes.txt

demultiplexing_only:
	$(EXE) --file1 demux_1.fastq --file2 demux_2.fastq --basename output_only_demux --barcode-list barcodes.txt --demultiplex-only
