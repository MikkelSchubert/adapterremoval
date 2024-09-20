EXE := adapterremoval3

TARGETS := single_end paired_end interleaved_input interleaved_output \
	interleaved_input_and_output phred_64_encoded merged multiple_adapters \
	identify_adapters demultiplexing demultiplexing_only

.PHONY: all ${TARGETS} clean

all: ${TARGETS}

clean:
	$(RM) output_*

single_end:
	$(EXE) --in-file1 reads_1.fastq --out-prefix output_single

paired_end:
	$(EXE) --in-file1 reads_1.fastq --in-file2 reads_2.fastq --out-prefix output_paired

interleaved_input:
	$(EXE) --interleaved-input --in-file1 interleaved.fastq --out-prefix output_interleaved_input

interleaved_output:
	$(EXE) --interleaved-output --in-file1 reads_1.fastq --in-file2 reads_2.fastq --out-prefix output_interleaved_output

interleaved_input_and_output:
	$(EXE) --interleaved --in-file1 interleaved.fastq --out-prefix output_interleaved_both

phred_64_encoded:
	$(EXE) --quality-format 64 --in-file1 reads_q64.fastq --out-prefix output_phred_64

merged:
	$(EXE) --merge-quality-max 50 --in-file1 reads_1.fastq --in-file2 reads_2.fastq --merge --out-prefix output_merged

multiple_adapters:
	$(EXE) --in-file1 reads_1.fastq --in-file2 reads_2.fastq --out-prefix output_multi --adapter-list adapters.txt

identify_adapters:
	$(EXE) --report-only --in-file1 reads_1.fastq --in-file2 reads_2.fastq --out-prefix output_adapters

demultiplexing:
	$(EXE) --in-file1 demux_1.fastq --in-file2 demux_2.fastq --out-prefix output_demux --barcode-list barcodes.txt

demultiplexing_only:
	$(EXE) --in-file1 demux_1.fastq --in-file2 demux_2.fastq --out-prefix output_only_demux --barcode-list barcodes.txt --demultiplex-only
