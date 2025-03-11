#!/usr/bin/python3
# /// script
# requires-python = ">=3.9"
# dependencies = [
#     "isal",
#     "tqdm",
#     "typed-argparse",
# ]
# [tool.ruff]
# target-version = "py39"
# [tool.uv]
# exclude-newer = "2025-02-25T00:00:00Z"
# ///
from __future__ import annotations

import json
import math
import sys
from collections.abc import Iterator, Sequence
from functools import partial
from itertools import chain, islice, repeat
from pathlib import Path
from random import Random
from typing import IO, Iterable, NoReturn, Optional, TypeVar, cast

import typed_argparse as tap
from isal import igzip
from tqdm import tqdm

T = TypeVar("T")


def parse_count(value: str) -> int:
    if value.isdigit():
        return int(value)
    elif value and value[-1] in "kmgKMG" and value[:-1].isdigit():
        count = int(value[:-1])
        if value[-1] in "kK":
            return count * 1_000
        elif value[-1] in "mM":
            return count * 1_000_000
        elif value[-1] in "gG":
            return count * 1_000_000_000
        raise NotImplementedError(value)
    else:
        raise ValueError(value)


def parse_rate(value: str) -> float:
    rate = float(value)
    if not (0.0 <= rate <= 1.0):
        raise ValueError(value)

    return rate


class Args(tap.TypedArgs):
    output_prefix: Optional[str] = tap.arg(
        help="Prefix for output filenames",
        positional=True,
    )

    # Optional is required for typed-argparse with older Python versions
    seed: Optional[int] = tap.arg(help="Seed used for RNG")  # noqa: UP007

    adapter1: str = tap.arg(
        default="AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG",
        help="Adapter added to the 5' end of the forward strand (read from 5' ...). Ns "
        "are replaced with random bases",
    )
    adapter2: str = tap.arg(
        default="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT",
        help="Adapter added to the 5' end of the reverse strand (read from 3' ...) Ns "
        "are replaced with random bases",
    )

    adapter_list: Optional[Path] = tap.arg(  # noqa: UP007
        help="Two-column table of forward/reverse adapters in the same format as "
        "expected by AdapterRemoval. Ns are replaced with random bases on load. The "
        "--adapter1 and --adapter2 options are ignored if this option is used"
    )
    barcode_list: Optional[Path] = tap.arg(  # noqa: UP007
        help="Three-column table of sample/forward/reverse barcodes in the same format "
        "as expected by AdapterRemoval. Ns are replaced with random bases on load"
    )
    reversible_barcodes: bool = tap.arg(
        help="Generate barcodes in both forward (user supplied) and reverse orientation"
    )

    frag_len_mu: int = tap.arg(
        default=150,
        help="Mean of normal distribution used to generate DNA fragment lengths",
    )
    frag_len_sigma: int = tap.arg(
        default=30,
        help="Variance of normal distribution used to generate DNA fragment lengths",
    )
    frag_len_min: int = tap.arg(
        default=0,
        type=lambda x: max(0, int(x)),
        help="Minimum length of generated DNA fragments",
    )
    frag_len_max: int = tap.arg(
        default=500,
        help="Maximum length of generated DNA fragments",
    )

    # Rates inspired by https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-0976-y
    reads_r1_sub_rate: float = tap.arg(
        default=0.0021,
        metavar="P",
        help="The absolute or relative (with templates) substitution rate for r1 reads",
    )
    reads_r2_sub_rate: float = tap.arg(
        default=0.0042,
        type=parse_rate,
        help="The absolute or relative (with templates) substitution rate for r2 reads",
    )
    reads_insert_rate: float = tap.arg(
        default=8.15e-6,
        type=parse_rate,
        help="The absolute or relative (with templates) insertion rate",
    )
    reads_deletion_rate: float = tap.arg(
        default=5.00e-06,
        type=parse_rate,
        help="The absolute or relative (with templates) deletion rate",
    )
    reads_indel_ext_rate: float = tap.arg(
        default=0.4,
        type=parse_rate,
        help="(Repeat) probability of extending an indel",
    )
    reads_len: int = tap.arg(
        default=150,
        help="Length of reads in bp; note that shorter reads may be generated if the "
        "quality-templates are shorter than this value",
    )

    pairs: int = tap.arg(
        default=10_000,
        type=parse_count,
        help="The number of read pairs to generate; units K, M, G are supported",
    )

    quality_1_template: Optional[Path] = tap.arg(
        metavar="FASTQ",
        help="Use Phred scores from this FASTQ file for mate 1 reads. The rate "
        "arguments above are used to determine relative rates of substitutions, "
        "insertions, and deletions for a given Phred score",
    )
    quality_2_template: Optional[Path] = tap.arg(
        metavar="FASTQ",
        help="Use Phred scores from this FASTQ file for mate 2 reads. The rate "
        "arguments above are used to determine relative rates of substitutions, "
        "insertions, and deletions for a given Phred score",
    )
    randomize_templates: bool = tap.arg(
        help="Randomly select N templates instead of the first N templates, when using "
        "--quality-1-template and/or --quality-2-template"
    )

    write_meta: Optional[Path] = tap.arg(
        metavar="PATH",
        help="Write arguments, barcods, and adapters to the specified PATH as JSON",
    )


def abort(*values: object) -> NoReturn:
    print(*values, file=sys.stderr)
    sys.exit(1)


def reservoir_sampling(
    items: Iterable[T],
    downsample_to: int,
    rng: Random,
) -> list[T]:
    if not isinstance(downsample_to, int):
        raise TypeError(f"downsample_to must be an int, not {downsample_to!r}")
    elif downsample_to < 0:
        raise ValueError(f"downsample_to must be >= 0, not {downsample_to}")

    items = iter(tqdm(items, unit_scale=True, unit=" read"))
    reservoir: list[T] = list(islice(items, downsample_to))
    for index, item in enumerate(items, start=downsample_to):
        index = rng.randint(0, index)
        if index < downsample_to:
            reservoir[index] = item

    return reservoir


def indexed_choice(rng: Random, items: Sequence[T]) -> tuple[int, T]:
    idx = rng.randrange(len(items))
    return (idx, items[idx])


def build_complement_table() -> str:
    # Pairs of complementa  ry bases and ambiguous bases
    table = ["N"] * 256
    for nuc_a, nuc_b in ["AT", "CG"]:
        table[ord(nuc_a)] = nuc_b
        table[ord(nuc_b)] = nuc_a

    return "".join(table)


# Table for complementing sequences using `str.translate`
NT_COMPLEMENTS = build_complement_table()


def reverse_complement(sequence: str) -> str:
    """Returns the reverse complement of a DNA sequence."""
    return sequence.translate(NT_COMPLEMENTS)[::-1]


def round_to_int(value: float) -> int:
    return int(round(value))


def get_indel_length(indel_ext_rate: float, limit: int, rng: Random) -> int:
    n = 1
    while n < limit and rng.random() < indel_ext_rate:
        n += 1

    return n


def prepare_tmpl(rng: Random, seq: str) -> str:
    if not seq:
        raise ValueError("empty sequence")

    result: list[str] = []
    for char in seq.upper():
        if char == "N":
            result.append(rng.choice("ACGT"))
        elif char in "ACGT":
            result.append(char)
        else:
            raise ValueError(result)

    return "".join(result)


def calculate_cumulative_error_rates(
    *,
    sub_rate: float,
    ins_rate: float,
    del_rate: float,
) -> dict[str, tuple[float, float, float]]:
    error_rate = sub_rate + ins_rate + del_rate
    sub_rate /= error_rate
    ins_rate /= error_rate
    del_rate /= error_rate

    error_rates: dict[str, tuple[float, float, float]] = {}
    for q in range(0, 94):
        error_rate = min(0.75, 10 ** (q / -10))
        error_rates[chr(33 + q)] = (
            error_rate * (del_rate + ins_rate + sub_rate),
            error_rate * (del_rate + ins_rate),
            error_rate * del_rate,
        )

    return error_rates


_SUBSTITUTIONS = {"A": "CGT", "C": "AGT", "G": "ACT", "T": "ACG"}


def mutate_sequence(
    *refseqs: str,
    rng: Random,
    length: int,
    padding: str,
    indel_ext_rate: float,
    qualities: str,
    error_rates: dict[str, tuple[float, float, float]],
) -> str:
    subs = _SUBSTITUTIONS
    reflen = sum(map(len, refseqs))
    refseq = chain(*refseqs, repeat(*padding))
    current = next(refseq)

    sequence: list[str] = []
    while len(sequence) < length:
        phred_33 = qualities[min(len(sequence), len(qualities) - 1)]
        sub_rate, ins_rate, del_rate = error_rates[phred_33]

        event = rng.random()
        # In order of (typical) probability
        if event >= sub_rate:  # normal read
            sequence.append(current)
            current = next(refseq)
        elif event >= ins_rate:  # substitution
            sequence.append(rng.choice(subs[current]))
            current = next(refseq)
        elif event >= del_rate:  # insertion
            limit = length - len(sequence)
            k = get_indel_length(indel_ext_rate, limit, rng)
            sequence.extend(rng.choices("ACGT", k=k))
        else:  # deletion
            k = get_indel_length(indel_ext_rate, reflen, rng)
            for _ in range(k):
                current = next(refseq)

    return "".join(sequence[:length])


def get_random_sequence(rng: Random, args: Args) -> str:
    length = round_to_int(
        rng.gauss(
            args.frag_len_mu,
            args.frag_len_sigma,
        )
    )

    length = max(
        args.frag_len_min,
        min(args.frag_len_max, length),
    )

    # Hard-coded GC content of 40%
    return "".join(rng.choices("AAACCGGTTT", k=length))


def get_random_read_pair(
    rng: Random,
    args: Args,
    qualities_1: str,
    qualities_2: str,
    error_rates_1: dict[str, tuple[float, float, float]],
    error_rates_2: dict[str, tuple[float, float, float]],
    adapters: list[tuple[str, str]],
    barcodes: list[tuple[str, str, str]],
) -> tuple[str, str, str]:
    adapter_idx, (adapter1, adapter2) = indexed_choice(rng, adapters)
    sequence = get_random_sequence(rng, args)

    sample, barcode_idx, orientation, barcode1, barcode2 = ".", ".", ".", "", ""
    if barcodes:
        barcode_idx, (sample, barcode1, barcode2) = indexed_choice(rng, barcodes)
        if args.reversible_barcodes and rng.getrandbits(1):
            barcode1, barcode2 = barcode2, barcode1
            orientation = "-"
        else:
            orientation = "+"

    forward = mutate_sequence(
        barcode1,
        sequence,
        reverse_complement(barcode2),
        adapter1,
        rng=rng,
        qualities=qualities_1,
        error_rates=error_rates_1,
        indel_ext_rate=args.reads_indel_ext_rate,
        length=min(args.reads_len, len(qualities_1)),
        padding="A",
    )

    reverse = mutate_sequence(
        barcode2,
        reverse_complement(sequence),
        reverse_complement(barcode1),
        adapter2,
        rng=rng,
        qualities=qualities_2,
        error_rates=error_rates_2,
        indel_ext_rate=args.reads_indel_ext_rate,
        length=min(args.reads_len, len(qualities_2)),
        padding="T",
    )

    meta = f"{len(sequence)}:{sample}:{barcode_idx}:{orientation}:{adapter_idx}"
    return (meta, forward, reverse)


def read_table(path: Path, *, columns: int) -> Iterator[list[str]]:
    with path.open("r", encoding="utf-8") as handle:
        rows = 0
        for line in handle:
            line = line.rstrip()
            if not line or line.startswith("#"):
                continue

            fields = line.split()
            if len(fields) != columns:
                abort(
                    f"Wrong number of columns in {path}; expected {columns} columns, "
                    f"but found {len(fields)}: {fields}"
                )

            rows += 1
            yield fields

    if not rows:
        abort(f"No non-empty rows found in {path!r}")


def read_adapters(
    *,
    rng: Random,
    path: Path | None,
    adapter1: str,
    adapter2: str,
) -> list[tuple[str, str]]:
    adapters = [(adapter1, adapter2)]
    if path is not None:
        adapters = read_table(path, columns=2)

    return [(prepare_tmpl(rng, a1), prepare_tmpl(rng, a2)) for a1, a2 in adapters]


def read_barcodes(
    *,
    rng: Random,
    path: Path | None,
) -> list[tuple[str, str, str]]:
    if path is None:
        return []

    def _validate_name(name: str) -> str:
        if not all((c.isalnum() or c == "_") for c in name):
            raise ValueError(f"invalid sample-name {name!r}")

        return name

    return [
        (_validate_name(name), prepare_tmpl(rng, bc1), prepare_tmpl(rng, bc2))
        for name, bc1, bc2 in read_table(path, columns=3)
    ]


def read_or_generate_qualities(
    template: Path | None,
    *,
    length: int,
    default_error_rate: float,
    downsample_to: int | None = None,
    rng: Random,
) -> Iterator[str]:
    if template is not None:
        it = read_fastq(template, length=length)
        if downsample_to is not None:
            it = iter(reservoir_sampling(it, downsample_to=downsample_to, rng=rng))

        return it
    else:
        return quality_template(default_error_rate, length=length)


def read_fastq(path: Path, *, length: int) -> Iterator[str]:
    with igzip.open(path, "rt", encoding="utf-8") as handle:
        for line in islice(cast(IO[str], handle), 3, None, 4):
            yield line.rstrip()[:length]


def quality_template(error: float, *, length: int) -> Iterator[str]:
    if not (0 <= error <= 1):
        raise ValueError(error)

    # Phred+33 encoded error rate
    line = chr(33 + int(-10 * math.log10(error))) * length

    return repeat(line)


def file_or_stdout(tmpl: str, value: str | None) -> IO[str]:
    if value is None:
        return sys.stdout

    return open(tmpl.format(value), "w", encoding="utf-8")


def metadata_to_json(
    args: Args,
    adapters: list[tuple[str, str]],
    barcodes: list[tuple[str, str, str]],
) -> str:
    def _serialize(o: object) -> object:
        if isinstance(o, Path):
            o = str(o)

        return o

    metadata = {
        "args": sys.argv[1:],
        "params": {
            f"--{key.replace('_', '-')}": _serialize(value)
            for key, value in vars(args).items()
        },
        "adapters": [
            {"adapter1": adapter1, "adapter2": adapter2}
            for adapter1, adapter2 in adapters
        ],
        "barcodes": [
            {"sample": sample, "barcode1": barcode1, "barcode2": barcode2}
            for sample, barcode1, barcode2 in barcodes
        ],
    }

    return json.dumps(metadata, indent=2)


def main(args: Args) -> None:
    rng = Random(args.seed)

    adapters = read_adapters(
        rng=rng,
        path=args.adapter_list,
        adapter1=args.adapter1,
        adapter2=args.adapter2,
    )

    barcodes = read_barcodes(rng=rng, path=args.barcode_list)

    error_rates_1 = calculate_cumulative_error_rates(
        sub_rate=args.reads_r1_sub_rate,
        ins_rate=args.reads_insert_rate,
        del_rate=args.reads_deletion_rate,
    )

    error_rates_2 = calculate_cumulative_error_rates(
        sub_rate=args.reads_r2_sub_rate,
        ins_rate=args.reads_insert_rate,
        del_rate=args.reads_deletion_rate,
    )

    with (
        file_or_stdout("{}_1.fastq", args.output_prefix) as out_1,
        file_or_stdout("{}_2.fastq", args.output_prefix) as out_2,
    ):
        indel_rate = args.reads_insert_rate + args.reads_deletion_rate

        _read_or_generate_qualities = partial(
            read_or_generate_qualities,
            rng=rng,
            length=args.reads_len,
            downsample_to=args.pairs if args.randomize_templates else None,
        )

        phred_1 = _read_or_generate_qualities(
            template=args.quality_1_template,
            default_error_rate=args.reads_r1_sub_rate + indel_rate,
        )

        phred_2 = _read_or_generate_qualities(
            template=args.quality_2_template,
            default_error_rate=args.reads_r2_sub_rate + indel_rate,
        )

        for idx in tqdm(range(1, args.pairs + 1), unit_scale=True, unit=" pair"):
            try:
                qualities_1 = next(phred_1)
                qualities_2 = next(phred_2)
            except StopIteration:
                break

            meta, seq_1, seq_2 = get_random_read_pair(
                rng=rng,
                args=args,
                qualities_1=qualities_1,
                qualities_2=qualities_2,
                error_rates_1=error_rates_1,
                error_rates_2=error_rates_2,
                adapters=adapters,
                barcodes=barcodes,
            )

            out_1.write("@%s:%s:1\n%s\n+\n%s\n" % (idx, meta, seq_1, qualities_1))
            out_2.write("@%s:%s:2\n%s\n+\n%s\n" % (idx, meta, seq_2, qualities_2))

    if args.write_meta is not None:
        args.write_meta.write_text(metadata_to_json(args, adapters, barcodes))


if __name__ == "__main__":
    try:
        tap.Parser(Args).bind(main).run()
    except KeyboardInterrupt:
        sys.exit(1)  # FIXME
    except BrokenPipeError:
        pass
