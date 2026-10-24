"""FASTA input and numerical ACGT representation.

This is the single FASTA reader for the project. Every script should read
sequence through ``read_fasta`` so that file handling (gzip, line endings,
blank lines) behaves the same everywhere.

Encoding
--------
Bases are coded A=0, C=1, G=2, T=3 (upper or lower case; soft-masking is
ignored). Every other character -- N, IUPAC ambiguity codes, gaps ('-', '.'),
anything else -- is coded ``INVALID`` (255).

Gaps: the rule used throughout the project
------------------------------------------
A pattern must never be matched, counted or simulated across a non-ACGT
character or across a record boundary. Bases on either side of an N block were
not adjacent in the genome, so joining them creates sequence that does not
exist.

The recommended representation keeps one code per original position:

    codes = encode(seq)                  # full length, INVALID at non-ACGT
    for start, seg in iter_segments(codes, min_len=L):
        ...                              # seg is a pure ACGT view; the genome
                                         # coordinate of seg[i] is start + i

This preserves exact 0-based coordinates within the record, so any hit can be
reported as ``record:position`` and inspected later.

Functions
---------
read_fasta(path)          stream records one at a time (plain or gzip)
read_fasta_list(path)     same, as a list (small files only)
encode(seq)               bytes -> full-length uint8 codes, INVALID at non-ACGT
iter_segments(x, min_len) yield (start, view) for each ACGT stretch >= min_len
acgt_segments(x, min_len) list of the ACGT arrays between gaps (no offsets)
acgt_runs(seq, min_run)   (codes, [(start, length), ...])  -- original API
acgt_stream(seq, keep_gaps=False)
                          one concatenated stream. By default non-ACGT bases
                          are DELETED, which joins sequence across gaps; kept
                          only for backward compatibility. With keep_gaps=True
                          each gap is collapsed to one GAP code (4), so a
                          caller can split on it -- but coordinates are lost,
                          and binary projections must mask GAP explicitly.
                          Prefer encode() + iter_segments().
"""
from __future__ import annotations

import gzip
from typing import Iterator, List, NamedTuple, Tuple, Union

import numpy as np

INVALID = 255   # code for any non-ACGT character in encode()
GAP = 4         # separator code used by acgt_stream(keep_gaps=True)

CODE = np.full(256, INVALID, np.uint8)
for _i, _c in enumerate(b'ACGT'):
    CODE[_c] = _i
for _i, _c in enumerate(b'acgt'):
    CODE[_c] = _i

_GZIP_MAGIC = b'\x1f\x8b'


# ------------------------------------------------------------------- reading

class FastaRecord(NamedTuple):
    """One FASTA record. Unpacks as ``header, seq`` for backward compatibility.

    header : bytes, the header line without the leading '>' and line ending
    seq    : bytes, the sequence with all line breaks and whitespace removed,
             original characters and case preserved
    """
    header: bytes
    seq: bytes

    @property
    def name(self) -> str:
        """Record ID: the header up to the first whitespace (e.g. 'chr1')."""
        parts = self.header.split(None, 1)
        return parts[0].decode('utf-8', 'replace') if parts else ''


def _open(path: str):
    """Open plain or gzip-compressed files, detected by content, not extension."""
    with open(path, 'rb') as fh:
        magic = fh.read(2)
    return gzip.open(path, 'rb') if magic == _GZIP_MAGIC else open(path, 'rb')


def read_fasta(path: str) -> Iterator[FastaRecord]:
    """Yield FASTA records one at a time.

    Only one record is held in memory, so a whole genome can be streamed.
    Handles gzip input, Windows (\\r\\n) line endings and blank lines.
    Raises ValueError if sequence appears before the first header.
    """
    header = None
    chunks: List[bytes] = []
    with _open(path) as fh:
        for line in fh:
            if line.startswith(b'>'):
                if header is not None:
                    yield FastaRecord(header, b''.join(chunks))
                header = line[1:].rstrip(b'\r\n')
                chunks = []
            else:
                s = b''.join(line.split())      # drops \r, \n, spaces, tabs
                if not s:
                    continue
                if header is None:
                    raise ValueError(f'{path}: sequence found before the first FASTA header')
                chunks.append(s)
    if header is not None:
        yield FastaRecord(header, b''.join(chunks))


def read_fasta_list(path: str) -> List[FastaRecord]:
    """All records as a list. Convenience for small files; use read_fasta for genomes."""
    return list(read_fasta(path))


# ------------------------------------------------------------------ encoding

def encode(seq: bytes) -> np.ndarray:
    """Full-length uint8 codes: A=0 C=1 G=2 T=3, INVALID (255) for anything else.

    Position i of the result is position i of the record (0-based).
    """
    return CODE[np.frombuffer(seq, dtype=np.uint8)]


def _as_codes(x: Union[bytes, np.ndarray]) -> np.ndarray:
    return encode(x) if isinstance(x, (bytes, bytearray, memoryview)) else np.asarray(x)


def _valid_runs(codes: np.ndarray, min_len: int) -> Tuple[np.ndarray, np.ndarray]:
    """Start and length arrays of maximal ACGT stretches with length >= min_len."""
    ok = (codes < 4).astype(np.int8)
    d = np.diff(np.concatenate(([0], ok, [0])))
    s = np.flatnonzero(d == 1)
    e = np.flatnonzero(d == -1)
    keep = (e - s) >= max(min_len, 1)
    return s[keep], (e - s)[keep]


def iter_segments(x: Union[bytes, np.ndarray], min_len: int = 1) -> Iterator[Tuple[int, np.ndarray]]:
    """Yield ``(start, segment)`` for every ACGT stretch of length >= min_len.

    x       : raw sequence bytes, or codes from encode()
    min_len : skip stretches too short to hold a window of this length
    segment : a view into the code array (no copy), values 0..3 only;
              its first base is at 0-based position ``start`` in the record
    """
    codes = _as_codes(x)
    starts, lens = _valid_runs(codes, min_len)
    for s, n in zip(starts.tolist(), lens.tolist()):
        yield s, codes[s:s + n]


def acgt_segments(x: Union[bytes, np.ndarray], min_len: int = 1) -> List[np.ndarray]:
    """The ACGT arrays between gaps, in order, without their offsets."""
    return [seg for _, seg in iter_segments(x, min_len)]


def acgt_runs(seq_bytes: bytes, min_run: int = 1):
    """Original API: ``(codes, [(start, length), ...])`` for ACGT stretches >= min_run."""
    codes = encode(seq_bytes)
    s, n = _valid_runs(codes, min_run)
    return codes, list(zip(s.tolist(), n.tolist()))


def acgt_stream(seq_bytes: bytes, keep_gaps: bool = False) -> np.ndarray:
    """One concatenated ACGT code stream.

    keep_gaps=False (default, original behaviour): non-ACGT characters are
        deleted, so bases on either side of a gap become adjacent. Kept for
        backward compatibility only -- do not use it for pattern matching.
    keep_gaps=True: each maximal block of non-ACGT characters becomes a single
        GAP code (4), including blocks at the start or end of the record.
        Callers must split on GAP or mask it; projections.ry/ws/mk would
        otherwise turn it into an ordinary 0 bit.

    Neither form preserves genome coordinates. Prefer encode() + iter_segments().
    """
    codes = encode(seq_bytes)
    valid = codes < 4
    if not keep_gaps:
        return codes[valid]
    # keep valid bases and the first position of each invalid block
    first_bad = ~valid & np.concatenate(([True], valid[:-1]))
    out = codes[valid | first_bad].copy()
    out[out == INVALID] = GAP
    return out
