var Jt=Object.defineProperty;var $t=(r,e,t)=>e in r?Jt(r,e,{enumerable:!0,configurable:!0,writable:!0,value:t}):r[e]=t;var M=(r,e,t)=>$t(r,typeof e!="symbol"?e+"":e,t);class ve{constructor(e){M(this,"size");M(this,"isPowerOfTwo");M(this,"_real");M(this,"_imag");M(this,"_scratch",null);M(this,"_rev",null);M(this,"_m",0);M(this,"_internalFFT",null);M(this,"_chirpReal",null);M(this,"_chirpImag",null);M(this,"_bReal",null);M(this,"_bImag",null);M(this,"_hanning",null);M(this,"_windowSumSq",0);this.size=e,this.isPowerOfTwo=(e&e-1)===0&&e>0,this._real=new Float32Array(e),this._imag=new Float32Array(e),this.isPowerOfTwo?this.initRadix2():this.initBluestein()}initRadix2(){const e=this.size,t=Math.log2(e);this._rev=new Uint32Array(e);for(let n=0;n<e;n++){let i=0,s=n;for(let a=0;a<t;a++)i=i<<1|s&1,s>>>=1;this._rev[n]=i}}initBluestein(){const e=this.size;this._m=Math.pow(2,Math.ceil(Math.log2(2*e-1))),this._internalFFT=new ve(this._m),this._chirpReal=new Float32Array(e),this._chirpImag=new Float32Array(e);for(let i=0;i<e;i++){const s=-Math.PI*(i*i)/e;this._chirpReal[i]=Math.cos(s),this._chirpImag[i]=Math.sin(s)}const t=new Float32Array(this._m),n=new Float32Array(this._m);for(let i=0;i<e;i++)t[i]=this._chirpReal[i],n[i]=-this._chirpImag[i];for(let i=1;i<e;i++)t[this._m-i]=t[i],n[this._m-i]=n[i];this._internalFFT.transform(t,n),this._bReal=new Float32Array(this._internalFFT._real),this._bImag=new Float32Array(this._internalFFT._imag)}initHanning(){if(this._hanning)return;const e=this.size;this._hanning=new Float32Array(e);let t=0;for(let n=0;n<e;n++){const i=.5*(1-Math.cos(2*Math.PI*n/(e-1)));this._hanning[n]=i,t+=i*i}this._windowSumSq=t}transform(e,t){this.isPowerOfTwo?this.transformRadix2(e,t):this.transformBluestein(e,t)}transformRadix2(e,t){const n=this.size,i=this._rev,s=this._real,a=this._imag;if(e===s)for(let o=0;o<n;o++){const c=i[o];if(o<c){const u=s[o],d=a[o];s[o]=s[c],a[o]=a[c],s[c]=u,a[c]=d}}else for(let o=0;o<n;o++){const c=i[o];s[o]=e[c],a[o]=t?t[c]:0}for(let o=2;o<=n;o*=2){const c=o/2,u=-2*Math.PI/o,d=Math.cos(u),l=Math.sin(u);for(let f=0;f<n;f+=o){let h=1,p=0;for(let y=0;y<c;y++){const w=f+y,b=f+y+c,C=h*s[b]-p*a[b],R=h*a[b]+p*s[b],S=s[w],k=a[w];s[w]=S+C,a[w]=k+R,s[b]=S-C,a[b]=k-R;const m=h*d-p*l,O=h*l+p*d;h=m,p=O}}}}transformBluestein(e,t){const n=this.size,i=this._m,s=this._internalFFT,a=s._real,o=s._imag;a.fill(0),o.fill(0);for(let l=0;l<n;l++){const f=e[l],h=t?t[l]:0,p=this._chirpReal[l],y=this._chirpImag[l];a[l]=f*p-h*y,o[l]=f*y+h*p}s.transformRadix2(a,o);for(let l=0;l<i;l++){const f=s._real[l],h=s._imag[l],p=this._bReal[l],y=this._bImag[l];s._real[l]=f*p-h*y,s._imag[l]=f*y+h*p}const c=s._real,u=s._imag;for(let l=0;l<i;l++)u[l]=-u[l];s.transformRadix2(c,u);const d=1/i;for(let l=0;l<n;l++){const f=s._real[l]*d,h=-s._imag[l]*d,p=this._chirpReal[l],y=this._chirpImag[l];this._real[l]=f*p-h*y,this._imag[l]=f*y+h*p}}calculateSpectrum(e,t,n=!1){const i=this.size;let s=0;for(let d=0;d<i;d++)s+=e[d];const a=s/i;this._scratch||(this._scratch=new Float32Array(i));const o=this._scratch;if(n){this.initHanning();const d=this._hanning;for(let l=0;l<i;l++)o[l]=(e[l]-a)*d[l]}else for(let d=0;d<i;d++)o[d]=e[d]-a;this.transform(o);const c=t.length;let u=1/i;n&&this._windowSumSq>0&&(u=1/this._windowSumSq);for(let d=0;d<c;d++){const l=this._real[d],f=this._imag[d];t[d]+=(l*l+f*f)*u}}calculateSpectrumWindow(e,t,n,i=!1){const s=this.size;let a=0;for(let l=0;l<s;l++)a+=e[t+l];const o=a/s;this._scratch||(this._scratch=new Float32Array(s));const c=this._scratch;if(i){this.initHanning();const l=this._hanning;for(let f=0;f<s;f++)c[f]=(e[t+f]-o)*l[f]}else for(let l=0;l<s;l++)c[l]=e[t+l]-o;this.transform(c);const u=n.length;let d=1/s;i&&this._windowSumSq>0&&(d=1/this._windowSumSq);for(let l=0;l<u;l++){const f=this._real[l],h=this._imag[l];n[l]+=(f*f+h*h)*d}}}const Qt={"Sony ILCE-7RM5":"0.82 -0.2976 -0.0719 -0.4296 1.2053 0.2532 -0.0429 0.1282 0.5774"};let Ne=null;async function Zt(r){return Ne||(Ne=(async()=>{if(typeof window.loadPyodide!="function")throw new Error("Pyodide missing: window.loadPyodide not found.");const e=await window.loadPyodide();return await e.loadPackage("numpy"),e})()),Ne}var en=`#!/usr/bin/env python3
"""Container parser for the Sony ARW6/LLVC3 raw strip.

This stops at packet records and control fields. Coefficient entropy decoding
lives in llvc3_entropy.py; this file is mostly the thing I dump to JSON and diff
against Imaging Edge traces when the packet framing looks suspicious.
"""

from __future__ import annotations

import argparse
import json
import struct
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any


RAW_STREAM_OFFSET = 0x200


def u16be(buf: bytes, off: int) -> int:
    return struct.unpack_from(">H", buf, off)[0]


def u16le(buf: bytes, off: int) -> int:
    return struct.unpack_from("<H", buf, off)[0]


def u32le(buf: bytes, off: int) -> int:
    return struct.unpack_from("<I", buf, off)[0]


def u32be(buf: bytes, off: int) -> int:
    return struct.unpack_from(">I", buf, off)[0]


class BitReader:
    """Small MSB-first reader for packet headers."""

    def __init__(self, data: bytes) -> None:
        self.data = data
        self.pos = 0

    def read(self, nbits: int) -> int:
        if nbits < 0:
            raise ValueError("negative bit count")
        out = 0
        for _ in range(nbits):
            if self.pos >= len(self.data) * 8:
                raise EOFError("packet bitstream exhausted")
            byte = self.data[self.pos >> 3]
            bit = (byte >> (7 - (self.pos & 7))) & 1
            out = (out << 1) | bit
            self.pos += 1
        return out


@dataclass
class TiffRawInfo:
    ifd_offset: int
    width: int
    height: int
    bits_per_sample: int
    compression: int
    photometric: int
    strip_offset: int
    strip_byte_count: int
    cfa_pattern: list[int] | None
    black_level_tag_0x7310: list[int] | None
    white_level: int | None
    default_crop_origin: list[int] | None
    default_crop_size: list[int] | None


@dataclass
class LlvcHeader:
    magic: str
    sequence_or_version: int
    coded_width: int
    coded_half_height: int
    logical_height: int
    decoded_bits: int
    component_count: int
    mode: int
    flags_low10: int


@dataclass
class DirectoryEntry:
    group: int
    index: int
    start: int
    length: int
    mode_in_decoder_trace: int | None = None


@dataclass
class PacketRecord:
    index: int
    byte_length: int
    selectors: list[int]
    payload_offset: int


@dataclass
class PacketHeader:
    group: int
    index: int
    stream_offset: int
    directory_length: int
    control_count: int
    extra_count: int
    tag4: int
    reserved4: int
    type2: int
    control_words: int
    block_count: int
    width_marker: int
    reserved8: int
    skipped_u8: list[int]
    header_bits: int
    control_bytes: int
    total_bytes: int
    payload_bytes_from_records: int
    validation: dict[str, bool]
    records: list[PacketRecord]
    first_records: list[PacketRecord]
    last_records: list[PacketRecord]


def parse_tiff_value(data: bytes, bo: str, typ: int, cnt: int, raw: bytes) -> Any:
    fmt_by_type = {1: "B", 3: "H", 4: "I", 8: "h", 9: "i"}
    if typ in fmt_by_type:
        fmt = bo + fmt_by_type[typ] * cnt
        size = struct.calcsize(fmt)
        val = list(struct.unpack(fmt, raw[:size]))
        return val[0] if len(val) == 1 else val
    if typ in (2, 7):
        return list(raw[:cnt])
    return raw.hex(" ")


def iter_tiff_ifds(data: bytes) -> list[tuple[int, dict[int, Any]]]:
    if data[:2] == b"II":
        bo = "<"
    elif data[:2] == b"MM":
        bo = ">"
    else:
        raise ValueError("input is not a TIFF/ARW file")

    def u16(off: int) -> int:
        return struct.unpack_from(bo + "H", data, off)[0]

    def u32(off: int) -> int:
        return struct.unpack_from(bo + "I", data, off)[0]

    type_sizes = {1: 1, 2: 1, 3: 2, 4: 4, 5: 8, 7: 1, 8: 2, 9: 4, 10: 8}
    stack = [u32(4)]
    seen: set[int] = set()
    out: list[tuple[int, dict[int, Any]]] = []
    while stack:
        off = stack.pop()
        if off in seen or off <= 0 or off >= len(data):
            continue
        seen.add(off)
        n = u16(off)
        tags: dict[int, Any] = {}
        for i in range(n):
            ent = off + 2 + i * 12
            tag = u16(ent)
            typ = u16(ent + 2)
            cnt = u32(ent + 4)
            value_area = ent + 8
            size = type_sizes.get(typ, 1) * cnt
            raw = data[value_area : value_area + 4] if size <= 4 else data[u32(value_area) : u32(value_area) + size]
            tags[tag] = parse_tiff_value(data, bo, typ, cnt, raw)
            if tag == 0x014A:
                offsets = tags[tag] if isinstance(tags[tag], list) else [tags[tag]]
                stack.extend(int(x) for x in offsets)
        next_ifd = u32(off + 2 + n * 12)
        if next_ifd:
            stack.append(next_ifd)
        out.append((off, tags))
    return out


def scalar_or_first(value: Any) -> int:
    if isinstance(value, list):
        return int(value[0])
    return int(value)


def find_raw_subifd(path: Path) -> tuple[TiffRawInfo, bytes]:
    data = path.read_bytes()
    for off, tags in iter_tiff_ifds(data):
        if tags.get(0x0103) == 32766 and tags.get(0x0106) == 32803:
            strip_offset = scalar_or_first(tags[0x0111])
            strip_len = scalar_or_first(tags[0x0117])
            info = TiffRawInfo(
                ifd_offset=off,
                width=int(tags[0x0100]),
                height=int(tags[0x0101]),
                bits_per_sample=int(tags[0x0102]),
                compression=int(tags[0x0103]),
                photometric=int(tags[0x0106]),
                strip_offset=strip_offset,
                strip_byte_count=strip_len,
                cfa_pattern=list(tags.get(0x828E, [])) or None,
                black_level_tag_0x7310=list(tags.get(0x7310, [])) if isinstance(tags.get(0x7310), list) else None,
                white_level=int(tags[0xC61D]) if 0xC61D in tags else None,
                default_crop_origin=list(tags.get(0xC61F, [])) if isinstance(tags.get(0xC61F), list) else None,
                default_crop_size=list(tags.get(0xC620, [])) if isinstance(tags.get(0xC620), list) else None,
            )
            return info, data[strip_offset : strip_offset + strip_len]
    raise ValueError("no ARW6 LLVC raw SubIFD found")


def parse_llvc_header(stream: bytes) -> LlvcHeader:
    word_c = u16be(stream, 0x0C)
    word_e = u16be(stream, 0x0E)
    return LlvcHeader(
        magic=stream[:4].decode("ascii", "replace"),
        sequence_or_version=u32le(stream, 0x04),
        coded_width=u16be(stream, 0x08),
        coded_half_height=u16be(stream, 0x0A),
        logical_height=u16be(stream, 0x0A) * 2,
        decoded_bits=(word_c >> 4) & 0x3F,
        component_count=word_e >> 13,
        mode=(word_e >> 10) & 0x03,
        flags_low10=word_e & 0x03FF,
    )


def initial_group_lengths(stream: bytes) -> list[int]:
    return [((u32be(stream, 0x10 + off) >> 4) & 0x0FFFFFF0) for off in (0, 3, 6, 9, 12)]


def parse_directory(stream: bytes) -> tuple[int, list[DirectoryEntry], list[dict[str, Any]]]:
    group_lengths = initial_group_lengths(stream)
    pos = 0x30
    base = 0
    entries: list[DirectoryEntry] = []
    groups: list[dict[str, Any]] = []
    for group, group_len in enumerate(group_lengths):
        n_entries = stream[pos] & 0x0F
        values = [(u32be(stream, pos + off) & 0x00FFFFFF) << 4 for off in (0, 3, 6, 9, 12)]
        consumed = 0x10
        if n_entries >= 5:
            values.extend((u32be(stream, pos + 0x10 + off) & 0x00FFFFFF) << 4 for off in (0, 3, 6, 9))
            consumed = 0x20
        local = 0
        for index in range(n_entries):
            entries.append(DirectoryEntry(group=group, index=index, start=0x80 + base + local, length=values[index]))
            local += values[index]
        groups.append(
            {
                "group": group,
                "declared_length": group_len,
                "entry_count": n_entries,
                "entry_lengths": values[:n_entries],
                "entry_sum": local,
                "directory_offset": pos,
                "directory_bytes": consumed,
                "sum_matches_declared": local == group_len,
            }
        )
        base += group_len
        pos += consumed
    return pos, entries, groups


def parse_packet(stream: bytes, entry: DirectoryEntry) -> PacketHeader:
    packet = stream[entry.start : entry.start + entry.length]
    br = BitReader(packet)
    control_count = br.read(16)
    extra_count = br.read(24)
    tag4 = br.read(4)
    reserved4 = br.read(4)
    type2 = br.read(2)
    control_words = br.read(6)
    block_count = br.read(16)
    width_marker = br.read(8)
    reserved8 = br.read(8)
    skipped = [br.read(8) for _ in range(5)]
    control_bytes = (control_count + 1) << 4
    total_bytes = (control_count + 1 + extra_count) << 4
    records: list[PacketRecord] = []
    cursor = control_bytes
    for i in range(block_count):
        byte_len = br.read(16)
        selectors = [br.read(4) for _ in range(type2)]
        records.append(PacketRecord(i, byte_len, selectors, cursor))
        cursor += byte_len
    # Native code consumes the 6-bit field while walking the header, then the
    # validation at 0x1aa717 compares the original 16-bit control_count with:
    # ceil((type2 + 4) * block_count * 4 / 128).
    formula = ((type2 + 4) * block_count * 4 + 0x7F) >> 7
    validation = {
        "directory_length_matches_total": entry.length == total_bytes,
        "tag4_is_4": tag4 == 4,
        "reserved4_is_0": reserved4 == 0,
        "type2_is_1_or_3": type2 in (1, 3),
        "reserved6_is_0": control_words == 0,
        "control_count_formula": control_count == formula,
        "block_count_le_300": block_count <= 300,
        "width_marker_is_0x10": width_marker == 0x10,
        "reserved8_is_0": reserved8 == 0,
        "skipped_bytes_are_0": all(x == 0 for x in skipped),
        "record_payload_within_total": cursor <= total_bytes,
    }
    return PacketHeader(
        group=entry.group,
        index=entry.index,
        stream_offset=entry.start,
        directory_length=entry.length,
        control_count=control_count,
        extra_count=extra_count,
        tag4=tag4,
        reserved4=reserved4,
        type2=type2,
        control_words=control_words,
        block_count=block_count,
        width_marker=width_marker,
        reserved8=reserved8,
        skipped_u8=skipped,
        header_bits=br.pos,
        control_bytes=control_bytes,
        total_bytes=total_bytes,
        payload_bytes_from_records=cursor,
        validation=validation,
        records=records,
        first_records=records[:8],
        last_records=records[-4:],
    )


def derive_metrics(raw_info: TiffRawInfo, strip_len: int, stream_len: int) -> dict[str, Any]:
    samples = raw_info.width * raw_info.height
    return {
        "samples": samples,
        "decoded_u16_bytes": samples * 2,
        "decoded_14bit_packed_bytes": samples * 14 / 8,
        "strip_bytes_per_sample": strip_len / samples,
        "strip_bits_per_sample": strip_len * 8 / samples,
        "stream_bits_per_sample_excluding_0x200_preamble": stream_len * 8 / samples,
        "compression_ratio_vs_u16": (samples * 2) / strip_len,
        "compression_ratio_vs_14bit_packed": (samples * 14 / 8) / strip_len,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("path", nargs="?", default="DSC00089.ARW")
    ap.add_argument("--out", default="out/reverse/llvc3_bitstream_probe.json")
    ns = ap.parse_args()

    raw_info, strip = find_raw_subifd(Path(ns.path))
    stream = strip[RAW_STREAM_OFFSET:]
    header = parse_llvc_header(stream)
    consumed, entries, groups = parse_directory(stream)
    packets = [parse_packet(stream, entry) for entry in entries]
    all_valid = all(all(p.validation.values()) for p in packets)
    type_counts: dict[str, int] = {}
    for p in packets:
        type_counts[str(p.type2)] = type_counts.get(str(p.type2), 0) + 1
    result = {
        "input": str(ns.path),
        "raw_subifd": asdict(raw_info),
        "llvc_header": asdict(header),
        "stream_offset_inside_raw_strip": RAW_STREAM_OFFSET,
        "strip_preamble_first_32": strip[:32].hex(" "),
        "directory": {
            "consumed_bytes_after_stream_header": consumed - 0x10,
            "packet_base_offset": 0x80,
            "groups": groups,
            "entries": [asdict(e) for e in entries],
        },
        "packets": [asdict(p) for p in packets],
        "summary": {
            "packet_count": len(packets),
            "packet_type_counts": type_counts,
            "all_packet_validations_pass": all_valid,
            "total_packet_bytes": sum(p.directory_length for p in packets),
            "payload_bytes_from_records": sum(p.payload_bytes_from_records for p in packets),
            "metrics": derive_metrics(raw_info, len(strip), len(stream)),
        },
    }
    out = Path(ns.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result["summary"], indent=2))


if __name__ == "__main__":
    main()
`,tn=`#!/usr/bin/env python3
"""Entropy-side notes for the Sony ARW6/LLVC3 stream.

The bit reader and 4-lane coefficient paths here came straight out of Imaging
Edge traces. I kept the code narrow on purpose: first replay a row, then a
packet, then let the higher-level decoder stitch the pieces together.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from llvc3_bitstream_probe import RAW_STREAM_OFFSET, find_raw_subifd


@dataclass
class NativeBitReader:
    """MSB-first 64-bit reader with the same state layout Sony uses.

    The offsets are worth keeping close by:

      +0x00 current big-endian 64-bit word
      +0x08 pointer to the last loaded word
      +0x10 bit position / remaining bits in current word
      +0x14 64-bit words left to load
      +0x1c status

    The pointer is initialized to eight bytes before the record payload. On
    underflow the decoder advances by eight bytes and loads a big-endian word.
    """

    data: bytes
    ptr: int = -8
    cur: int = 0
    bit: int = 0
    words_left: int = 0
    status: int = 0

    @classmethod
    def for_record(cls, packet: bytes, payload_offset: int, byte_length: int) -> "NativeBitReader":
        words = (((byte_length + 7) // 8) + 1) & ~1
        return cls(data=packet, ptr=payload_offset - 8, words_left=words)

    def _load_next_word(self) -> int:
        self.words_left -= 1
        if self.words_left < 0:
            self.status = 2
            return 0
        self.ptr += 8
        chunk = self.data[self.ptr : self.ptr + 8]
        if len(chunk) < 8:
            chunk = chunk + b"\\x00" * (8 - len(chunk))
        self.cur = int.from_bytes(chunk, "big")
        self.bit = 64
        return self.cur

    def read_bits(self, nbits: int) -> int:
        if nbits == 0:
            return 0
        if nbits < 0:
            raise ValueError("negative bit count")
        out = 0
        remaining = nbits
        while remaining > 0:
            if self.bit <= 0:
                self._load_next_word()
                if self.status:
                    return out << remaining
            take = min(remaining, self.bit)
            self.bit -= take
            out = (out << take) | ((self.cur >> self.bit) & ((1 << take) - 1))
            remaining -= take
        return out

    def skip_zero_words_to_one(self, initial_bit_count: int) -> int:
        """Count leading zeros until the next one, in the 0x1a9080 style."""

        total = 0
        while True:
            if self.bit <= 0:
                self._load_next_word()
                if self.status:
                    return total
            if self.cur == 0:
                total += self.bit
                self.bit = 0
                continue
            # Count zeros among currently available MSB-side bits.
            while self.bit > 0:
                self.bit -= 1
                if (self.cur >> self.bit) & 1:
                    return total
                total += 1
                if total >= initial_bit_count and initial_bit_count > 0:
                    return total

    def read_unary_zeros_plus_one(self) -> int:
        zeros = 0
        while not self.status:
            bit = self.read_bits(1)
            if self.status:
                break
            if bit:
                return zeros + 1
            zeros += 1
        return zeros + 1


def split_lanes(packed: int, width: int) -> list[int]:
    mask = (1 << width) - 1
    return [
        (packed >> (3 * width)) & mask,
        (packed >> (2 * width)) & mask,
        (packed >> width) & mask,
        packed & mask,
    ]


def magnitude4(br: NativeBitReader, width: int, shift: int) -> list[int]:
    """Magnitude expansion from 0x1a8b00."""

    if width <= 0:
        raw = [0, 0, 0, 0]
    else:
        raw = split_lanes(br.read_bits(width * 4), width)
    if shift <= 0:
        return raw
    out: list[int] = []
    for x in raw:
        if x <= 0:
            out.append(0)
        else:
            # Visible SIMD shape: ((2*x + 1) << (shift - 1)) - (x & 1).
            out.append(((2 * x + 1) << (shift - 1)) - (x & 1))
    return out


def apply_sign4(br: NativeBitReader, coeffs: Iterable[int]) -> list[int]:
    """Apply one sign bit to each positive lane; see 0x1a8dd0."""

    out: list[int] = []
    for x in coeffs:
        if x > 0:
            bit = br.read_bits(1)
            out.append(x - 2 * x * bit)
        else:
            out.append(x)
    return out


def update_width(br: NativeBitReader, width: int) -> int:
    """Adaptive width update shared by 0x1a9080 and 0x1ac060.

    Prefix structure from traces and branch shape:

      0      keep width
      10 U   increase by unary U, where U is zeros+1 terminated by a one
      11 U   decrease by unary U, clipped at zero
    """

    if br.read_bits(1) == 0 or br.status:
        return width
    if br.read_bits(1) == 0 or br.status:
        return width + br.read_unary_zeros_plus_one()

    # Decrement path: if the run would cross zero, native code consumes only
    # width-1 zeros and leaves the terminating one for the next state; that is
    # the 0x1ac4ff/0x1ac50b saturation branch.
    for zeros in range(max(0, width - 1)):
        if br.read_bits(1):
            return width - (zeros + 1)
    return 0


def read_initial_width(br: NativeBitReader) -> int:
    """Read the row's first adaptive width with the 0x1a9080 zero state."""

    return update_width(br, 0)


def read_zero_run(br: NativeBitReader, remaining: int) -> int:
    """Read the zero-group run used when the adaptive width is zero."""

    if remaining <= 1:
        return remaining
    max_prefix = (remaining - 1).bit_length()
    zeros = 0
    while zeros < max_prefix:
        if br.read_bits(1):
            break
        if br.status:
            return remaining
        zeros += 1
    else:
        return remaining
    base = 1 << zeros
    if base >= remaining:
        return remaining
    extra = br.read_bits(zeros) if zeros else 0
    run = base + extra
    return min(run, remaining)


def packet_row_multiplier(group: int) -> int:
    """Number of output rows produced by one packet record for this scale."""

    if group < 0:
        raise ValueError("negative packet group")
    if group == 4:
        return 8
    if group == 0:
        return 1
    return 1 << (group - 1)


def read_width_after_zero_run(br: NativeBitReader) -> int:
    """Positive width code after a zero run: leading zeros plus the one bit."""

    return br.read_unary_zeros_plus_one()


def alt4(br: NativeBitReader, width: int, shift: int, has_next: bool = True) -> tuple[list[int], int]:
    """Current reconstruction of the alternate 0x1ac060 4-lane path.

    It shares the magnitude transform with 0x1a8b00, then signs only positive
    lanes. The width predictor agrees with the small-width traces I have, but
    wider rows still deserve spot checks.
    """

    coeffs = magnitude4(br, width, shift)
    next_width = width
    if has_next:
        next_width = update_width(br, width)
    return apply_sign4(br, coeffs), next_width


def decode_record_component(
    br: NativeBitReader, groups: int, shift: int = 0
) -> tuple[list[int], list[int], int]:
    """Decode one component from the current bitreader position.

    \`groups\` is the number of 4-lane coefficient groups. The packet selector
    nibble is passed as \`shift\`; the initial adaptive width is read from the
    record payload itself by the 0x1a9080 prefix reader.
    """

    width = read_initial_width(br)
    initial_width = width
    coeffs: list[int] = []
    widths: list[int] = []
    gi = 0
    while gi < groups:
        if br.status or width > 0x13:
            br.status = br.status or 1
            coeffs.extend([0, 0, 0, 0] * (groups - gi))
            widths.extend([width] * (groups - gi))
            break
        if width == 0:
            run = read_zero_run(br, groups - gi)
            coeffs.extend([0, 0, 0, 0] * run)
            widths.extend([0] * run)
            gi += run
            if gi >= groups:
                break
            width = read_width_after_zero_run(br)
            continue
        vals, width = alt4(br, width, shift=shift, has_next=gi + 1 < groups)
        coeffs.extend(vals)
        widths.append(width)
        gi += 1
    return coeffs, widths, initial_width


def load_packet(arw: Path, group: int, index: int) -> tuple[bytes, dict]:
    from llvc3_bitstream_probe import parse_directory, parse_llvc_header, parse_packet

    raw_info, strip = find_raw_subifd(arw)
    stream = strip[RAW_STREAM_OFFSET:]
    parse_llvc_header(stream)
    _consumed, entries, _groups = parse_directory(stream)
    entry = next(e for e in entries if e.group == group and e.index == index)
    packet_info = parse_packet(stream, entry)
    packet = stream[entry.start : entry.start + entry.length]
    info = json.loads(json.dumps(packet_info, default=lambda o: o.__dict__))
    info["raw_width"] = raw_info.width
    info["raw_height"] = raw_info.height
    return packet, info


def replay_row(arw: Path, group: int, index: int, row: int, groups: int = 16) -> dict[str, object]:
    packet, info = load_packet(arw, group, index)
    rec = info["records"][row]
    coeffs, widths, br, initial_width = decode_type1_row(packet, rec, groups)
    return {
        "packet": {"group": group, "index": index},
        "row": row,
        "record": rec,
        "selector": rec["selectors"][0] if rec["selectors"] else 0,
        "initial_width": initial_width,
        "coeffs": coeffs,
        "next_widths": widths,
        "bitreader": {"ptr": br.ptr, "bit": br.bit, "words_left": br.words_left, "status": br.status},
    }


def decode_type1_row(packet: bytes, rec: dict, groups: int) -> tuple[list[int], list[int], NativeBitReader, int]:
    """Decode one type-1 record into 4-lane signed coefficients."""

    br = NativeBitReader.for_record(packet, rec["payload_offset"], rec["byte_length"])
    shift = rec["selectors"][0] if rec["selectors"] else 0
    coeffs, widths, initial_width = decode_record_component(br, groups, shift)
    return coeffs, widths, br, initial_width


def decode_record_components(
    packet: bytes, rec: dict, groups: int, components: int, row_multiplier: int = 1
) -> tuple[list[list[int]], list[dict[str, int]]]:
    """Decode all components stored in one packet record.

    Type-3 records share one payload bitreader across three component streams.
    The native block decoder calls the same entropy routine three times with
    selector nibbles from the control record, so the reader stays live between
    components.
    """

    br = NativeBitReader.for_record(packet, rec["payload_offset"], rec["byte_length"])
    rows: list[list[int]] = []
    states: list[dict[str, int]] = []
    selectors = rec["selectors"] or []
    for ci in range(components):
        shift = selectors[ci] if ci < len(selectors) else 0
        for _ in range(row_multiplier):
            if rec["byte_length"] <= 0:
                coeffs = [0] * (groups * 4)
            else:
                coeffs, _widths, _initial = decode_record_component(br, groups, shift)
            rows.append(coeffs)
            states.append({"ptr": br.ptr, "bit": br.bit, "words_left": br.words_left, "status": br.status})
    return rows, states


def infer_packet_width(group: int, packet_type: int, raw_width: int = 7040) -> int:
    """Entropy row width in coefficients for this ARW6 sample's packet groups.

    Type-3 packets carry wavelet detail subbands. Their entropy rows are
    half-width relative to the reconstructed scale reported by the block
    object; the horizontal synthesis stage doubles that later.
    """

    if packet_type == 1:
        return raw_width // 2 if group == 4 else raw_width // 16
    if packet_type == 3:
        if 1 <= group <= 3:
            return raw_width // (1 << (5 - group))
    raise ValueError(f"width inference not yet known for group {group}, type {packet_type}")


def decode_packet_components(arw: Path, group: int, index: int, out_prefix: Path | None = None) -> dict[str, object]:
    """Decode a type-1 or type-3 packet into one or three int32 component arrays."""

    packet, info = load_packet(arw, group, index)
    packet_type = info["type2"]
    components = 1 if packet_type == 1 else 3
    width = infer_packet_width(group, packet_type, int(info.get("raw_width", 7040)))
    row_multiplier = packet_row_multiplier(group)
    groups_per_row = (width + 3) // 4
    planes: list[list[list[int]]] = [[] for _ in range(components)]
    final_states: list[list[dict[str, int]]] = []
    for rec in info["records"]:
        comp_rows, states = decode_record_components(packet, rec, groups_per_row, components, row_multiplier)
        for ci in range(components):
            start = ci * row_multiplier
            end = start + row_multiplier
            planes[ci].extend(row[:width] for row in comp_rows[start:end])
        final_states.append(states)
    outs: list[str] = []
    if out_prefix:
        import numpy as np

        out_prefix.parent.mkdir(parents=True, exist_ok=True)
        for ci, rows in enumerate(planes):
            arr = np.asarray(rows, dtype=np.int32)
            path = out_prefix.with_name(f"{out_prefix.name}_c{ci}.bin")
            arr.tofile(path)
            outs.append(str(path))
    return {
        "packet": {"group": group, "index": index, "type": packet_type},
        "shape": [len(planes[0]), width],
        "row_multiplier": row_multiplier,
        "components": components,
        "first_row_first16": [plane[0][:16] for plane in planes],
        "last_nonempty_row_first16": [
            next((plane[i][:16] for i in range(len(plane) - 1, -1, -1) if any(plane[i])), []) for plane in planes
        ],
        "final_states_tail": final_states[-5:],
        "outs": outs,
    }


def decode_packet_arrays(arw: Path, group: int, index: int) -> tuple[list["object"], dict[str, object]]:
    """Decode a packet and return its int32 component arrays in memory."""

    import numpy as np

    packet, info = load_packet(arw, group, index)
    packet_type = info["type2"]
    components = 1 if packet_type == 1 else 3
    width = infer_packet_width(group, packet_type, int(info.get("raw_width", 7040)))
    row_multiplier = packet_row_multiplier(group)
    groups_per_row = (width + 3) // 4
    planes: list[list[list[int]]] = [[] for _ in range(components)]
    final_states: list[list[dict[str, int]]] = []
    for rec in info["records"]:
        comp_rows, states = decode_record_components(packet, rec, groups_per_row, components, row_multiplier)
        for ci in range(components):
            start = ci * row_multiplier
            end = start + row_multiplier
            planes[ci].extend(row[:width] for row in comp_rows[start:end])
        final_states.append(states)

    arrays = [np.asarray(rows, dtype=np.int32) for rows in planes]
    meta = {
        "packet": {"group": group, "index": index, "type": packet_type},
        "shape": [int(arrays[0].shape[0]), int(arrays[0].shape[1])],
        "row_multiplier": row_multiplier,
        "components": components,
        "final_state": final_states[-1] if final_states else [],
    }
    return arrays, meta


def decode_type1_packet(arw: Path, group: int, index: int, out: Path | None = None) -> dict[str, object]:
    """Decode an independently parsed type-1 packet into int32 coefficient rows."""

    packet, info = load_packet(arw, group, index)
    if info["type2"] != 1:
        raise ValueError(f"packet g{group}i{index} is type {info['type2']}, not type 1")
    width = infer_packet_width(group, info["type2"], int(info.get("raw_width", 7040)))
    groups_per_row = (width + 3) // 4
    rows: list[list[int]] = []
    final_states: list[dict[str, int]] = []
    for rec in info["records"]:
        if rec["byte_length"] <= 0:
            coeffs = [0] * (groups_per_row * 4)
            br = NativeBitReader.for_record(packet, rec["payload_offset"], rec["byte_length"])
        else:
            coeffs, _widths, br, _initial = decode_type1_row(packet, rec, groups_per_row)
        rows.append(coeffs[:width])
        final_states.append({"ptr": br.ptr, "bit": br.bit, "words_left": br.words_left, "status": br.status})
    if out:
        import numpy as np

        arr = np.asarray(rows, dtype=np.int32)
        out.parent.mkdir(parents=True, exist_ok=True)
        arr.tofile(out)
    return {
        "packet": {"group": group, "index": index, "type": info["type2"]},
        "shape": [len(rows), width],
        "first_row_first16": rows[0][:16],
        "last_nonempty_row_first16": next((rows[i][:16] for i in range(len(rows) - 1, -1, -1) if any(rows[i])), []),
        "final_states_tail": final_states[-5:],
        "out": str(out) if out else "",
    }


def integrate_type1_coefficients(coeffs: "object", dc_offset: int) -> "object":
    """Apply the row postprocess after the entropy call at 0x1a7fb3.

    For each row:

      acc0 = int16(coeff[0]) * 2
      out[0] = acc0 >> 1
      acc_i = acc_{i-1} + int16(coeff[i]) * 2
      out[i] = acc_i >> 1

    Applying this to group0/index0 and adding 2048 reproduces the native v4 c0
    lowpass plane.
    """

    import numpy as np

    c = np.asarray(coeffs, dtype=np.int32)
    out = np.empty_like(c, dtype=np.int32)
    signed = c.astype(np.int16).astype(np.int32)
    acc = signed[:, 0] * 2
    out[:, 0] = acc >> 1
    for x in range(1, c.shape[1]):
        acc = acc + signed[:, x] * 2
        out[:, x] = acc >> 1
    return (out + dc_offset).astype(np.int32)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("arw", nargs="?", default="DSC00089.ARW")
    ap.add_argument("--group", type=int, default=0)
    ap.add_argument("--index", type=int, default=2)
    ap.add_argument("--row", type=int, default=280)
    ap.add_argument("--groups", type=int, default=16)
    ap.add_argument("--packet", action="store_true", help="decode the whole type-1 packet instead of one row")
    ap.add_argument("--components", action="store_true", help="decode a type-1/type-3 packet into component files")
    ap.add_argument("--out", default="")
    ns = ap.parse_args()
    if ns.components:
        result = decode_packet_components(Path(ns.arw), ns.group, ns.index, Path(ns.out) if ns.out else None)
    elif ns.packet:
        result = decode_type1_packet(Path(ns.arw), ns.group, ns.index, Path(ns.out) if ns.out else None)
    else:
        result = replay_row(Path(ns.arw), ns.group, ns.index, ns.row, ns.groups)
    text = json.dumps(result, indent=2)
    if ns.out and not ns.packet and not ns.components:
        Path(ns.out).write_text(text, encoding="utf-8")
    print(text)


if __name__ == "__main__":
    main()
`,nn=`#!/usr/bin/env python3
"""Integer helpers from the LLVC3 reversing notes.

Closer to a lab notebook than a polished codec module. The small functions
below are named around the traces they came from, then reused by the pure
decoder once a stage lines up.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


INTERNAL_BIAS = 2048


def clamp_sample(x: np.ndarray, bits: int) -> np.ndarray:
    return np.clip(x, 0, (1 << bits) - 1)


def merge_average_detail(lo0: np.ndarray, lo1: np.ndarray, detail: np.ndarray, bits: int = 16) -> np.ndarray:
    """Inverse merge from Edit.exe RVA 0x1a9a40.

    The scalar path is simple enough to leave here as a breadcrumb:

        avg = (lo0 + lo1) >> 1
        sample = avg + 2 * signed_detail
        sample = clamp(sample, 0, (1 << bits) - 1)

    Native output buffers are 16-bit, even when the TIFF tags describe 14-bit
    sensor samples, so the helper returns uint16 too.
    """

    x = ((lo0.astype(np.int32) + lo1.astype(np.int32)) >> 1) + 2 * detail.astype(np.int32)
    return clamp_sample(x, bits).astype(np.uint16)


def add_detail(low: np.ndarray, detail: np.ndarray, bits: int = 16) -> np.ndarray:
    """Lowpass plus signed detail, from the 0x1aafd0 helper family."""

    x = low.astype(np.int32) + detail.astype(np.int32)
    return clamp_sample(x, bits).astype(np.uint16)


def add_double_detail(low: np.ndarray, detail: np.ndarray, bits: int = 16) -> np.ndarray:
    """Lowpass plus 2*detail, seen around the 0x1ab2b0 helpers."""

    x = low.astype(np.int32) + 2 * detail.astype(np.int32)
    return clamp_sample(x, bits).astype(np.uint16)


def sony_inv53_1d(low: np.ndarray, high: np.ndarray, axis: int) -> np.ndarray:
    """One-axis inverse 5/3 lifting in Sony's signed working domain."""

    lo = np.asarray(low, dtype=np.int32)
    hi = np.asarray(high, dtype=np.int32)
    if axis == 0:
        hi_prev = np.vstack([hi[:1], hi[:-1]])
        lo2 = lo - ((hi_prev + hi + 2) >> 2)
        lo_next = np.vstack([lo2[1:], lo2[-1:]])
        hi2 = hi + ((lo2 + lo_next) >> 1)
        out = np.empty((lo.shape[0] * 2, lo.shape[1]), dtype=np.int32)
        out[0::2] = lo2
        out[1::2] = hi2
        return out
    if axis == 1:
        hi_prev = np.concatenate([hi[:, :1], hi[:, :-1]], axis=1)
        lo2 = lo - ((hi_prev + hi + 2) >> 2)
        lo_next = np.concatenate([lo2[:, 1:], lo2[:, -1:]], axis=1)
        hi2 = hi + ((lo2 + lo_next) >> 1)
        out = np.empty((lo.shape[0], lo.shape[1] * 2), dtype=np.int32)
        out[:, 0::2] = lo2
        out[:, 1::2] = hi2
        return out
    raise ValueError("axis must be 0 or 1")


def synthesize_llvc3_level(ll: np.ndarray, sub0: np.ndarray, sub1: np.ndarray, sub2: np.ndarray) -> np.ndarray:
    """Synthesize one LLVC3 scale from LL plus three detail subbands.

    Packet component mapping, as verified against Imaging Edge for this ARW6
    sample:

    * sub0: horizontal detail for the low vertical branch (HL)
    * sub1: vertical detail for the low horizontal branch (LH)
    * sub2: diagonal detail (HH)

    The annoying bit is the extra flush row. At the bottom edge Sony feeds it
    across the branches: LH's final row comes from sub0[-1], HH's final row from
    sub1[-1].
    """

    ll_i = np.asarray(ll, dtype=np.int32)
    h = ll_i.shape[0]
    if ll_i.shape[1] != sub0.shape[1] or sub0.shape != sub1.shape or sub1.shape != sub2.shape:
        raise ValueError(f"unexpected subband shapes: ll={ll_i.shape}, sub0={sub0.shape}, sub1={sub1.shape}, sub2={sub2.shape}")
    if sub0.shape[0] < h + 1:
        raise ValueError(f"subbands need one flush row: ll={ll_i.shape}, sub={sub0.shape}")

    lh = np.empty((h, ll_i.shape[1]), dtype=np.int32)
    lh[:-1] = sub1[1:h]
    lh[-1] = sub0[h]

    hh = np.empty((h, ll_i.shape[1]), dtype=np.int32)
    hh[:-1] = sub2[1:h]
    hh[-1] = sub1[h]

    low_horizontal = sony_inv53_1d(ll_i, lh, axis=0)
    high_horizontal = sony_inv53_1d(sub0[:h], hh, axis=0)
    return sony_inv53_1d(low_horizontal, high_horizontal, axis=1)


def trunc_div2(x: np.ndarray) -> np.ndarray:
    """Integer division by two with C/C++ truncation toward zero."""

    a = np.asarray(x, dtype=np.int32)
    return np.where(a >= 0, a // 2, -((-a) // 2)).astype(np.int32)


def synthesize_llvc3_level_stride(
    ll: np.ndarray,
    sub0: np.ndarray,
    sub1: np.ndarray,
    sub2: np.ndarray,
    edge_rows: int,
    bottom_hh_extra: np.ndarray | None = None,
    edge_mode: str = "even",
) -> np.ndarray:
    """Same synthesis, with the larger line-flush padding used above group 1.

    \`edge_rows\` is the row multiplier minus one for the higher scales:
    group1 -> 0, group2 -> 1, group3 -> 2.  The group1 case falls back to the
    smaller helper above.
    """

    if edge_rows == 0:
        return synthesize_llvc3_level(ll, sub0, sub1, sub2)

    ll_i = np.asarray(ll, dtype=np.int32)
    h, w = ll_i.shape
    if sub0.shape != sub1.shape or sub1.shape != sub2.shape or sub0.shape[1] != w:
        raise ValueError(f"unexpected subband shapes: ll={ll_i.shape}, sub0={sub0.shape}, sub1={sub1.shape}, sub2={sub2.shape}")
    if sub0.shape[0] < h + edge_rows * 2:
        raise ValueError(f"not enough line-flush rows: ll={ll_i.shape}, sub={sub0.shape}, edge_rows={edge_rows}")

    hl = np.empty((h, w), dtype=np.int32)
    hl[:edge_rows] = sub0[:edge_rows]
    hl[edge_rows : h - edge_rows] = sub0[2 * edge_rows : h]
    hl[h - edge_rows :] = sub0[h : h + edge_rows]

    lh = np.empty((h, w), dtype=np.int32)
    lh[:edge_rows] = sub0[edge_rows : 2 * edge_rows]
    lh[edge_rows : h - edge_rows] = sub1[2 * edge_rows : h]
    lh[h - edge_rows :] = sub0[h + edge_rows : h + 2 * edge_rows]

    hh = np.empty((h, w), dtype=np.int32)

    def edge_fn(x: np.ndarray) -> np.ndarray:
        xi = x.astype(np.int32)
        signed_half_step = np.where(xi > 0, 1, np.where(xi < 0, -1, 0))
        if edge_mode == "even":
            mask = ((xi & 1) == 0) & (xi != 0)
        elif edge_mode == "odd":
            mask = (xi & 1) != 0
        else:
            raise ValueError(f"unknown LLVC3 edge_mode {edge_mode!r}")
        return (2 * xi + np.where(mask, signed_half_step, 0)).astype(np.int32)

    hh[:edge_rows] = edge_fn(sub1[:edge_rows])
    hh[edge_rows : h - edge_rows] = sub2[np.arange(edge_rows, h - edge_rows) + edge_rows]
    bottom_hh = edge_fn(sub1[h : h + edge_rows])
    if bottom_hh_extra is not None:
        extra = np.asarray(bottom_hh_extra, dtype=np.int32)
        if extra.shape != bottom_hh.shape:
            raise ValueError(f"bottom_hh_extra shape {extra.shape} != bottom edge {bottom_hh.shape}")
        bottom_hh = bottom_hh + extra
    hh[h - edge_rows :] = bottom_hh

    low_horizontal = sony_inv53_1d(ll_i, lh, axis=0)
    high_horizontal = sony_inv53_1d(hl, hh, axis=0)
    return sony_inv53_1d(low_horizontal, high_horizontal, axis=1)


def synthesize_llvc3_final_green(ll: np.ndarray, detail: np.ndarray) -> np.ndarray:
    """Final CFA-green reconstruction, from the 0x1ab570 path.

    Not the same 2-D 5/3 inverse used by groups 1..3. It expands the half-width
    green lowpass into both RGGB green sites. The row offsets look odd because
    the native line buffer keeps four guard rows at the top and a few latency
    rows at the bottom.
    """

    ll_i = np.asarray(ll, dtype=np.int32)
    det = np.asarray(detail, dtype=np.int32)
    h, w = ll_i.shape
    if det.shape[1] != w or det.shape[0] < h + 4:
        raise ValueError(f"unexpected final green shapes: ll={ll_i.shape}, detail={det.shape}")

    selected = np.empty((h, w), dtype=np.int32)
    top = min(4, h)
    selected[:top] = det[:top]
    if h > 4:
        selected[4:] = det[8 : 8 + (h - 4)]

    odd_green = np.empty((h, w), dtype=np.int32)
    for y in range(h):
        cur = selected[y]
        prev = selected[y - 1] if y > 0 else cur
        pred = np.empty(w, dtype=np.int32)
        pred[:-1] = (cur[1:] + prev[:-1] + cur[:-1] + prev[1:]) >> 2
        pred[-1] = ((prev[-1] + cur[-1]) * 2) >> 2
        odd_green[y] = ((2 * ll_i[y] - pred) >> 1).astype(np.int32)

    even_green = np.empty((h, w), dtype=np.int32)
    for y in range(h):
        cur = odd_green[y]
        nxt = odd_green[y + 1] if y + 1 < h else cur
        even_green[y, 0] = selected[y, 0] + (((cur[0] + nxt[0]) * 2) >> 2)
        if w > 1:
            even_green[y, 1:] = selected[y, 1:] + ((cur[:-1] + nxt[:-1] + nxt[1:] + cur[1:]) >> 2)

    out = np.empty((h, w * 2), dtype=np.int32)
    out[:, 0::2] = even_green
    out[:, 1::2] = odd_green
    return out


def finalize_llvc3_color_planes(
    v1_green: np.ndarray, v1_red: np.ndarray, v1_blue: np.ndarray, full_green: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Apply the final color-plane relation for the decoded output.

    Red and blue keep the group-3 residuals. Only the green predictor is swapped
    for the final CFA-green average:

        residual = (v1_color - v1_green) / 2
        v0_color = average(final_green_pair) + 2 * residual
    """

    g = np.asarray(v1_green, dtype=np.int32)
    r = np.asarray(v1_red, dtype=np.int32)
    b = np.asarray(v1_blue, dtype=np.int32)
    fg = np.asarray(full_green, dtype=np.int32)
    if fg.shape != (g.shape[0], g.shape[1] * 2):
        raise ValueError(f"unexpected final green shape: v1={g.shape}, full={fg.shape}")
    if r.shape != g.shape or b.shape != g.shape:
        raise ValueError(f"unexpected v1 color shapes: green={g.shape}, red={r.shape}, blue={b.shape}")

    gavg = (fg[:, 0::2] + fg[:, 1::2]) >> 1
    red_residual = (r - g) // 2
    blue_residual = (b - g) // 2
    return fg, gavg + 2 * red_residual, gavg + 2 * blue_residual


def signed_to_sample(x: np.ndarray, bits: int = 16, bias: int = INTERNAL_BIAS) -> np.ndarray:
    """Convert LLVC3 signed internal rows to Sony's unsigned output samples."""

    return clamp_sample(np.asarray(x, dtype=np.int32) + bias, bits).astype(np.uint16)


def apply_sample_lut(code_samples: np.ndarray, lut: np.ndarray) -> np.ndarray:
    """Map unsigned LLVC3 code-domain samples through a Sony sample LUT."""

    table = np.asarray(lut, dtype=np.uint16).reshape(-1)
    if table.size == 0:
        raise ValueError("sample LUT is empty")
    code = np.clip(np.asarray(code_samples, dtype=np.int32), 0, table.size - 1)
    return table[code].astype(np.uint16)


def clamp_signed_to_code_range(x: np.ndarray, max_code: int = 4095, bias: int = INTERNAL_BIAS) -> np.ndarray:
    """Clamp signed LLVC3 rows to Sony's 12-bit code range, then return signed rows."""

    return np.clip(np.asarray(x, dtype=np.int32) + bias, 0, max_code).astype(np.int32) - bias


def lifting_predict_detail(detail: np.ndarray, a: np.ndarray, b: np.ndarray, c: np.ndarray, d: np.ndarray) -> np.ndarray:
    """Prediction/update kernel from the 0x1ab570 scalar path.

    The four neighbor names are still placeholder-ish, but the integer
    operation itself is clear in the trace:

        pred = (a + b + c + d) >> 2
        out = (2 * detail - pred) >> 1
    """

    pred = (a.astype(np.int32) + b.astype(np.int32) + c.astype(np.int32) + d.astype(np.int32)) >> 2
    return ((2 * detail.astype(np.int32) - pred) >> 1).astype(np.int32)


def recombine_rggb(c0: np.ndarray, c1: np.ndarray, c2: np.ndarray) -> np.ndarray:
    """Recombine decoded LLVC planes into the TIFF-declared RGGB Bayer mosaic."""

    half_h, width = c0.shape
    if c1.shape != (half_h, width // 2) or c2.shape != (half_h, width // 2):
        raise ValueError(f"unexpected plane shapes: {c0.shape}, {c1.shape}, {c2.shape}")
    out = np.empty((half_h * 2, width), dtype=np.uint16)
    out[0::2, 0::2] = c1
    out[0::2, 1::2] = c0[:, 1::2]
    out[1::2, 0::2] = c0[:, 0::2]
    out[1::2, 1::2] = c2
    return out


def selftest() -> dict[str, object]:
    lo0 = np.array([1000, 1002, 10, 65530], dtype=np.uint16)
    lo1 = np.array([1002, 1004, 10, 65530], dtype=np.uint16)
    detail = np.array([0, 1, -20, 20], dtype=np.int32)
    merged = merge_average_detail(lo0, lo1, detail, bits=16)
    added = add_detail(lo0, detail, bits=16)
    doubled = add_double_detail(lo0, detail, bits=16)
    pred = lifting_predict_detail(
        detail,
        np.array([4, 8, 12, 16]),
        np.array([4, 8, 12, 16]),
        np.array([4, 8, 12, 16]),
        np.array([4, 8, 12, 16]),
    )
    return {
        "merge_average_detail": merged.tolist(),
        "add_detail": added.tolist(),
        "add_double_detail": doubled.tolist(),
        "lifting_predict_detail": pred.tolist(),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--out", default="")
    ns = ap.parse_args()
    result = selftest()
    text = json.dumps(result, indent=2)
    if ns.out:
        Path(ns.out).write_text(text, encoding="utf-8")
    print(text)


if __name__ == "__main__":
    main()
`;const ot=512,rn=Qt["Sony ILCE-7RM5"].split(/\s+/).map(Number).filter(Number.isFinite);let Be=null;function sn(r){if(r.byteLength<8)return null;const e=r.getUint16(0,!1);return e===18761?!0:e===19789?!1:null}function Tt(r,e,t){const n=Math.min(r.length,e+t);let i="";for(let s=e;s<n;s++){const a=r[s];if(a===0)break;i+=String.fromCharCode(a)}return i.trim()}function U(r,e,t,n,i){const s=t===1||t===2||t===7?1:t===3||t===8?2:t===4||t===9?4:0;if(!s)return[];const a=s*n,o=a<=4?i:r.getUint32(i,e);if(o<0||o+a>r.byteLength)return[];const c=[];for(let u=0;u<n;u++){const d=o+u*s;t===1||t===2||t===7?c.push(r.getUint8(d)):t===3?c.push(r.getUint16(d,e)):t===8?c.push(r.getInt16(d,e)):t===4?c.push(r.getUint32(d,e)):t===9&&c.push(r.getInt32(d,e))}return c}function lt(r,e,t,n,i,s){if(n!==2||i<=0)return"";const a=i<=4?s:e.getUint32(s,t);return a<0||a>=r.length?"":Tt(r,a,i)}function Rt(r){const e=new Uint8Array(r),t=new DataView(r),n=sn(t);if(n===null||t.getUint16(2,n)!==42)return null;const s=l=>t.getUint16(l,n),a=l=>t.getUint32(l,n),o=[a(4)],c=new Set;let u="",d="";for(;o.length;){const l=o.pop();if(c.has(l)||l<=0||l+2>t.byteLength)continue;c.add(l);const f=s(l);if(l+2+f*12+4>t.byteLength)continue;const h=new Map;for(let S=0;S<f;S++){const k=l+2+S*12,m=s(k),O=s(k+2),fe=a(k+4),V=k+8;h.set(m,{type:O,count:fe,valueOffset:V})}const p=h.get(271),y=h.get(272);p&&!u&&(u=lt(e,t,n,p.type,p.count,p.valueOffset)),y&&!d&&(d=lt(e,t,n,y.type,y.count,y.valueOffset));const w=h.get(330);if(w){const S=U(t,n,w.type,w.count,w.valueOffset);for(const k of S)o.push(k)}const b=h.get(259),C=h.get(262);if(b&&C){const S=U(t,n,b.type,b.count,b.valueOffset)[0],k=U(t,n,C.type,C.count,C.valueOffset)[0];if(S===32766&&k===32803){const m=U(t,n,h.get(256).type,h.get(256).count,h.get(256).valueOffset)[0],O=U(t,n,h.get(257).type,h.get(257).count,h.get(257).valueOffset)[0],fe=U(t,n,h.get(258).type,h.get(258).count,h.get(258).valueOffset)[0],V=U(t,n,h.get(273).type,h.get(273).count,h.get(273).valueOffset)[0],oe=U(t,n,h.get(279).type,h.get(279).count,h.get(279).valueOffset)[0],le=h.get(33422)?U(t,n,h.get(33422).type,h.get(33422).count,h.get(33422).valueOffset):[0,1,1,2],N=h.get(29456)?U(t,n,h.get(29456).type,h.get(29456).count,h.get(29456).valueOffset):[],ye=h.get(50717)?U(t,n,h.get(50717).type,h.get(50717).count,h.get(50717).valueOffset)[0]:16383,ce=h.get(50719)?U(t,n,h.get(50719).type,h.get(50719).count,h.get(50719).valueOffset):[],_=h.get(50720)?U(t,n,h.get(50720).type,h.get(50720).count,h.get(50720).valueOffset):[];if(V+ot+16>e.length||V+oe>e.length)return null;const E=V+ot,Q=Tt(e,E,4),D=e[E+8]<<8|e[E+9],A=e[E+10]<<8|e[E+11],j=e[E+12]<<8|e[E+13],be=e[E+14]<<8|e[E+15],Wt=j>>4&63,Xt=be>>13,Kt=be>>10&3;if(Q!=="A000"||D!==m||A*2!==O||Wt!==16||Xt!==3||Kt!==3)return null;const qt=[Number(N[0]??512),Number(N[1]??N[0]??512),Number(N[2]??N[0]??512),Number(N[3]??N[0]??512)],Yt=le.slice(0,4).map(at=>at===0?"R":at===2?"B":"G").join("")||"RGGB";return{width:m,height:O,bitsPerSample:fe,compression:S,photometric:k,blackLevel:qt,whiteLevel:Number(ye||16383),cfaPattern:Yt,defaultCropOrigin:ce.length>=2?[Number(ce[0]),Number(ce[1])]:void 0,defaultCropSize:_.length>=2?[Number(_[0]),Number(_[1])]:void 0,make:u||"SONY",model:d||"ILCE-7M5"}}}const R=a(l+2+f*12);R&&o.push(R)}return null}async function an(r){return Be||(Be=(async()=>{const e=await Zt();return e.__jtrSonyCrawHqDecoderReady||(await e.FS.mkdirTree("/sony_craw_hq"),await e.FS.writeFile("/sony_craw_hq/llvc3_bitstream_probe.py",en),await e.FS.writeFile("/sony_craw_hq/llvc3_entropy.py",tn),await e.FS.writeFile("/sony_craw_hq/llvc3_math.py",nn),await e.runPythonAsync(`
import sys
if "/sony_craw_hq" not in sys.path:
    sys.path.insert(0, "/sony_craw_hq")
from pathlib import Path
import numpy as np
from llvc3_bitstream_probe import find_raw_subifd, parse_llvc_header, RAW_STREAM_OFFSET
from llvc3_entropy import decode_packet_arrays, integrate_type1_coefficients
from llvc3_math import apply_sample_lut, clamp_signed_to_code_range, recombine_rggb, signed_to_sample
from llvc3_math import finalize_llvc3_color_planes, synthesize_llvc3_final_green, synthesize_llvc3_level_stride

def jtr_decode_signed_planes(arw):
    raw_info, _strip = find_raw_subifd(arw)
    low_rows = raw_info.height // 16

    g0, _meta = decode_packet_arrays(arw, 0, 0)
    green = integrate_type1_coefficients(g0[0][:low_rows], 2048) - 2048

    r0, _meta = decode_packet_arrays(arw, 0, 1)
    red_residual = integrate_type1_coefficients(r0[0][:low_rows], 0)

    b0, _meta = decode_packet_arrays(arw, 0, 2)
    blue_residual = integrate_type1_coefficients(b0[0][:low_rows], 0)

    for group, edge_rows in ((1, 0), (2, 1), (3, 2)):
        old_green = green
        old_red_residual = red_residual
        old_blue_residual = blue_residual

        planes, _meta = decode_packet_arrays(arw, group, 0)
        green = synthesize_llvc3_level_stride(old_green, planes[0], planes[1], planes[2], edge_rows)

        planes, _meta = decode_packet_arrays(arw, group, 1)
        edge_mode = "odd" if group == 3 else "even"
        red_residual = synthesize_llvc3_level_stride(
            old_red_residual, planes[0], planes[1], planes[2], edge_rows, edge_mode=edge_mode
        )

        planes, _meta = decode_packet_arrays(arw, group, 2)
        blue_residual = synthesize_llvc3_level_stride(
            old_blue_residual, planes[0], planes[1], planes[2], edge_rows, edge_mode=edge_mode
        )

    g4, _meta = decode_packet_arrays(arw, 4, 0)
    full_green = synthesize_llvc3_final_green(green, g4[0])
    v1_red = green + 2 * red_residual
    v1_blue = green + 2 * blue_residual
    return finalize_llvc3_color_planes(green, v1_red, v1_blue, full_green)

def jtr_decode_sony_craw_hq(arw_bytes, lut_bytes=None):
    path = Path("/tmp/jtr_sony_craw_hq_input.arw")
    path.write_bytes(bytes(arw_bytes))
    raw_info, strip = find_raw_subifd(path)
    header = parse_llvc_header(strip[RAW_STREAM_OFFSET:])
    if header.coded_width != raw_info.width or header.logical_height != raw_info.height or header.component_count != 3:
        raise ValueError(f"unexpected ARW6/LLVC3 header/raw mismatch: raw={raw_info}, header={header}")
    if raw_info.width % 16 or raw_info.height % 16:
        raise ValueError(f"decoder expects dimensions divisible by 16, got {raw_info.width}x{raw_info.height}")
    signed_c0, signed_c1, signed_c2 = jtr_decode_signed_planes(path)
    if lut_bytes is not None:
        lut = np.frombuffer(bytes(lut_bytes), dtype="<u2")
        if lut.size:
            if lut.size < 65536:
                lut = np.pad(lut, (0, 65536 - lut.size), constant_values=int(lut[-1]))
            lut = lut[:65536].astype(np.uint16)
            sample_c0 = apply_sample_lut(signed_to_sample(clamp_signed_to_code_range(signed_c0)), lut)
            sample_c1 = apply_sample_lut(signed_to_sample(clamp_signed_to_code_range(signed_c1)), lut)
            sample_c2 = apply_sample_lut(signed_to_sample(clamp_signed_to_code_range(signed_c2)), lut)
        else:
            sample_c0 = signed_to_sample(signed_c0)
            sample_c1 = signed_to_sample(signed_c1)
            sample_c2 = signed_to_sample(signed_c2)
    else:
        sample_c0 = signed_to_sample(signed_c0)
        sample_c1 = signed_to_sample(signed_c1)
        sample_c2 = signed_to_sample(signed_c2)
    raw = recombine_rggb(sample_c0, sample_c1, sample_c2)
    return raw.astype("<u2", copy=False).tobytes()
`),e.__jtrSonyCrawHqDecoderReady=!0),e})()),Be}function on(r){return Rt(r)}async function ln(r,e){const t=typeof performance<"u"?performance.now():Date.now(),n=Rt(r);if(!n)return null;const i=await an(),s=typeof performance<"u"?performance.now():Date.now(),a=new Uint8Array(r),o=await fetch(new URL("/assets/sony_llvc3_static_lut4096_padded_u16-FsVBk-IV.bin",import.meta.url));if(!o.ok)throw new Error(`Failed to load Sony LLVC3 sample LUT: HTTP ${o.status}`);const c=new Uint8Array(await o.arrayBuffer()),u=typeof performance<"u"?performance.now():Date.now();i.globals.set("jtr_sony_arw_bytes",a),i.globals.set("jtr_sony_lut_bytes",c);const d=await i.runPythonAsync("jtr_decode_sony_craw_hq(jtr_sony_arw_bytes.to_py(), jtr_sony_lut_bytes.to_py())"),l=typeof performance<"u"?performance.now():Date.now(),f=d.toJs();typeof d.destroy=="function"&&d.destroy(),i.globals.delete("jtr_sony_arw_bytes"),i.globals.delete("jtr_sony_lut_bytes");const h=new Uint8Array(f.byteLength);h.set(f);const p=new Uint16Array(h.buffer),y=typeof performance<"u"?performance.now():Date.now();if(p.length!==n.width*n.height)throw new Error(`Sony cRAW HQ decoded size mismatch: got ${p.length}, expected ${n.width*n.height}`);const w=n.model||"ILCE-7M5",b=w.startsWith("Sony ")?w:`Sony ${w}`,C=w==="ILCE-7M5"?rn:null,R={make:n.make||"SONY",model:w,camera_make:n.make||"SONY",camera_model:w,UniqueCameraModel:b,sourceFormat:"Sony cRAW HQ / LLVC3",sonyCrawHq:n,color_desc:n.cfaPattern,black_level_per_channel:n.blackLevel,white_level:n.whiteLevel,color_matrix:C&&C.length===9?C:void 0,idata:{filters:2492765332,colors:3},color_data:{cblack_rawpy_style:n.blackLevel,dng_levels:{dng_cblack:n.blackLevel,dng_whitelevel:n.whiteLevel}}},S=typeof performance<"u"?performance.now():Date.now();return console.info("[Sony cRAW HQ] decode timings",{width:n.width,height:n.height,pyodideReadyMs:Math.round(s-t),lutLoadMs:Math.round(u-s),llvc3DecodeMs:Math.round(l-u),copyMs:Math.round(y-l),totalMs:Math.round(S-t)}),{rawImageData:{data:p,width:n.width,height:n.height,bayerPattern:n.cfaPattern,blackLevels:n.blackLevel,whiteLevel:n.whiteLevel,metadata:R,isThreePlane:!1,isXTrans:!1},info:n}}async function cn(r,e){return ln(r)}var De=typeof self<"u"?self:global;const _e=typeof navigator<"u",dn=_e&&typeof HTMLImageElement>"u",Pe=!(typeof global>"u"||typeof process>"u"||!process.versions||!process.versions.node),Me=De.Buffer,we=De.BigInt,Fe=!!Me,hn=r=>r;function Ie(r,e=hn){if(Pe)try{return typeof require=="function"?Promise.resolve(e(require(r))):import(r).then(e)}catch{console.warn(`Couldn't load ${r}`)}}let Qe=De.fetch;const un=r=>Qe=r;if(!De.fetch){const r=Ie("http",n=>n),e=Ie("https",n=>n),t=(n,{headers:i}={})=>new Promise(async(s,a)=>{let{port:o,hostname:c,pathname:u,protocol:d,search:l}=new URL(n);const f={method:"GET",hostname:c,path:encodeURI(u)+l,headers:i};o!==""&&(f.port=Number(o));const h=(d==="https:"?await e:await r).request(f,p=>{if(p.statusCode===301||p.statusCode===302){let y=new URL(p.headers.location,n).toString();return t(y,{headers:i}).then(s).catch(a)}s({status:p.statusCode,arrayBuffer:()=>new Promise(y=>{let w=[];p.on("data",b=>w.push(b)),p.on("end",()=>y(Buffer.concat(w)))})})});h.on("error",a),h.end()});un(t)}function g(r,e,t){return e in r?Object.defineProperty(r,e,{value:t,enumerable:!0,configurable:!0,writable:!0}):r[e]=t,r}const Ae=r=>Ot(r)?void 0:r,fn=r=>r!==void 0;function Ot(r){return r===void 0||(r instanceof Map?r.size===0:Object.values(r).filter(fn).length===0)}function I(r){let e=new Error(r);throw delete e.stack,e}function Z(r){return(r=function(e){for(;e.endsWith("\0");)e=e.slice(0,-1);return e}(r).trim())===""?void 0:r}function Ke(r){let e=function(t){let n=0;return t.ifd0.enabled&&(n+=1024),t.exif.enabled&&(n+=2048),t.makerNote&&(n+=2048),t.userComment&&(n+=1024),t.gps.enabled&&(n+=512),t.interop.enabled&&(n+=100),t.ifd1.enabled&&(n+=1024),n+2048}(r);return r.jfif.enabled&&(e+=50),r.xmp.enabled&&(e+=2e4),r.iptc.enabled&&(e+=14e3),r.icc.enabled&&(e+=6e3),e}const qe=r=>String.fromCharCode.apply(null,r),ct=typeof TextDecoder<"u"?new TextDecoder("utf-8"):void 0;function Dt(r){return ct?ct.decode(r):Fe?Buffer.from(r).toString("utf8"):decodeURIComponent(escape(qe(r)))}class L{static from(e,t){return e instanceof this&&e.le===t?e:new L(e,void 0,void 0,t)}constructor(e,t=0,n,i){if(typeof i=="boolean"&&(this.le=i),Array.isArray(e)&&(e=new Uint8Array(e)),e===0)this.byteOffset=0,this.byteLength=0;else if(e instanceof ArrayBuffer){n===void 0&&(n=e.byteLength-t);let s=new DataView(e,t,n);this._swapDataView(s)}else if(e instanceof Uint8Array||e instanceof DataView||e instanceof L){n===void 0&&(n=e.byteLength-t),(t+=e.byteOffset)+n>e.byteOffset+e.byteLength&&I("Creating view outside of available memory in ArrayBuffer");let s=new DataView(e.buffer,t,n);this._swapDataView(s)}else if(typeof e=="number"){let s=new DataView(new ArrayBuffer(e));this._swapDataView(s)}else I("Invalid input argument for BufferView: "+e)}_swapArrayBuffer(e){this._swapDataView(new DataView(e))}_swapBuffer(e){this._swapDataView(new DataView(e.buffer,e.byteOffset,e.byteLength))}_swapDataView(e){this.dataView=e,this.buffer=e.buffer,this.byteOffset=e.byteOffset,this.byteLength=e.byteLength}_lengthToEnd(e){return this.byteLength-e}set(e,t,n=L){return e instanceof DataView||e instanceof L?e=new Uint8Array(e.buffer,e.byteOffset,e.byteLength):e instanceof ArrayBuffer&&(e=new Uint8Array(e)),e instanceof Uint8Array||I("BufferView.set(): Invalid data argument."),this.toUint8().set(e,t),new n(this,t,e.byteLength)}subarray(e,t){return t=t||this._lengthToEnd(e),new L(this,e,t)}toUint8(){return new Uint8Array(this.buffer,this.byteOffset,this.byteLength)}getUint8Array(e,t){return new Uint8Array(this.buffer,this.byteOffset+e,t)}getString(e=0,t=this.byteLength){return Dt(this.getUint8Array(e,t))}getLatin1String(e=0,t=this.byteLength){let n=this.getUint8Array(e,t);return qe(n)}getUnicodeString(e=0,t=this.byteLength){const n=[];for(let i=0;i<t&&e+i<this.byteLength;i+=2)n.push(this.getUint16(e+i));return qe(n)}getInt8(e){return this.dataView.getInt8(e)}getUint8(e){return this.dataView.getUint8(e)}getInt16(e,t=this.le){return this.dataView.getInt16(e,t)}getInt32(e,t=this.le){return this.dataView.getInt32(e,t)}getUint16(e,t=this.le){return this.dataView.getUint16(e,t)}getUint32(e,t=this.le){return this.dataView.getUint32(e,t)}getFloat32(e,t=this.le){return this.dataView.getFloat32(e,t)}getFloat64(e,t=this.le){return this.dataView.getFloat64(e,t)}getFloat(e,t=this.le){return this.dataView.getFloat32(e,t)}getDouble(e,t=this.le){return this.dataView.getFloat64(e,t)}getUintBytes(e,t,n){switch(t){case 1:return this.getUint8(e,n);case 2:return this.getUint16(e,n);case 4:return this.getUint32(e,n);case 8:return this.getUint64&&this.getUint64(e,n)}}getUint(e,t,n){switch(t){case 8:return this.getUint8(e,n);case 16:return this.getUint16(e,n);case 32:return this.getUint32(e,n);case 64:return this.getUint64&&this.getUint64(e,n)}}toString(e){return this.dataView.toString(e,this.constructor.name)}ensureChunk(){}}function Ye(r,e){I(`${r} '${e}' was not loaded, try using full build of exifr.`)}class Ze extends Map{constructor(e){super(),this.kind=e}get(e,t){return this.has(e)||Ye(this.kind,e),t&&(e in t||function(n,i){I(`Unknown ${n} '${i}'.`)}(this.kind,e),t[e].enabled||Ye(this.kind,e)),super.get(e)}keyList(){return Array.from(this.keys())}}var G=new Ze("file parser"),P=new Ze("segment parser"),H=new Ze("file reader");function pn(r,e){return typeof r=="string"?dt(r,e):_e&&!dn&&r instanceof HTMLImageElement?dt(r.src,e):r instanceof Uint8Array||r instanceof ArrayBuffer||r instanceof DataView?new L(r):_e&&r instanceof Blob?Je(r,e,"blob",ie):void I("Invalid input argument")}function dt(r,e){return(t=r).startsWith("data:")||t.length>1e4?$e(r,e,"base64"):Pe&&r.includes("://")?Je(r,e,"url",ne):Pe?$e(r,e,"fs"):_e?Je(r,e,"url",ne):void I("Invalid input argument");var t}async function Je(r,e,t,n){return H.has(t)?$e(r,e,t):n?async function(i,s){let a=await s(i);return new L(a)}(r,n):void I(`Parser ${t} is not loaded`)}async function $e(r,e,t){let n=new(H.get(t))(r,e);return await n.read(),n}const ne=r=>Qe(r).then(e=>e.arrayBuffer()),ie=r=>new Promise((e,t)=>{let n=new FileReader;n.onloadend=()=>e(n.result||new ArrayBuffer),n.onerror=t,n.readAsArrayBuffer(r)});class gn extends Map{get tagKeys(){return this.allKeys||(this.allKeys=Array.from(this.keys())),this.allKeys}get tagValues(){return this.allValues||(this.allValues=Array.from(this.values())),this.allValues}}function v(r,e,t){let n=new gn;for(let[i,s]of t)n.set(i,s);if(Array.isArray(e))for(let i of e)r.set(i,n);else r.set(e,n);return n}function re(r,e,t){let n,i=r.get(e);for(n of t)i.set(n[0],n[1])}const T=new Map,z=new Map,J=new Map,X=["chunked","firstChunkSize","firstChunkSizeNode","firstChunkSizeBrowser","chunkSize","chunkLimit"],he=["jfif","xmp","icc","iptc","ihdr"],se=["tiff",...he],x=["ifd0","ifd1","exif","gps","interop"],K=[...se,...x],q=["makerNote","userComment"],ue=["translateKeys","translateValues","reviveValues","multiSegment"],Y=[...ue,"sanitize","mergeOutput","silentErrors"];class Mt{get translate(){return this.translateKeys||this.translateValues||this.reviveValues}}class pe extends Mt{get needed(){return this.enabled||this.deps.size>0}constructor(e,t,n,i){if(super(),g(this,"enabled",!1),g(this,"skip",new Set),g(this,"pick",new Set),g(this,"deps",new Set),g(this,"translateKeys",!1),g(this,"translateValues",!1),g(this,"reviveValues",!1),this.key=e,this.enabled=t,this.parse=this.enabled,this.applyInheritables(i),this.canBeFiltered=x.includes(e),this.canBeFiltered&&(this.dict=T.get(e)),n!==void 0)if(Array.isArray(n))this.parse=this.enabled=!0,this.canBeFiltered&&n.length>0&&this.translateTagSet(n,this.pick);else if(typeof n=="object"){if(this.enabled=!0,this.parse=n.parse!==!1,this.canBeFiltered){let{pick:s,skip:a}=n;s&&s.length>0&&this.translateTagSet(s,this.pick),a&&a.length>0&&this.translateTagSet(a,this.skip)}this.applyInheritables(n)}else n===!0||n===!1?this.parse=this.enabled=n:I(`Invalid options argument: ${n}`)}applyInheritables(e){let t,n;for(t of ue)n=e[t],n!==void 0&&(this[t]=n)}translateTagSet(e,t){if(this.dict){let n,i,{tagKeys:s,tagValues:a}=this.dict;for(n of e)typeof n=="string"?(i=a.indexOf(n),i===-1&&(i=s.indexOf(Number(n))),i!==-1&&t.add(Number(s[i]))):t.add(n)}else for(let n of e)t.add(n)}finalizeFilters(){!this.enabled&&this.deps.size>0?(this.enabled=!0,Te(this.pick,this.deps)):this.enabled&&this.pick.size>0&&Te(this.pick,this.deps)}}var F={jfif:!1,tiff:!0,xmp:!1,icc:!1,iptc:!1,ifd0:!0,ifd1:!1,exif:!0,gps:!0,interop:!1,ihdr:void 0,makerNote:!1,userComment:!1,multiSegment:!1,skip:[],pick:[],translateKeys:!0,translateValues:!0,reviveValues:!0,sanitize:!0,mergeOutput:!0,silentErrors:!0,chunked:!0,firstChunkSize:void 0,firstChunkSizeNode:512,firstChunkSizeBrowser:65536,chunkSize:65536,chunkLimit:5},ht=new Map;class ae extends Mt{static useCached(e){let t=ht.get(e);return t!==void 0||(t=new this(e),ht.set(e,t)),t}constructor(e){super(),e===!0?this.setupFromTrue():e===void 0?this.setupFromUndefined():Array.isArray(e)?this.setupFromArray(e):typeof e=="object"?this.setupFromObject(e):I(`Invalid options argument ${e}`),this.firstChunkSize===void 0&&(this.firstChunkSize=_e?this.firstChunkSizeBrowser:this.firstChunkSizeNode),this.mergeOutput&&(this.ifd1.enabled=!1),this.filterNestedSegmentTags(),this.traverseTiffDependencyTree(),this.checkLoadedPlugins()}setupFromUndefined(){let e;for(e of X)this[e]=F[e];for(e of Y)this[e]=F[e];for(e of q)this[e]=F[e];for(e of K)this[e]=new pe(e,F[e],void 0,this)}setupFromTrue(){let e;for(e of X)this[e]=F[e];for(e of Y)this[e]=F[e];for(e of q)this[e]=!0;for(e of K)this[e]=new pe(e,!0,void 0,this)}setupFromArray(e){let t;for(t of X)this[t]=F[t];for(t of Y)this[t]=F[t];for(t of q)this[t]=F[t];for(t of K)this[t]=new pe(t,!1,void 0,this);this.setupGlobalFilters(e,void 0,x)}setupFromObject(e){let t;for(t of(x.ifd0=x.ifd0||x.image,x.ifd1=x.ifd1||x.thumbnail,Object.assign(this,e),X))this[t]=ze(e[t],F[t]);for(t of Y)this[t]=ze(e[t],F[t]);for(t of q)this[t]=ze(e[t],F[t]);for(t of se)this[t]=new pe(t,F[t],e[t],this);for(t of x)this[t]=new pe(t,F[t],e[t],this.tiff);this.setupGlobalFilters(e.pick,e.skip,x,K),e.tiff===!0?this.batchEnableWithBool(x,!0):e.tiff===!1?this.batchEnableWithUserValue(x,e):Array.isArray(e.tiff)?this.setupGlobalFilters(e.tiff,void 0,x):typeof e.tiff=="object"&&this.setupGlobalFilters(e.tiff.pick,e.tiff.skip,x)}batchEnableWithBool(e,t){for(let n of e)this[n].enabled=t}batchEnableWithUserValue(e,t){for(let n of e){let i=t[n];this[n].enabled=i!==!1&&i!==void 0}}setupGlobalFilters(e,t,n,i=n){if(e&&e.length){for(let a of i)this[a].enabled=!1;let s=ut(e,n);for(let[a,o]of s)Te(this[a].pick,o),this[a].enabled=!0}else if(t&&t.length){let s=ut(t,n);for(let[a,o]of s)Te(this[a].skip,o)}}filterNestedSegmentTags(){let{ifd0:e,exif:t,xmp:n,iptc:i,icc:s}=this;this.makerNote?t.deps.add(37500):t.skip.add(37500),this.userComment?t.deps.add(37510):t.skip.add(37510),n.enabled||e.skip.add(700),i.enabled||e.skip.add(33723),s.enabled||e.skip.add(34675)}traverseTiffDependencyTree(){let{ifd0:e,exif:t,gps:n,interop:i}=this;i.needed&&(t.deps.add(40965),e.deps.add(40965)),t.needed&&e.deps.add(34665),n.needed&&e.deps.add(34853),this.tiff.enabled=x.some(s=>this[s].enabled===!0)||this.makerNote||this.userComment;for(let s of x)this[s].finalizeFilters()}get onlyTiff(){return!he.map(e=>this[e].enabled).some(e=>e===!0)&&this.tiff.enabled}checkLoadedPlugins(){for(let e of se)this[e].enabled&&!P.has(e)&&Ye("segment parser",e)}}function ut(r,e){let t,n,i,s,a=[];for(i of e){for(s of(t=T.get(i),n=[],t))(r.includes(s[0])||r.includes(s[1]))&&n.push(s[0]);n.length&&a.push([i,n])}return a}function ze(r,e){return r!==void 0?r:e!==void 0?e:void 0}function Te(r,e){for(let t of e)r.add(t)}g(ae,"default",F);class ${constructor(e){g(this,"parsers",{}),g(this,"output",{}),g(this,"errors",[]),g(this,"pushToErrors",t=>this.errors.push(t)),this.options=ae.useCached(e)}async read(e){this.file=await pn(e,this.options)}setup(){if(this.fileParser)return;let{file:e}=this,t=e.getUint16(0);for(let[n,i]of G)if(i.canHandle(e,t))return this.fileParser=new i(this.options,this.file,this.parsers),e[n]=!0;this.file.close&&this.file.close(),I("Unknown file format")}async parse(){let{output:e,errors:t}=this;return this.setup(),this.options.silentErrors?(await this.executeParsers().catch(this.pushToErrors),t.push(...this.fileParser.errors)):await this.executeParsers(),this.file.close&&this.file.close(),this.options.silentErrors&&t.length>0&&(e.errors=t),Ae(e)}async executeParsers(){let{output:e}=this;await this.fileParser.parse();let t=Object.values(this.parsers).map(async n=>{let i=await n.parse();n.assignToOutput(e,i)});this.options.silentErrors&&(t=t.map(n=>n.catch(this.pushToErrors))),await Promise.all(t)}async extractThumbnail(){this.setup();let{options:e,file:t}=this,n=P.get("tiff",e);var i;if(t.tiff?i={start:0,type:"tiff"}:t.jpeg&&(i=await this.fileParser.getOrFindSegment("tiff")),i===void 0)return;let s=await this.fileParser.ensureSegmentChunk(i),a=this.parsers.tiff=new n(s,e,t),o=await a.extractThumbnail();return t.close&&t.close(),o}}async function Le(r,e){let t=new $(e);return await t.read(r),t.parse()}var mn=Object.freeze({__proto__:null,parse:Le,Exifr:$,fileParsers:G,segmentParsers:P,fileReaders:H,tagKeys:T,tagValues:z,tagRevivers:J,createDictionary:v,extendDictionary:re,fetchUrlAsArrayBuffer:ne,readBlobAsArrayBuffer:ie,chunkedProps:X,otherSegments:he,segments:se,tiffBlocks:x,segmentsAndBlocks:K,tiffExtractables:q,inheritables:ue,allFormatters:Y,Options:ae});class Ee{constructor(e,t,n){g(this,"errors",[]),g(this,"ensureSegmentChunk",async i=>{let s=i.start,a=i.size||65536;if(this.file.chunked)if(this.file.available(s,a))i.chunk=this.file.subarray(s,a);else try{i.chunk=await this.file.readChunk(s,a)}catch(o){I(`Couldn't read segment: ${JSON.stringify(i)}. ${o.message}`)}else this.file.byteLength>s+a?i.chunk=this.file.subarray(s,a):i.size===void 0?i.chunk=this.file.subarray(s):I("Segment unreachable: "+JSON.stringify(i));return i.chunk}),this.extendOptions&&this.extendOptions(e),this.options=e,this.file=t,this.parsers=n}injectSegment(e,t){this.options[e].enabled&&this.createParser(e,t)}createParser(e,t){let n=new(P.get(e))(t,this.options,this.file);return this.parsers[e]=n}createParsers(e){for(let t of e){let{type:n,chunk:i}=t,s=this.options[n];if(s&&s.enabled){let a=this.parsers[n];a&&a.append||a||this.createParser(n,i)}}}async readSegments(e){let t=e.map(this.ensureSegmentChunk);await Promise.all(t)}}class B{static findPosition(e,t){let n=e.getUint16(t+2)+2,i=typeof this.headerLength=="function"?this.headerLength(e,t,n):this.headerLength,s=t+i,a=n-i;return{offset:t,length:n,headerLength:i,start:s,size:a,end:s+a}}static parse(e,t={}){return new this(e,new ae({[this.type]:t}),e).parse()}normalizeInput(e){return e instanceof L?e:new L(e)}constructor(e,t={},n){g(this,"errors",[]),g(this,"raw",new Map),g(this,"handleError",i=>{if(!this.options.silentErrors)throw i;this.errors.push(i.message)}),this.chunk=this.normalizeInput(e),this.file=n,this.type=this.constructor.type,this.globalOptions=this.options=t,this.localOptions=t[this.type],this.canTranslate=this.localOptions&&this.localOptions.translate}translate(){this.canTranslate&&(this.translated=this.translateBlock(this.raw,this.type))}get output(){return this.translated?this.translated:this.raw?Object.fromEntries(this.raw):void 0}translateBlock(e,t){let n=J.get(t),i=z.get(t),s=T.get(t),a=this.options[t],o=a.reviveValues&&!!n,c=a.translateValues&&!!i,u=a.translateKeys&&!!s,d={};for(let[l,f]of e)o&&n.has(l)?f=n.get(l)(f):c&&i.has(l)&&(f=this.translateValue(f,i.get(l))),u&&s.has(l)&&(l=s.get(l)||l),d[l]=f;return d}translateValue(e,t){return t[e]||t.DEFAULT||e}assignToOutput(e,t){this.assignObjectToOutput(e,this.constructor.type,t)}assignObjectToOutput(e,t,n){if(this.globalOptions.mergeOutput)return Object.assign(e,n);e[t]?Object.assign(e[t],n):e[t]=n}}g(B,"headerLength",4),g(B,"type",void 0),g(B,"multiSegment",!1),g(B,"canHandle",()=>!1);function _n(r){return r===192||r===194||r===196||r===219||r===221||r===218||r===254}function yn(r){return r>=224&&r<=239}function bn(r,e,t){for(let[n,i]of P)if(i.canHandle(r,e,t))return n}class ft extends Ee{constructor(...e){super(...e),g(this,"appSegments",[]),g(this,"jpegSegments",[]),g(this,"unknownSegments",[])}static canHandle(e,t){return t===65496}async parse(){await this.findAppSegments(),await this.readSegments(this.appSegments),this.mergeMultiSegments(),this.createParsers(this.mergedAppSegments||this.appSegments)}setupSegmentFinderArgs(e){e===!0?(this.findAll=!0,this.wanted=new Set(P.keyList())):(e=e===void 0?P.keyList().filter(t=>this.options[t].enabled):e.filter(t=>this.options[t].enabled&&P.has(t)),this.findAll=!1,this.remaining=new Set(e),this.wanted=new Set(e)),this.unfinishedMultiSegment=!1}async findAppSegments(e=0,t){this.setupSegmentFinderArgs(t);let{file:n,findAll:i,wanted:s,remaining:a}=this;if(!i&&this.file.chunked&&(i=Array.from(s).some(o=>{let c=P.get(o),u=this.options[o];return c.multiSegment&&u.multiSegment}),i&&await this.file.readWhole()),e=this.findAppSegmentsInRange(e,n.byteLength),!this.options.onlyTiff&&n.chunked){let o=!1;for(;a.size>0&&!o&&(n.canReadNextChunk||this.unfinishedMultiSegment);){let{nextChunkOffset:c}=n,u=this.appSegments.some(d=>!this.file.available(d.offset||d.start,d.length||d.size));if(o=e>c&&!u?!await n.readNextChunk(e):!await n.readNextChunk(c),(e=this.findAppSegmentsInRange(e,n.byteLength))===void 0)return}}}findAppSegmentsInRange(e,t){t-=2;let n,i,s,a,o,c,{file:u,findAll:d,wanted:l,remaining:f,options:h}=this;for(;e<t;e++)if(u.getUint8(e)===255){if(n=u.getUint8(e+1),yn(n)){if(i=u.getUint16(e+2),s=bn(u,e,i),s&&l.has(s)&&(a=P.get(s),o=a.findPosition(u,e),c=h[s],o.type=s,this.appSegments.push(o),!d&&(a.multiSegment&&c.multiSegment?(this.unfinishedMultiSegment=o.chunkNumber<o.chunkCount,this.unfinishedMultiSegment||f.delete(s)):f.delete(s),f.size===0)))break;h.recordUnknownSegments&&(o=B.findPosition(u,e),o.marker=n,this.unknownSegments.push(o)),e+=i+1}else if(_n(n)){if(i=u.getUint16(e+2),n===218&&h.stopAfterSos!==!1)return;h.recordJpegSegments&&this.jpegSegments.push({offset:e,length:i,marker:n}),e+=i+1}}return e}mergeMultiSegments(){if(!this.appSegments.some(t=>t.multiSegment))return;let e=function(t,n){let i,s,a,o=new Map;for(let c=0;c<t.length;c++)i=t[c],s=i[n],o.has(s)?a=o.get(s):o.set(s,a=[]),a.push(i);return Array.from(o)}(this.appSegments,"type");this.mergedAppSegments=e.map(([t,n])=>{let i=P.get(t,this.options);return i.handleMultiSegments?{type:t,chunk:i.handleMultiSegments(n)}:n[0]})}getSegment(e){return this.appSegments.find(t=>t.type===e)}async getOrFindSegment(e){let t=this.getSegment(e);return t===void 0&&(await this.findAppSegments(0,[e]),t=this.getSegment(e)),t}}g(ft,"type","jpeg"),G.set("jpeg",ft);const wn=[void 0,1,1,2,4,8,1,1,2,4,8,4,8,4];class Sn extends B{parseHeader(){var e=this.chunk.getUint16();e===18761?this.le=!0:e===19789&&(this.le=!1),this.chunk.le=this.le,this.headerParsed=!0}parseTags(e,t,n=new Map){let{pick:i,skip:s}=this.options[t];i=new Set(i);let a=i.size>0,o=s.size===0,c=this.chunk.getUint16(e);e+=2;for(let u=0;u<c;u++){let d=this.chunk.getUint16(e);if(a){if(i.has(d)&&(n.set(d,this.parseTag(e,d,t)),i.delete(d),i.size===0))break}else!o&&s.has(d)||n.set(d,this.parseTag(e,d,t));e+=12}return n}parseTag(e,t,n){let{chunk:i}=this,s=i.getUint16(e+2),a=i.getUint32(e+4),o=wn[s];if(o*a<=4?e+=8:e=i.getUint32(e+8),(s<1||s>13)&&I(`Invalid TIFF value type. block: ${n.toUpperCase()}, tag: ${t.toString(16)}, type: ${s}, offset ${e}`),e>i.byteLength&&I(`Invalid TIFF value offset. block: ${n.toUpperCase()}, tag: ${t.toString(16)}, type: ${s}, offset ${e} is outside of chunk size ${i.byteLength}`),s===1)return i.getUint8Array(e,a);if(s===2)return Z(i.getString(e,a));if(s===7)return i.getUint8Array(e,a);if(a===1)return this.parseTagValue(s,e);{let c=new(function(d){switch(d){case 1:return Uint8Array;case 3:return Uint16Array;case 4:return Uint32Array;case 5:return Array;case 6:return Int8Array;case 8:return Int16Array;case 9:return Int32Array;case 10:return Array;case 11:return Float32Array;case 12:return Float64Array;default:return Array}}(s))(a),u=o;for(let d=0;d<a;d++)c[d]=this.parseTagValue(s,e),e+=u;return c}}parseTagValue(e,t){let{chunk:n}=this;switch(e){case 1:return n.getUint8(t);case 3:return n.getUint16(t);case 4:return n.getUint32(t);case 5:return n.getUint32(t)/n.getUint32(t+4);case 6:return n.getInt8(t);case 8:return n.getInt16(t);case 9:return n.getInt32(t);case 10:return n.getInt32(t)/n.getInt32(t+4);case 11:return n.getFloat(t);case 12:return n.getDouble(t);case 13:return n.getUint32(t);default:I(`Invalid tiff type ${e}`)}}}class Ve extends Sn{static canHandle(e,t){return e.getUint8(t+1)===225&&e.getUint32(t+4)===1165519206&&e.getUint16(t+8)===0}async parse(){this.parseHeader();let{options:e}=this;return e.ifd0.enabled&&await this.parseIfd0Block(),e.exif.enabled&&await this.safeParse("parseExifBlock"),e.gps.enabled&&await this.safeParse("parseGpsBlock"),e.interop.enabled&&await this.safeParse("parseInteropBlock"),e.ifd1.enabled&&await this.safeParse("parseThumbnailBlock"),this.createOutput()}safeParse(e){let t=this[e]();return t.catch!==void 0&&(t=t.catch(this.handleError)),t}findIfd0Offset(){this.ifd0Offset===void 0&&(this.ifd0Offset=this.chunk.getUint32(4))}findIfd1Offset(){if(this.ifd1Offset===void 0){this.findIfd0Offset();let e=this.chunk.getUint16(this.ifd0Offset),t=this.ifd0Offset+2+12*e;this.ifd1Offset=this.chunk.getUint32(t)}}parseBlock(e,t){let n=new Map;return this[t]=n,this.parseTags(e,t,n),n}async parseIfd0Block(){if(this.ifd0)return;let{file:e}=this;this.findIfd0Offset(),this.ifd0Offset<8&&I("Malformed EXIF data"),!e.chunked&&this.ifd0Offset>e.byteLength&&I(`IFD0 offset points to outside of file.
this.ifd0Offset: ${this.ifd0Offset}, file.byteLength: ${e.byteLength}`),e.tiff&&await e.ensureChunk(this.ifd0Offset,Ke(this.options));let t=this.parseBlock(this.ifd0Offset,"ifd0");return t.size!==0?(this.exifOffset=t.get(34665),this.interopOffset=t.get(40965),this.gpsOffset=t.get(34853),this.xmp=t.get(700),this.iptc=t.get(33723),this.icc=t.get(34675),this.options.sanitize&&(t.delete(34665),t.delete(40965),t.delete(34853),t.delete(700),t.delete(33723),t.delete(34675)),t):void 0}async parseExifBlock(){if(this.exif||(this.ifd0||await this.parseIfd0Block(),this.exifOffset===void 0))return;this.file.tiff&&await this.file.ensureChunk(this.exifOffset,Ke(this.options));let e=this.parseBlock(this.exifOffset,"exif");return this.interopOffset||(this.interopOffset=e.get(40965)),this.makerNote=e.get(37500),this.userComment=e.get(37510),this.options.sanitize&&(e.delete(40965),e.delete(37500),e.delete(37510)),this.unpack(e,41728),this.unpack(e,41729),e}unpack(e,t){let n=e.get(t);n&&n.length===1&&e.set(t,n[0])}async parseGpsBlock(){if(this.gps||(this.ifd0||await this.parseIfd0Block(),this.gpsOffset===void 0))return;let e=this.parseBlock(this.gpsOffset,"gps");return e&&e.has(2)&&e.has(4)&&(e.set("latitude",pt(...e.get(2),e.get(1))),e.set("longitude",pt(...e.get(4),e.get(3)))),e}async parseInteropBlock(){if(!this.interop&&(this.ifd0||await this.parseIfd0Block(),this.interopOffset!==void 0||this.exif||await this.parseExifBlock(),this.interopOffset!==void 0))return this.parseBlock(this.interopOffset,"interop")}async parseThumbnailBlock(e=!1){if(!this.ifd1&&!this.ifd1Parsed&&(!this.options.mergeOutput||e))return this.findIfd1Offset(),this.ifd1Offset>0&&(this.parseBlock(this.ifd1Offset,"ifd1"),this.ifd1Parsed=!0),this.ifd1}async extractThumbnail(){if(this.headerParsed||this.parseHeader(),this.ifd1Parsed||await this.parseThumbnailBlock(!0),this.ifd1===void 0)return;let e=this.ifd1.get(513),t=this.ifd1.get(514);return this.chunk.getUint8Array(e,t)}get image(){return this.ifd0}get thumbnail(){return this.ifd1}createOutput(){let e,t,n,i={};for(t of x)if(e=this[t],!Ot(e))if(n=this.canTranslate?this.translateBlock(e,t):Object.fromEntries(e),this.options.mergeOutput){if(t==="ifd1")continue;Object.assign(i,n)}else i[t]=n;return this.makerNote&&(i.makerNote=this.makerNote),this.userComment&&(i.userComment=this.userComment),i}assignToOutput(e,t){if(this.globalOptions.mergeOutput)Object.assign(e,t);else for(let[n,i]of Object.entries(t))this.assignObjectToOutput(e,n,i)}}function pt(r,e,t,n){var i=r+e/60+t/3600;return n!=="S"&&n!=="W"||(i*=-1),i}g(Ve,"type","tiff"),g(Ve,"headerLength",10),P.set("tiff",Ve);var kn=Object.freeze({__proto__:null,default:mn,Exifr:$,fileParsers:G,segmentParsers:P,fileReaders:H,tagKeys:T,tagValues:z,tagRevivers:J,createDictionary:v,extendDictionary:re,fetchUrlAsArrayBuffer:ne,readBlobAsArrayBuffer:ie,chunkedProps:X,otherSegments:he,segments:se,tiffBlocks:x,segmentsAndBlocks:K,tiffExtractables:q,inheritables:ue,allFormatters:Y,Options:ae,parse:Le});const et={ifd0:!1,ifd1:!1,exif:!1,gps:!1,interop:!1,sanitize:!1,reviveValues:!0,translateKeys:!1,translateValues:!1,mergeOutput:!1},tt=Object.assign({},et,{firstChunkSize:4e4,gps:[1,2,3,4]});async function Ft(r){let e=new $(tt);await e.read(r);let t=await e.parse();if(t&&t.gps){let{latitude:n,longitude:i}=t.gps;return{latitude:n,longitude:i}}}const nt=Object.assign({},et,{tiff:!1,ifd1:!0,mergeOutput:!1});async function Lt(r){let e=new $(nt);await e.read(r);let t=await e.extractThumbnail();return t&&Fe?Me.from(t):t}async function Et(r){let e=await this.thumbnail(r);if(e!==void 0){let t=new Blob([e]);return URL.createObjectURL(t)}}const it=Object.assign({},et,{firstChunkSize:4e4,ifd0:[274]});async function rt(r){let e=new $(it);await e.read(r);let t=await e.parse();if(t&&t.ifd0)return t.ifd0[274]}const st=Object.freeze({1:{dimensionSwapped:!1,scaleX:1,scaleY:1,deg:0,rad:0},2:{dimensionSwapped:!1,scaleX:-1,scaleY:1,deg:0,rad:0},3:{dimensionSwapped:!1,scaleX:1,scaleY:1,deg:180,rad:180*Math.PI/180},4:{dimensionSwapped:!1,scaleX:-1,scaleY:1,deg:180,rad:180*Math.PI/180},5:{dimensionSwapped:!0,scaleX:1,scaleY:-1,deg:90,rad:90*Math.PI/180},6:{dimensionSwapped:!0,scaleX:1,scaleY:1,deg:90,rad:90*Math.PI/180},7:{dimensionSwapped:!0,scaleX:1,scaleY:-1,deg:270,rad:270*Math.PI/180},8:{dimensionSwapped:!0,scaleX:1,scaleY:1,deg:270,rad:270*Math.PI/180}});let ee=!0,te=!0;if(typeof navigator=="object"){let r=navigator.userAgent;if(r.includes("iPad")||r.includes("iPhone")){let e=r.match(/OS (\d+)_(\d+)/);if(e){let[,t,n]=e;ee=Number(t)+.1*Number(n)<13.4,te=!1}}else if(r.includes("OS X 10")){let[,e]=r.match(/OS X 10[_.](\d+)/);ee=te=Number(e)<15}if(r.includes("Chrome/")){let[,e]=r.match(/Chrome\/(\d+)/);ee=te=Number(e)<81}else if(r.includes("Firefox/")){let[,e]=r.match(/Firefox\/(\d+)/);ee=te=Number(e)<77}}async function Ut(r){let e=await rt(r);return Object.assign({canvas:ee,css:te},st[e])}class Cn extends L{constructor(...e){super(...e),g(this,"ranges",new xn),this.byteLength!==0&&this.ranges.add(0,this.byteLength)}_tryExtend(e,t,n){if(e===0&&this.byteLength===0&&n){let i=new DataView(n.buffer||n,n.byteOffset,n.byteLength);this._swapDataView(i)}else{let i=e+t;if(i>this.byteLength){let{dataView:s}=this._extend(i);this._swapDataView(s)}}}_extend(e){let t;t=Fe?Me.allocUnsafe(e):new Uint8Array(e);let n=new DataView(t.buffer,t.byteOffset,t.byteLength);return t.set(new Uint8Array(this.buffer,this.byteOffset,this.byteLength),0),{uintView:t,dataView:n}}subarray(e,t,n=!1){return t=t||this._lengthToEnd(e),n&&this._tryExtend(e,t),this.ranges.add(e,t),super.subarray(e,t)}set(e,t,n=!1){n&&this._tryExtend(t,e.byteLength,e);let i=super.set(e,t);return this.ranges.add(t,i.byteLength),i}async ensureChunk(e,t){this.chunked&&(this.ranges.available(e,t)||await this.readChunk(e,t))}available(e,t){return this.ranges.available(e,t)}}class xn{constructor(){g(this,"list",[])}get length(){return this.list.length}add(e,t,n=0){let i=e+t,s=this.list.filter(a=>gt(e,a.offset,i)||gt(e,a.end,i));if(s.length>0){e=Math.min(e,...s.map(o=>o.offset)),i=Math.max(i,...s.map(o=>o.end)),t=i-e;let a=s.shift();a.offset=e,a.length=t,a.end=i,this.list=this.list.filter(o=>!s.includes(o))}else this.list.push({offset:e,length:t,end:i})}available(e,t){let n=e+t;return this.list.some(i=>i.offset<=e&&n<=i.end)}}function gt(r,e,t){return r<=e&&e<=t}class Ue extends Cn{constructor(e,t){super(0),g(this,"chunksRead",0),this.input=e,this.options=t}async readWhole(){this.chunked=!1,await this.readChunk(this.nextChunkOffset)}async readChunked(){this.chunked=!0,await this.readChunk(0,this.options.firstChunkSize)}async readNextChunk(e=this.nextChunkOffset){if(this.fullyRead)return this.chunksRead++,!1;let t=this.options.chunkSize,n=await this.readChunk(e,t);return!!n&&n.byteLength===t}async readChunk(e,t){if(this.chunksRead++,(t=this.safeWrapAddress(e,t))!==0)return this._readChunk(e,t)}safeWrapAddress(e,t){return this.size!==void 0&&e+t>this.size?Math.max(0,this.size-e):t}get nextChunkOffset(){if(this.ranges.list.length!==0)return this.ranges.list[0].length}get canReadNextChunk(){return this.chunksRead<this.options.chunkLimit}get fullyRead(){return this.size!==void 0&&this.nextChunkOffset===this.size}read(){return this.options.chunked?this.readChunked():this.readWhole()}close(){}}H.set("blob",class extends Ue{async readWhole(){this.chunked=!1;let r=await ie(this.input);this._swapArrayBuffer(r)}readChunked(){return this.chunked=!0,this.size=this.input.size,super.readChunked()}async _readChunk(r,e){let t=e?r+e:void 0,n=this.input.slice(r,t),i=await ie(n);return this.set(i,r,!0)}});var vn=Object.freeze({__proto__:null,default:kn,Exifr:$,fileParsers:G,segmentParsers:P,fileReaders:H,tagKeys:T,tagValues:z,tagRevivers:J,createDictionary:v,extendDictionary:re,fetchUrlAsArrayBuffer:ne,readBlobAsArrayBuffer:ie,chunkedProps:X,otherSegments:he,segments:se,tiffBlocks:x,segmentsAndBlocks:K,tiffExtractables:q,inheritables:ue,allFormatters:Y,Options:ae,parse:Le,gpsOnlyOptions:tt,gps:Ft,thumbnailOnlyOptions:nt,thumbnail:Lt,thumbnailUrl:Et,orientationOnlyOptions:it,orientation:rt,rotations:st,get rotateCanvas(){return ee},get rotateCss(){return te},rotation:Ut});H.set("url",class extends Ue{async readWhole(){this.chunked=!1;let r=await ne(this.input);r instanceof ArrayBuffer?this._swapArrayBuffer(r):r instanceof Uint8Array&&this._swapBuffer(r)}async _readChunk(r,e){let t=e?r+e-1:void 0,n=this.options.httpHeaders||{};(r||t)&&(n.range=`bytes=${[r,t].join("-")}`);let i=await Qe(this.input,{headers:n}),s=await i.arrayBuffer(),a=s.byteLength;if(i.status!==416)return a!==e&&(this.size=r+a),this.set(s,r,!0)}});L.prototype.getUint64=function(r){let e=this.getUint32(r),t=this.getUint32(r+4);return e<1048575?e<<32|t:typeof we!==void 0?(console.warn("Using BigInt because of type 64uint but JS can only handle 53b numbers."),we(e)<<we(32)|we(t)):void I("Trying to read 64b value but JS can only handle 53b numbers.")};class Pn extends Ee{parseBoxes(e=0){let t=[];for(;e<this.file.byteLength-4;){let n=this.parseBoxHead(e);if(t.push(n),n.length===0)break;e+=n.length}return t}parseSubBoxes(e){e.boxes=this.parseBoxes(e.start)}findBox(e,t){return e.boxes===void 0&&this.parseSubBoxes(e),e.boxes.find(n=>n.kind===t)}parseBoxHead(e){let t=this.file.getUint32(e),n=this.file.getString(e+4,4),i=e+8;return t===1&&(t=this.file.getUint64(e+8),i+=8),{offset:e,length:t,kind:n,start:i}}parseBoxFullHead(e){if(e.version!==void 0)return;let t=this.file.getUint32(e.start);e.version=t>>24,e.start+=4}}class Nt extends Pn{static canHandle(e,t){if(t!==0)return!1;let n=e.getUint16(2);if(n>50)return!1;let i=16,s=[];for(;i<n;)s.push(e.getString(i,4)),i+=4;return s.includes(this.type)}async parse(){let e=this.file.getUint32(0),t=this.parseBoxHead(e);for(;t.kind!=="meta";)e+=t.length,await this.file.ensureChunk(e,16),t=this.parseBoxHead(e);await this.file.ensureChunk(t.offset,t.length),this.parseBoxFullHead(t),this.parseSubBoxes(t),this.options.icc.enabled&&await this.findIcc(t),this.options.tiff.enabled&&await this.findExif(t)}async registerSegment(e,t,n){await this.file.ensureChunk(t,n);let i=this.file.subarray(t,n);this.createParser(e,i)}async findIcc(e){let t=this.findBox(e,"iprp");if(t===void 0)return;let n=this.findBox(t,"ipco");if(n===void 0)return;let i=this.findBox(n,"colr");i!==void 0&&await this.registerSegment("icc",i.offset+12,i.length)}async findExif(e){let t=this.findBox(e,"iinf");if(t===void 0)return;let n=this.findBox(e,"iloc");if(n===void 0)return;let i=this.findExifLocIdInIinf(t),s=this.findExtentInIloc(n,i);if(s===void 0)return;let[a,o]=s;await this.file.ensureChunk(a,o);let c=4+this.file.getUint32(a);a+=c,o-=c,await this.registerSegment("tiff",a,o)}findExifLocIdInIinf(e){this.parseBoxFullHead(e);let t,n,i,s,a=e.start,o=this.file.getUint16(a);for(a+=2;o--;){if(t=this.parseBoxHead(a),this.parseBoxFullHead(t),n=t.start,t.version>=2&&(i=t.version===3?4:2,s=this.file.getString(n+i+2,4),s==="Exif"))return this.file.getUintBytes(n,i);a+=t.length}}get8bits(e){let t=this.file.getUint8(e);return[t>>4,15&t]}findExtentInIloc(e,t){this.parseBoxFullHead(e);let n=e.start,[i,s]=this.get8bits(n++),[a,o]=this.get8bits(n++),c=e.version===2?4:2,u=e.version===1||e.version===2?2:0,d=o+i+s,l=e.version===2?4:2,f=this.file.getUintBytes(n,l);for(n+=l;f--;){let h=this.file.getUintBytes(n,c);n+=c+u+2+a;let p=this.file.getUint16(n);if(n+=2,h===t)return p>1&&console.warn(`ILOC box has more than one extent but we're only processing one
Please create an issue at https://github.com/MikeKovarik/exifr with this file`),[this.file.getUintBytes(n+o,i),this.file.getUintBytes(n+o+i,s)];n+=p*d}}}class Bt extends Nt{}g(Bt,"type","heic");class mt extends Nt{}g(mt,"type","avif"),G.set("heic",Bt),G.set("avif",mt),v(T,["ifd0","ifd1"],[[256,"ImageWidth"],[257,"ImageHeight"],[258,"BitsPerSample"],[259,"Compression"],[262,"PhotometricInterpretation"],[270,"ImageDescription"],[271,"Make"],[272,"Model"],[273,"StripOffsets"],[274,"Orientation"],[277,"SamplesPerPixel"],[278,"RowsPerStrip"],[279,"StripByteCounts"],[282,"XResolution"],[283,"YResolution"],[284,"PlanarConfiguration"],[296,"ResolutionUnit"],[301,"TransferFunction"],[305,"Software"],[306,"ModifyDate"],[315,"Artist"],[316,"HostComputer"],[317,"Predictor"],[318,"WhitePoint"],[319,"PrimaryChromaticities"],[513,"ThumbnailOffset"],[514,"ThumbnailLength"],[529,"YCbCrCoefficients"],[530,"YCbCrSubSampling"],[531,"YCbCrPositioning"],[532,"ReferenceBlackWhite"],[700,"ApplicationNotes"],[33432,"Copyright"],[33723,"IPTC"],[34665,"ExifIFD"],[34675,"ICC"],[34853,"GpsIFD"],[330,"SubIFD"],[40965,"InteropIFD"],[40091,"XPTitle"],[40092,"XPComment"],[40093,"XPAuthor"],[40094,"XPKeywords"],[40095,"XPSubject"]]),v(T,"exif",[[33434,"ExposureTime"],[33437,"FNumber"],[34850,"ExposureProgram"],[34852,"SpectralSensitivity"],[34855,"ISO"],[34858,"TimeZoneOffset"],[34859,"SelfTimerMode"],[34864,"SensitivityType"],[34865,"StandardOutputSensitivity"],[34866,"RecommendedExposureIndex"],[34867,"ISOSpeed"],[34868,"ISOSpeedLatitudeyyy"],[34869,"ISOSpeedLatitudezzz"],[36864,"ExifVersion"],[36867,"DateTimeOriginal"],[36868,"CreateDate"],[36873,"GooglePlusUploadCode"],[36880,"OffsetTime"],[36881,"OffsetTimeOriginal"],[36882,"OffsetTimeDigitized"],[37121,"ComponentsConfiguration"],[37122,"CompressedBitsPerPixel"],[37377,"ShutterSpeedValue"],[37378,"ApertureValue"],[37379,"BrightnessValue"],[37380,"ExposureCompensation"],[37381,"MaxApertureValue"],[37382,"SubjectDistance"],[37383,"MeteringMode"],[37384,"LightSource"],[37385,"Flash"],[37386,"FocalLength"],[37393,"ImageNumber"],[37394,"SecurityClassification"],[37395,"ImageHistory"],[37396,"SubjectArea"],[37500,"MakerNote"],[37510,"UserComment"],[37520,"SubSecTime"],[37521,"SubSecTimeOriginal"],[37522,"SubSecTimeDigitized"],[37888,"AmbientTemperature"],[37889,"Humidity"],[37890,"Pressure"],[37891,"WaterDepth"],[37892,"Acceleration"],[37893,"CameraElevationAngle"],[40960,"FlashpixVersion"],[40961,"ColorSpace"],[40962,"ExifImageWidth"],[40963,"ExifImageHeight"],[40964,"RelatedSoundFile"],[41483,"FlashEnergy"],[41486,"FocalPlaneXResolution"],[41487,"FocalPlaneYResolution"],[41488,"FocalPlaneResolutionUnit"],[41492,"SubjectLocation"],[41493,"ExposureIndex"],[41495,"SensingMethod"],[41728,"FileSource"],[41729,"SceneType"],[41730,"CFAPattern"],[41985,"CustomRendered"],[41986,"ExposureMode"],[41987,"WhiteBalance"],[41988,"DigitalZoomRatio"],[41989,"FocalLengthIn35mmFormat"],[41990,"SceneCaptureType"],[41991,"GainControl"],[41992,"Contrast"],[41993,"Saturation"],[41994,"Sharpness"],[41996,"SubjectDistanceRange"],[42016,"ImageUniqueID"],[42032,"OwnerName"],[42033,"SerialNumber"],[42034,"LensInfo"],[42035,"LensMake"],[42036,"LensModel"],[42037,"LensSerialNumber"],[42080,"CompositeImage"],[42081,"CompositeImageCount"],[42082,"CompositeImageExposureTimes"],[42240,"Gamma"],[59932,"Padding"],[59933,"OffsetSchema"],[65e3,"OwnerName"],[65001,"SerialNumber"],[65002,"Lens"],[65100,"RawFile"],[65101,"Converter"],[65102,"WhiteBalance"],[65105,"Exposure"],[65106,"Shadows"],[65107,"Brightness"],[65108,"Contrast"],[65109,"Saturation"],[65110,"Sharpness"],[65111,"Smoothness"],[65112,"MoireFilter"],[40965,"InteropIFD"]]),v(T,"gps",[[0,"GPSVersionID"],[1,"GPSLatitudeRef"],[2,"GPSLatitude"],[3,"GPSLongitudeRef"],[4,"GPSLongitude"],[5,"GPSAltitudeRef"],[6,"GPSAltitude"],[7,"GPSTimeStamp"],[8,"GPSSatellites"],[9,"GPSStatus"],[10,"GPSMeasureMode"],[11,"GPSDOP"],[12,"GPSSpeedRef"],[13,"GPSSpeed"],[14,"GPSTrackRef"],[15,"GPSTrack"],[16,"GPSImgDirectionRef"],[17,"GPSImgDirection"],[18,"GPSMapDatum"],[19,"GPSDestLatitudeRef"],[20,"GPSDestLatitude"],[21,"GPSDestLongitudeRef"],[22,"GPSDestLongitude"],[23,"GPSDestBearingRef"],[24,"GPSDestBearing"],[25,"GPSDestDistanceRef"],[26,"GPSDestDistance"],[27,"GPSProcessingMethod"],[28,"GPSAreaInformation"],[29,"GPSDateStamp"],[30,"GPSDifferential"],[31,"GPSHPositioningError"]]),v(z,["ifd0","ifd1"],[[274,{1:"Horizontal (normal)",2:"Mirror horizontal",3:"Rotate 180",4:"Mirror vertical",5:"Mirror horizontal and rotate 270 CW",6:"Rotate 90 CW",7:"Mirror horizontal and rotate 90 CW",8:"Rotate 270 CW"}],[296,{1:"None",2:"inches",3:"cm"}]]);let me=v(z,"exif",[[34850,{0:"Not defined",1:"Manual",2:"Normal program",3:"Aperture priority",4:"Shutter priority",5:"Creative program",6:"Action program",7:"Portrait mode",8:"Landscape mode"}],[37121,{0:"-",1:"Y",2:"Cb",3:"Cr",4:"R",5:"G",6:"B"}],[37383,{0:"Unknown",1:"Average",2:"CenterWeightedAverage",3:"Spot",4:"MultiSpot",5:"Pattern",6:"Partial",255:"Other"}],[37384,{0:"Unknown",1:"Daylight",2:"Fluorescent",3:"Tungsten (incandescent light)",4:"Flash",9:"Fine weather",10:"Cloudy weather",11:"Shade",12:"Daylight fluorescent (D 5700 - 7100K)",13:"Day white fluorescent (N 4600 - 5400K)",14:"Cool white fluorescent (W 3900 - 4500K)",15:"White fluorescent (WW 3200 - 3700K)",17:"Standard light A",18:"Standard light B",19:"Standard light C",20:"D55",21:"D65",22:"D75",23:"D50",24:"ISO studio tungsten",255:"Other"}],[37385,{0:"Flash did not fire",1:"Flash fired",5:"Strobe return light not detected",7:"Strobe return light detected",9:"Flash fired, compulsory flash mode",13:"Flash fired, compulsory flash mode, return light not detected",15:"Flash fired, compulsory flash mode, return light detected",16:"Flash did not fire, compulsory flash mode",24:"Flash did not fire, auto mode",25:"Flash fired, auto mode",29:"Flash fired, auto mode, return light not detected",31:"Flash fired, auto mode, return light detected",32:"No flash function",65:"Flash fired, red-eye reduction mode",69:"Flash fired, red-eye reduction mode, return light not detected",71:"Flash fired, red-eye reduction mode, return light detected",73:"Flash fired, compulsory flash mode, red-eye reduction mode",77:"Flash fired, compulsory flash mode, red-eye reduction mode, return light not detected",79:"Flash fired, compulsory flash mode, red-eye reduction mode, return light detected",89:"Flash fired, auto mode, red-eye reduction mode",93:"Flash fired, auto mode, return light not detected, red-eye reduction mode",95:"Flash fired, auto mode, return light detected, red-eye reduction mode"}],[41495,{1:"Not defined",2:"One-chip color area sensor",3:"Two-chip color area sensor",4:"Three-chip color area sensor",5:"Color sequential area sensor",7:"Trilinear sensor",8:"Color sequential linear sensor"}],[41728,{1:"Film Scanner",2:"Reflection Print Scanner",3:"Digital Camera"}],[41729,{1:"Directly photographed"}],[41985,{0:"Normal",1:"Custom",2:"HDR (no original saved)",3:"HDR (original saved)",4:"Original (for HDR)",6:"Panorama",7:"Portrait HDR",8:"Portrait"}],[41986,{0:"Auto",1:"Manual",2:"Auto bracket"}],[41987,{0:"Auto",1:"Manual"}],[41990,{0:"Standard",1:"Landscape",2:"Portrait",3:"Night",4:"Other"}],[41991,{0:"None",1:"Low gain up",2:"High gain up",3:"Low gain down",4:"High gain down"}],[41996,{0:"Unknown",1:"Macro",2:"Close",3:"Distant"}],[42080,{0:"Unknown",1:"Not a Composite Image",2:"General Composite Image",3:"Composite Image Captured While Shooting"}]]);const _t={1:"No absolute unit of measurement",2:"Inch",3:"Centimeter"};me.set(37392,_t),me.set(41488,_t);const Ge={0:"Normal",1:"Low",2:"High"};function yt(r){return typeof r=="object"&&r.length!==void 0?r[0]:r}function bt(r){let e=Array.from(r).slice(1);return e[1]>15&&(e=e.map(t=>String.fromCharCode(t))),e[2]!=="0"&&e[2]!==0||e.pop(),e.join(".")}function He(r){if(typeof r=="string"){var[e,t,n,i,s,a]=r.trim().split(/[-: ]/g).map(Number),o=new Date(e,t-1,n);return Number.isNaN(i)||Number.isNaN(s)||Number.isNaN(a)||(o.setHours(i),o.setMinutes(s),o.setSeconds(a)),Number.isNaN(+o)?r:o}}function ge(r){if(typeof r=="string")return r;let e=[];if(r[1]===0&&r[r.length-1]===0)for(let t=0;t<r.length;t+=2)e.push(wt(r[t+1],r[t]));else for(let t=0;t<r.length;t+=2)e.push(wt(r[t],r[t+1]));return Z(String.fromCodePoint(...e))}function wt(r,e){return r<<8|e}me.set(41992,Ge),me.set(41993,Ge),me.set(41994,Ge),v(J,["ifd0","ifd1"],[[50827,function(r){return typeof r!="string"?Dt(r):r}],[306,He],[40091,ge],[40092,ge],[40093,ge],[40094,ge],[40095,ge]]),v(J,"exif",[[40960,bt],[36864,bt],[36867,He],[36868,He],[40962,yt],[40963,yt]]),v(J,"gps",[[0,r=>Array.from(r).join(".")],[7,r=>Array.from(r).join(":")]]);class je extends B{static canHandle(e,t){return e.getUint8(t+1)===225&&e.getUint32(t+4)===1752462448&&e.getString(t+4,20)==="http://ns.adobe.com/"}static headerLength(e,t){return e.getString(t+4,34)==="http://ns.adobe.com/xmp/extension/"?79:33}static findPosition(e,t){let n=super.findPosition(e,t);return n.multiSegment=n.extended=n.headerLength===79,n.multiSegment?(n.chunkCount=e.getUint8(t+72),n.chunkNumber=e.getUint8(t+76),e.getUint8(t+77)!==0&&n.chunkNumber++):(n.chunkCount=1/0,n.chunkNumber=-1),n}static handleMultiSegments(e){return e.map(t=>t.chunk.getString()).join("")}normalizeInput(e){return typeof e=="string"?e:L.from(e).getString()}parse(e=this.chunk){if(!this.localOptions.parse)return e;e=function(s){let a={},o={};for(let c of Ht)a[c]=[],o[c]=0;return s.replace(Rn,(c,u,d)=>{if(u==="<"){let l=++o[d];return a[d].push(l),`${c}#${l}`}return`${c}#${a[d].pop()}`})}(e);let t=de.findAll(e,"rdf","Description");t.length===0&&t.push(new de("rdf","Description",void 0,e));let n,i={};for(let s of t)for(let a of s.properties)n=Tn(a.ns,i),zt(a,n);return function(s){let a;for(let o in s)a=s[o]=Ae(s[o]),a===void 0&&delete s[o];return Ae(s)}(i)}assignToOutput(e,t){if(this.localOptions.parse)for(let[n,i]of Object.entries(t))switch(n){case"tiff":this.assignObjectToOutput(e,"ifd0",i);break;case"exif":this.assignObjectToOutput(e,"exif",i);break;case"xmlns":break;default:this.assignObjectToOutput(e,n,i)}else e.xmp=t}}g(je,"type","xmp"),g(je,"multiSegment",!0),P.set("xmp",je);class Re{static findAll(e){return Vt(e,/([a-zA-Z0-9-]+):([a-zA-Z0-9-]+)=("[^"]*"|'[^']*')/gm).map(Re.unpackMatch)}static unpackMatch(e){let t=e[1],n=e[2],i=e[3].slice(1,-1);return i=Gt(i),new Re(t,n,i)}constructor(e,t,n){this.ns=e,this.name=t,this.value=n}serialize(){return this.value}}class de{static findAll(e,t,n){if(t!==void 0||n!==void 0){t=t||"[\\w\\d-]+",n=n||"[\\w\\d-]+";var i=new RegExp(`<(${t}):(${n})(#\\d+)?((\\s+?[\\w\\d-:]+=("[^"]*"|'[^']*'))*\\s*)(\\/>|>([\\s\\S]*?)<\\/\\1:\\2\\3>)`,"gm")}else i=/<([\w\d-]+):([\w\d-]+)(#\d+)?((\s+?[\w\d-:]+=("[^"]*"|'[^']*'))*\s*)(\/>|>([\s\S]*?)<\/\1:\2\3>)/gm;return Vt(e,i).map(de.unpackMatch)}static unpackMatch(e){let t=e[1],n=e[2],i=e[4],s=e[8];return new de(t,n,i,s)}constructor(e,t,n,i){this.ns=e,this.name=t,this.attrString=n,this.innerXml=i,this.attrs=Re.findAll(n),this.children=de.findAll(i),this.value=this.children.length===0?Gt(i):void 0,this.properties=[...this.attrs,...this.children]}get isPrimitive(){return this.value!==void 0&&this.attrs.length===0&&this.children.length===0}get isListContainer(){return this.children.length===1&&this.children[0].isList}get isList(){let{ns:e,name:t}=this;return e==="rdf"&&(t==="Seq"||t==="Bag"||t==="Alt")}get isListItem(){return this.ns==="rdf"&&this.name==="li"}serialize(){if(this.properties.length===0&&this.value===void 0)return;if(this.isPrimitive)return this.value;if(this.isListContainer)return this.children[0].serialize();if(this.isList)return An(this.children.map(In));if(this.isListItem&&this.children.length===1&&this.attrs.length===0)return this.children[0].serialize();let e={};for(let t of this.properties)zt(t,e);return this.value!==void 0&&(e.value=this.value),Ae(e)}}function zt(r,e){let t=r.serialize();t!==void 0&&(e[r.name]=t)}var In=r=>r.serialize(),An=r=>r.length===1?r[0]:r,Tn=(r,e)=>e[r]?e[r]:e[r]={};function Vt(r,e){let t,n=[];if(!r)return n;for(;(t=e.exec(r))!==null;)n.push(t);return n}function Gt(r){if(function(n){return n==null||n==="null"||n==="undefined"||n===""||n.trim()===""}(r))return;let e=Number(r);if(!Number.isNaN(e))return e;let t=r.toLowerCase();return t==="true"||t!=="false"&&r.trim()}const Ht=["rdf:li","rdf:Seq","rdf:Bag","rdf:Alt","rdf:Description"],Rn=new RegExp(`(<|\\/)(${Ht.join("|")})`,"g");var On=Object.freeze({__proto__:null,default:vn,Exifr:$,fileParsers:G,segmentParsers:P,fileReaders:H,tagKeys:T,tagValues:z,tagRevivers:J,createDictionary:v,extendDictionary:re,fetchUrlAsArrayBuffer:ne,readBlobAsArrayBuffer:ie,chunkedProps:X,otherSegments:he,segments:se,tiffBlocks:x,segmentsAndBlocks:K,tiffExtractables:q,inheritables:ue,allFormatters:Y,Options:ae,parse:Le,gpsOnlyOptions:tt,gps:Ft,thumbnailOnlyOptions:nt,thumbnail:Lt,thumbnailUrl:Et,orientationOnlyOptions:it,orientation:rt,rotations:st,get rotateCanvas(){return ee},get rotateCss(){return te},rotation:Ut});let St=Ie("fs",r=>r.promises);H.set("fs",class extends Ue{async readWhole(){this.chunked=!1,this.fs=await St;let r=await this.fs.readFile(this.input);this._swapBuffer(r)}async readChunked(){this.chunked=!0,this.fs=await St,await this.open(),await this.readChunk(0,this.options.firstChunkSize)}async open(){this.fh===void 0&&(this.fh=await this.fs.open(this.input,"r"),this.size=(await this.fh.stat(this.input)).size)}async _readChunk(r,e){this.fh===void 0&&await this.open(),r+e>this.size&&(e=this.size-r);var t=this.subarray(r,e,!0);return await this.fh.read(t.dataView,0,e,r),t}async close(){if(this.fh){let r=this.fh;this.fh=void 0,await r.close()}}});H.set("base64",class extends Ue{constructor(...r){super(...r),this.input=this.input.replace(/^data:([^;]+);base64,/gim,""),this.size=this.input.length/4*3,this.input.endsWith("==")?this.size-=2:this.input.endsWith("=")&&(this.size-=1)}async _readChunk(r,e){let t,n,i=this.input;r===void 0?(r=0,t=0,n=0):(t=4*Math.floor(r/3),n=r-t/4*3),e===void 0&&(e=this.size);let s=r+e,a=t+4*Math.ceil(s/3);i=i.slice(t,a);let o=Math.min(e,this.size-r);if(Fe){let c=Me.from(i,"base64").slice(n,n+o);return this.set(c,r,!0)}{let c=this.subarray(r,o,!0),u=atob(i),d=c.toUint8();for(let l=0;l<o;l++)d[l]=u.charCodeAt(n+l);return c}}});class kt extends Ee{static canHandle(e,t){return t===18761||t===19789}extendOptions(e){let{ifd0:t,xmp:n,iptc:i,icc:s}=e;n.enabled&&t.deps.add(700),i.enabled&&t.deps.add(33723),s.enabled&&t.deps.add(34675),t.finalizeFilters()}async parse(){let{tiff:e,xmp:t,iptc:n,icc:i}=this.options;if(e.enabled||t.enabled||n.enabled||i.enabled){let s=Math.max(Ke(this.options),this.options.chunkSize);await this.file.ensureChunk(0,s),this.createParser("tiff",this.file),this.parsers.tiff.parseHeader(),await this.parsers.tiff.parseIfd0Block(),this.adaptTiffPropAsSegment("xmp"),this.adaptTiffPropAsSegment("iptc"),this.adaptTiffPropAsSegment("icc")}}adaptTiffPropAsSegment(e){if(this.parsers.tiff[e]){let t=this.parsers.tiff[e];this.injectSegment(e,t)}}}g(kt,"type","tiff"),G.set("tiff",kt);let Dn=Ie("zlib");const Mn=["ihdr","iccp","text","itxt","exif"];class Ct extends Ee{constructor(...e){super(...e),g(this,"catchError",t=>this.errors.push(t)),g(this,"metaChunks",[]),g(this,"unknownChunks",[])}static canHandle(e,t){return t===35152&&e.getUint32(0)===2303741511&&e.getUint32(4)===218765834}async parse(){let{file:e}=this;await this.findPngChunksInRange(8,e.byteLength),await this.readSegments(this.metaChunks),this.findIhdr(),this.parseTextChunks(),await this.findExif().catch(this.catchError),await this.findXmp().catch(this.catchError),await this.findIcc().catch(this.catchError)}async findPngChunksInRange(e,t){let{file:n}=this;for(;e<t;){let i=n.getUint32(e),s=n.getUint32(e+4),a=n.getString(e+4,4).toLowerCase(),o=i+4+4+4,c={type:a,offset:e,length:o,start:e+4+4,size:i,marker:s};Mn.includes(a)?this.metaChunks.push(c):this.unknownChunks.push(c),e+=o}}parseTextChunks(){let e=this.metaChunks.filter(t=>t.type==="text");for(let t of e){let[n,i]=this.file.getString(t.start,t.size).split("\0");this.injectKeyValToIhdr(n,i)}}injectKeyValToIhdr(e,t){let n=this.parsers.ihdr;n&&n.raw.set(e,t)}findIhdr(){let e=this.metaChunks.find(t=>t.type==="ihdr");e&&this.options.ihdr.enabled!==!1&&this.createParser("ihdr",e.chunk)}async findExif(){let e=this.metaChunks.find(t=>t.type==="exif");e&&this.injectSegment("tiff",e.chunk)}async findXmp(){let e=this.metaChunks.filter(t=>t.type==="itxt");for(let t of e)t.chunk.getString(0,17)==="XML:com.adobe.xmp"&&this.injectSegment("xmp",t.chunk)}async findIcc(){let e=this.metaChunks.find(o=>o.type==="iccp");if(!e)return;let{chunk:t}=e,n=t.getUint8Array(0,81),i=0;for(;i<80&&n[i]!==0;)i++;let s=i+2,a=t.getString(0,i);if(this.injectKeyValToIhdr("ProfileName",a),Pe){let o=await Dn,c=t.getUint8Array(s);c=o.inflateSync(c),this.injectSegment("icc",c)}}}g(Ct,"type","png"),G.set("png",Ct),v(T,"interop",[[1,"InteropIndex"],[2,"InteropVersion"],[4096,"RelatedImageFileFormat"],[4097,"RelatedImageWidth"],[4098,"RelatedImageHeight"]]),re(T,"ifd0",[[11,"ProcessingSoftware"],[254,"SubfileType"],[255,"OldSubfileType"],[263,"Thresholding"],[264,"CellWidth"],[265,"CellLength"],[266,"FillOrder"],[269,"DocumentName"],[280,"MinSampleValue"],[281,"MaxSampleValue"],[285,"PageName"],[286,"XPosition"],[287,"YPosition"],[290,"GrayResponseUnit"],[297,"PageNumber"],[321,"HalftoneHints"],[322,"TileWidth"],[323,"TileLength"],[332,"InkSet"],[337,"TargetPrinter"],[18246,"Rating"],[18249,"RatingPercent"],[33550,"PixelScale"],[34264,"ModelTransform"],[34377,"PhotoshopSettings"],[50706,"DNGVersion"],[50707,"DNGBackwardVersion"],[50708,"UniqueCameraModel"],[50709,"LocalizedCameraModel"],[50736,"DNGLensInfo"],[50739,"ShadowScale"],[50740,"DNGPrivateData"],[33920,"IntergraphMatrix"],[33922,"ModelTiePoint"],[34118,"SEMInfo"],[34735,"GeoTiffDirectory"],[34736,"GeoTiffDoubleParams"],[34737,"GeoTiffAsciiParams"],[50341,"PrintIM"],[50721,"ColorMatrix1"],[50722,"ColorMatrix2"],[50723,"CameraCalibration1"],[50724,"CameraCalibration2"],[50725,"ReductionMatrix1"],[50726,"ReductionMatrix2"],[50727,"AnalogBalance"],[50728,"AsShotNeutral"],[50729,"AsShotWhiteXY"],[50730,"BaselineExposure"],[50731,"BaselineNoise"],[50732,"BaselineSharpness"],[50734,"LinearResponseLimit"],[50735,"CameraSerialNumber"],[50741,"MakerNoteSafety"],[50778,"CalibrationIlluminant1"],[50779,"CalibrationIlluminant2"],[50781,"RawDataUniqueID"],[50827,"OriginalRawFileName"],[50828,"OriginalRawFileData"],[50831,"AsShotICCProfile"],[50832,"AsShotPreProfileMatrix"],[50833,"CurrentICCProfile"],[50834,"CurrentPreProfileMatrix"],[50879,"ColorimetricReference"],[50885,"SRawType"],[50898,"PanasonicTitle"],[50899,"PanasonicTitle2"],[50931,"CameraCalibrationSig"],[50932,"ProfileCalibrationSig"],[50933,"ProfileIFD"],[50934,"AsShotProfileName"],[50936,"ProfileName"],[50937,"ProfileHueSatMapDims"],[50938,"ProfileHueSatMapData1"],[50939,"ProfileHueSatMapData2"],[50940,"ProfileToneCurve"],[50941,"ProfileEmbedPolicy"],[50942,"ProfileCopyright"],[50964,"ForwardMatrix1"],[50965,"ForwardMatrix2"],[50966,"PreviewApplicationName"],[50967,"PreviewApplicationVersion"],[50968,"PreviewSettingsName"],[50969,"PreviewSettingsDigest"],[50970,"PreviewColorSpace"],[50971,"PreviewDateTime"],[50972,"RawImageDigest"],[50973,"OriginalRawFileDigest"],[50981,"ProfileLookTableDims"],[50982,"ProfileLookTableData"],[51043,"TimeCodes"],[51044,"FrameRate"],[51058,"TStop"],[51081,"ReelName"],[51089,"OriginalDefaultFinalSize"],[51090,"OriginalBestQualitySize"],[51091,"OriginalDefaultCropSize"],[51105,"CameraLabel"],[51107,"ProfileHueSatMapEncoding"],[51108,"ProfileLookTableEncoding"],[51109,"BaselineExposureOffset"],[51110,"DefaultBlackRender"],[51111,"NewRawImageDigest"],[51112,"RawToPreviewGain"]]);let xt=[[273,"StripOffsets"],[279,"StripByteCounts"],[288,"FreeOffsets"],[289,"FreeByteCounts"],[291,"GrayResponseCurve"],[292,"T4Options"],[293,"T6Options"],[300,"ColorResponseUnit"],[320,"ColorMap"],[324,"TileOffsets"],[325,"TileByteCounts"],[326,"BadFaxLines"],[327,"CleanFaxData"],[328,"ConsecutiveBadFaxLines"],[330,"SubIFD"],[333,"InkNames"],[334,"NumberofInks"],[336,"DotRange"],[338,"ExtraSamples"],[339,"SampleFormat"],[340,"SMinSampleValue"],[341,"SMaxSampleValue"],[342,"TransferRange"],[343,"ClipPath"],[344,"XClipPathUnits"],[345,"YClipPathUnits"],[346,"Indexed"],[347,"JPEGTables"],[351,"OPIProxy"],[400,"GlobalParametersIFD"],[401,"ProfileType"],[402,"FaxProfile"],[403,"CodingMethods"],[404,"VersionYear"],[405,"ModeNumber"],[433,"Decode"],[434,"DefaultImageColor"],[435,"T82Options"],[437,"JPEGTables"],[512,"JPEGProc"],[515,"JPEGRestartInterval"],[517,"JPEGLosslessPredictors"],[518,"JPEGPointTransforms"],[519,"JPEGQTables"],[520,"JPEGDCTables"],[521,"JPEGACTables"],[559,"StripRowCounts"],[999,"USPTOMiscellaneous"],[18247,"XP_DIP_XML"],[18248,"StitchInfo"],[28672,"SonyRawFileType"],[28688,"SonyToneCurve"],[28721,"VignettingCorrection"],[28722,"VignettingCorrParams"],[28724,"ChromaticAberrationCorrection"],[28725,"ChromaticAberrationCorrParams"],[28726,"DistortionCorrection"],[28727,"DistortionCorrParams"],[29895,"SonyCropTopLeft"],[29896,"SonyCropSize"],[32781,"ImageID"],[32931,"WangTag1"],[32932,"WangAnnotation"],[32933,"WangTag3"],[32934,"WangTag4"],[32953,"ImageReferencePoints"],[32954,"RegionXformTackPoint"],[32955,"WarpQuadrilateral"],[32956,"AffineTransformMat"],[32995,"Matteing"],[32996,"DataType"],[32997,"ImageDepth"],[32998,"TileDepth"],[33300,"ImageFullWidth"],[33301,"ImageFullHeight"],[33302,"TextureFormat"],[33303,"WrapModes"],[33304,"FovCot"],[33305,"MatrixWorldToScreen"],[33306,"MatrixWorldToCamera"],[33405,"Model2"],[33421,"CFARepeatPatternDim"],[33422,"CFAPattern2"],[33423,"BatteryLevel"],[33424,"KodakIFD"],[33445,"MDFileTag"],[33446,"MDScalePixel"],[33447,"MDColorTable"],[33448,"MDLabName"],[33449,"MDSampleInfo"],[33450,"MDPrepDate"],[33451,"MDPrepTime"],[33452,"MDFileUnits"],[33589,"AdventScale"],[33590,"AdventRevision"],[33628,"UIC1Tag"],[33629,"UIC2Tag"],[33630,"UIC3Tag"],[33631,"UIC4Tag"],[33918,"IntergraphPacketData"],[33919,"IntergraphFlagRegisters"],[33921,"INGRReserved"],[34016,"Site"],[34017,"ColorSequence"],[34018,"IT8Header"],[34019,"RasterPadding"],[34020,"BitsPerRunLength"],[34021,"BitsPerExtendedRunLength"],[34022,"ColorTable"],[34023,"ImageColorIndicator"],[34024,"BackgroundColorIndicator"],[34025,"ImageColorValue"],[34026,"BackgroundColorValue"],[34027,"PixelIntensityRange"],[34028,"TransparencyIndicator"],[34029,"ColorCharacterization"],[34030,"HCUsage"],[34031,"TrapIndicator"],[34032,"CMYKEquivalent"],[34152,"AFCP_IPTC"],[34232,"PixelMagicJBIGOptions"],[34263,"JPLCartoIFD"],[34306,"WB_GRGBLevels"],[34310,"LeafData"],[34687,"TIFF_FXExtensions"],[34688,"MultiProfiles"],[34689,"SharedData"],[34690,"T88Options"],[34732,"ImageLayer"],[34750,"JBIGOptions"],[34856,"Opto-ElectricConvFactor"],[34857,"Interlace"],[34908,"FaxRecvParams"],[34909,"FaxSubAddress"],[34910,"FaxRecvTime"],[34929,"FedexEDR"],[34954,"LeafSubIFD"],[37387,"FlashEnergy"],[37388,"SpatialFrequencyResponse"],[37389,"Noise"],[37390,"FocalPlaneXResolution"],[37391,"FocalPlaneYResolution"],[37392,"FocalPlaneResolutionUnit"],[37397,"ExposureIndex"],[37398,"TIFF-EPStandardID"],[37399,"SensingMethod"],[37434,"CIP3DataFile"],[37435,"CIP3Sheet"],[37436,"CIP3Side"],[37439,"StoNits"],[37679,"MSDocumentText"],[37680,"MSPropertySetStorage"],[37681,"MSDocumentTextPosition"],[37724,"ImageSourceData"],[40965,"InteropIFD"],[40976,"SamsungRawPointersOffset"],[40977,"SamsungRawPointersLength"],[41217,"SamsungRawByteOrder"],[41218,"SamsungRawUnknown"],[41484,"SpatialFrequencyResponse"],[41485,"Noise"],[41489,"ImageNumber"],[41490,"SecurityClassification"],[41491,"ImageHistory"],[41494,"TIFF-EPStandardID"],[41995,"DeviceSettingDescription"],[42112,"GDALMetadata"],[42113,"GDALNoData"],[44992,"ExpandSoftware"],[44993,"ExpandLens"],[44994,"ExpandFilm"],[44995,"ExpandFilterLens"],[44996,"ExpandScanner"],[44997,"ExpandFlashLamp"],[46275,"HasselbladRawImage"],[48129,"PixelFormat"],[48130,"Transformation"],[48131,"Uncompressed"],[48132,"ImageType"],[48256,"ImageWidth"],[48257,"ImageHeight"],[48258,"WidthResolution"],[48259,"HeightResolution"],[48320,"ImageOffset"],[48321,"ImageByteCount"],[48322,"AlphaOffset"],[48323,"AlphaByteCount"],[48324,"ImageDataDiscard"],[48325,"AlphaDataDiscard"],[50215,"OceScanjobDesc"],[50216,"OceApplicationSelector"],[50217,"OceIDNumber"],[50218,"OceImageLogic"],[50255,"Annotations"],[50459,"HasselbladExif"],[50547,"OriginalFileName"],[50560,"USPTOOriginalContentType"],[50656,"CR2CFAPattern"],[50710,"CFAPlaneColor"],[50711,"CFALayout"],[50712,"LinearizationTable"],[50713,"BlackLevelRepeatDim"],[50714,"BlackLevel"],[50715,"BlackLevelDeltaH"],[50716,"BlackLevelDeltaV"],[50717,"WhiteLevel"],[50718,"DefaultScale"],[50719,"DefaultCropOrigin"],[50720,"DefaultCropSize"],[50733,"BayerGreenSplit"],[50737,"ChromaBlurRadius"],[50738,"AntiAliasStrength"],[50752,"RawImageSegmentation"],[50780,"BestQualityScale"],[50784,"AliasLayerMetadata"],[50829,"ActiveArea"],[50830,"MaskedAreas"],[50935,"NoiseReductionApplied"],[50974,"SubTileBlockSize"],[50975,"RowInterleaveFactor"],[51008,"OpcodeList1"],[51009,"OpcodeList2"],[51022,"OpcodeList3"],[51041,"NoiseProfile"],[51114,"CacheVersion"],[51125,"DefaultUserCrop"],[51157,"NikonNEFInfo"],[65024,"KdcIFD"]];re(T,"ifd0",xt),re(T,"exif",xt),v(z,"gps",[[23,{M:"Magnetic North",T:"True North"}],[25,{K:"Kilometers",M:"Miles",N:"Nautical Miles"}]]);class We extends B{static canHandle(e,t){return e.getUint8(t+1)===224&&e.getUint32(t+4)===1246120262&&e.getUint8(t+8)===0}parse(){return this.parseTags(),this.translate(),this.output}parseTags(){this.raw=new Map([[0,this.chunk.getUint16(0)],[2,this.chunk.getUint8(2)],[3,this.chunk.getUint16(3)],[5,this.chunk.getUint16(5)],[7,this.chunk.getUint8(7)],[8,this.chunk.getUint8(8)]])}}g(We,"type","jfif"),g(We,"headerLength",9),P.set("jfif",We),v(T,"jfif",[[0,"JFIFVersion"],[2,"ResolutionUnit"],[3,"XResolution"],[5,"YResolution"],[7,"ThumbnailWidth"],[8,"ThumbnailHeight"]]);class vt extends B{parse(){return this.parseTags(),this.translate(),this.output}parseTags(){this.raw=new Map([[0,this.chunk.getUint32(0)],[4,this.chunk.getUint32(4)],[8,this.chunk.getUint8(8)],[9,this.chunk.getUint8(9)],[10,this.chunk.getUint8(10)],[11,this.chunk.getUint8(11)],[12,this.chunk.getUint8(12)],...Array.from(this.raw)])}}g(vt,"type","ihdr"),P.set("ihdr",vt),v(T,"ihdr",[[0,"ImageWidth"],[4,"ImageHeight"],[8,"BitDepth"],[9,"ColorType"],[10,"Compression"],[11,"Filter"],[12,"Interlace"]]),v(z,"ihdr",[[9,{0:"Grayscale",2:"RGB",3:"Palette",4:"Grayscale with Alpha",6:"RGB with Alpha",DEFAULT:"Unknown"}],[10,{0:"Deflate/Inflate",DEFAULT:"Unknown"}],[11,{0:"Adaptive",DEFAULT:"Unknown"}],[12,{0:"Noninterlaced",1:"Adam7 Interlace",DEFAULT:"Unknown"}]]);class Ce extends B{static canHandle(e,t){return e.getUint8(t+1)===226&&e.getUint32(t+4)===1229144927}static findPosition(e,t){let n=super.findPosition(e,t);return n.chunkNumber=e.getUint8(t+16),n.chunkCount=e.getUint8(t+17),n.multiSegment=n.chunkCount>1,n}static handleMultiSegments(e){return function(t){let n=function(i){let s=i[0].constructor,a=0;for(let u of i)a+=u.length;let o=new s(a),c=0;for(let u of i)o.set(u,c),c+=u.length;return o}(t.map(i=>i.chunk.toUint8()));return new L(n)}(e)}parse(){return this.raw=new Map,this.parseHeader(),this.parseTags(),this.translate(),this.output}parseHeader(){let{raw:e}=this;this.chunk.byteLength<84&&I("ICC header is too short");for(let[t,n]of Object.entries(Fn)){t=parseInt(t,10);let i=n(this.chunk,t);i!=="\0\0\0\0"&&e.set(t,i)}}parseTags(){let e,t,n,i,s,{raw:a}=this,o=this.chunk.getUint32(128),c=132,u=this.chunk.byteLength;for(;o--;){if(e=this.chunk.getString(c,4),t=this.chunk.getUint32(c+4),n=this.chunk.getUint32(c+8),i=this.chunk.getString(t,4),t+n>u)return void console.warn("reached the end of the first ICC chunk. Enable options.tiff.multiSegment to read all ICC segments.");s=this.parseTag(i,t,n),s!==void 0&&s!=="\0\0\0\0"&&a.set(e,s),c+=12}}parseTag(e,t,n){switch(e){case"desc":return this.parseDesc(t);case"mluc":return this.parseMluc(t);case"text":return this.parseText(t,n);case"sig ":return this.parseSig(t)}if(!(t+n>this.chunk.byteLength))return this.chunk.getUint8Array(t,n)}parseDesc(e){let t=this.chunk.getUint32(e+8)-1;return Z(this.chunk.getString(e+12,t))}parseText(e,t){return Z(this.chunk.getString(e+8,t-8))}parseSig(e){return Z(this.chunk.getString(e+8,4))}parseMluc(e){let{chunk:t}=this,n=t.getUint32(e+8),i=t.getUint32(e+12),s=e+16,a=[];for(let o=0;o<n;o++){let c=t.getString(s+0,2),u=t.getString(s+2,2),d=t.getUint32(s+4),l=t.getUint32(s+8)+e,f=Z(t.getUnicodeString(l,d));a.push({lang:c,country:u,text:f}),s+=i}return n===1?a[0].text:a}translateValue(e,t){return typeof e=="string"?t[e]||t[e.toLowerCase()]||e:t[e]||e}}g(Ce,"type","icc"),g(Ce,"multiSegment",!0),g(Ce,"headerLength",18);const Fn={4:W,8:function(r,e){return[r.getUint8(e),r.getUint8(e+1)>>4,r.getUint8(e+1)%16].map(t=>t.toString(10)).join(".")},12:W,16:W,20:W,24:function(r,e){const t=r.getUint16(e),n=r.getUint16(e+2)-1,i=r.getUint16(e+4),s=r.getUint16(e+6),a=r.getUint16(e+8),o=r.getUint16(e+10);return new Date(Date.UTC(t,n,i,s,a,o))},36:W,40:W,48:W,52:W,64:(r,e)=>r.getUint32(e),80:W};function W(r,e){return Z(r.getString(e,4))}P.set("icc",Ce),v(T,"icc",[[4,"ProfileCMMType"],[8,"ProfileVersion"],[12,"ProfileClass"],[16,"ColorSpaceData"],[20,"ProfileConnectionSpace"],[24,"ProfileDateTime"],[36,"ProfileFileSignature"],[40,"PrimaryPlatform"],[44,"CMMFlags"],[48,"DeviceManufacturer"],[52,"DeviceModel"],[56,"DeviceAttributes"],[64,"RenderingIntent"],[68,"ConnectionSpaceIlluminant"],[80,"ProfileCreator"],[84,"ProfileID"],["Header","ProfileHeader"],["MS00","WCSProfiles"],["bTRC","BlueTRC"],["bXYZ","BlueMatrixColumn"],["bfd","UCRBG"],["bkpt","MediaBlackPoint"],["calt","CalibrationDateTime"],["chad","ChromaticAdaptation"],["chrm","Chromaticity"],["ciis","ColorimetricIntentImageState"],["clot","ColorantTableOut"],["clro","ColorantOrder"],["clrt","ColorantTable"],["cprt","ProfileCopyright"],["crdi","CRDInfo"],["desc","ProfileDescription"],["devs","DeviceSettings"],["dmdd","DeviceModelDesc"],["dmnd","DeviceMfgDesc"],["dscm","ProfileDescriptionML"],["fpce","FocalPlaneColorimetryEstimates"],["gTRC","GreenTRC"],["gXYZ","GreenMatrixColumn"],["gamt","Gamut"],["kTRC","GrayTRC"],["lumi","Luminance"],["meas","Measurement"],["meta","Metadata"],["mmod","MakeAndModel"],["ncl2","NamedColor2"],["ncol","NamedColor"],["ndin","NativeDisplayInfo"],["pre0","Preview0"],["pre1","Preview1"],["pre2","Preview2"],["ps2i","PS2RenderingIntent"],["ps2s","PostScript2CSA"],["psd0","PostScript2CRD0"],["psd1","PostScript2CRD1"],["psd2","PostScript2CRD2"],["psd3","PostScript2CRD3"],["pseq","ProfileSequenceDesc"],["psid","ProfileSequenceIdentifier"],["psvm","PS2CRDVMSize"],["rTRC","RedTRC"],["rXYZ","RedMatrixColumn"],["resp","OutputResponse"],["rhoc","ReflectionHardcopyOrigColorimetry"],["rig0","PerceptualRenderingIntentGamut"],["rig2","SaturationRenderingIntentGamut"],["rpoc","ReflectionPrintOutputColorimetry"],["sape","SceneAppearanceEstimates"],["scoe","SceneColorimetryEstimates"],["scrd","ScreeningDesc"],["scrn","Screening"],["targ","CharTarget"],["tech","Technology"],["vcgt","VideoCardGamma"],["view","ViewingConditions"],["vued","ViewingCondDesc"],["wtpt","MediaWhitePoint"]]);const Se={"4d2p":"Erdt Systems",AAMA:"Aamazing Technologies",ACER:"Acer",ACLT:"Acolyte Color Research",ACTI:"Actix Sytems",ADAR:"Adara Technology",ADBE:"Adobe",ADI:"ADI Systems",AGFA:"Agfa Graphics",ALMD:"Alps Electric",ALPS:"Alps Electric",ALWN:"Alwan Color Expertise",AMTI:"Amiable Technologies",AOC:"AOC International",APAG:"Apago",APPL:"Apple Computer",AST:"AST","AT&T":"AT&T",BAEL:"BARBIERI electronic",BRCO:"Barco NV",BRKP:"Breakpoint",BROT:"Brother",BULL:"Bull",BUS:"Bus Computer Systems","C-IT":"C-Itoh",CAMR:"Intel",CANO:"Canon",CARR:"Carroll Touch",CASI:"Casio",CBUS:"Colorbus PL",CEL:"Crossfield",CELx:"Crossfield",CGS:"CGS Publishing Technologies International",CHM:"Rochester Robotics",CIGL:"Colour Imaging Group, London",CITI:"Citizen",CL00:"Candela",CLIQ:"Color IQ",CMCO:"Chromaco",CMiX:"CHROMiX",COLO:"Colorgraphic Communications",COMP:"Compaq",COMp:"Compeq/Focus Technology",CONR:"Conrac Display Products",CORD:"Cordata Technologies",CPQ:"Compaq",CPRO:"ColorPro",CRN:"Cornerstone",CTX:"CTX International",CVIS:"ColorVision",CWC:"Fujitsu Laboratories",DARI:"Darius Technology",DATA:"Dataproducts",DCP:"Dry Creek Photo",DCRC:"Digital Contents Resource Center, Chung-Ang University",DELL:"Dell Computer",DIC:"Dainippon Ink and Chemicals",DICO:"Diconix",DIGI:"Digital","DL&C":"Digital Light & Color",DPLG:"Doppelganger",DS:"Dainippon Screen",DSOL:"DOOSOL",DUPN:"DuPont",EPSO:"Epson",ESKO:"Esko-Graphics",ETRI:"Electronics and Telecommunications Research Institute",EVER:"Everex Systems",EXAC:"ExactCODE",Eizo:"Eizo",FALC:"Falco Data Products",FF:"Fuji Photo Film",FFEI:"FujiFilm Electronic Imaging",FNRD:"Fnord Software",FORA:"Fora",FORE:"Forefront Technology",FP:"Fujitsu",FPA:"WayTech Development",FUJI:"Fujitsu",FX:"Fuji Xerox",GCC:"GCC Technologies",GGSL:"Global Graphics Software",GMB:"Gretagmacbeth",GMG:"GMG",GOLD:"GoldStar Technology",GOOG:"Google",GPRT:"Giantprint",GTMB:"Gretagmacbeth",GVC:"WayTech Development",GW2K:"Sony",HCI:"HCI",HDM:"Heidelberger Druckmaschinen",HERM:"Hermes",HITA:"Hitachi America",HP:"Hewlett-Packard",HTC:"Hitachi",HiTi:"HiTi Digital",IBM:"IBM",IDNT:"Scitex",IEC:"Hewlett-Packard",IIYA:"Iiyama North America",IKEG:"Ikegami Electronics",IMAG:"Image Systems",IMI:"Ingram Micro",INTC:"Intel",INTL:"N/A (INTL)",INTR:"Intra Electronics",IOCO:"Iocomm International Technology",IPS:"InfoPrint Solutions Company",IRIS:"Scitex",ISL:"Ichikawa Soft Laboratory",ITNL:"N/A (ITNL)",IVM:"IVM",IWAT:"Iwatsu Electric",Idnt:"Scitex",Inca:"Inca Digital Printers",Iris:"Scitex",JPEG:"Joint Photographic Experts Group",JSFT:"Jetsoft Development",JVC:"JVC Information Products",KART:"Scitex",KFC:"KFC Computek Components",KLH:"KLH Computers",KMHD:"Konica Minolta",KNCA:"Konica",KODA:"Kodak",KYOC:"Kyocera",Kart:"Scitex",LCAG:"Leica",LCCD:"Leeds Colour",LDAK:"Left Dakota",LEAD:"Leading Technology",LEXM:"Lexmark International",LINK:"Link Computer",LINO:"Linotronic",LITE:"Lite-On",Leaf:"Leaf",Lino:"Linotronic",MAGC:"Mag Computronic",MAGI:"MAG Innovision",MANN:"Mannesmann",MICN:"Micron Technology",MICR:"Microtek",MICV:"Microvitec",MINO:"Minolta",MITS:"Mitsubishi Electronics America",MITs:"Mitsuba",MNLT:"Minolta",MODG:"Modgraph",MONI:"Monitronix",MONS:"Monaco Systems",MORS:"Morse Technology",MOTI:"Motive Systems",MSFT:"Microsoft",MUTO:"MUTOH INDUSTRIES",Mits:"Mitsubishi Electric",NANA:"NANAO",NEC:"NEC",NEXP:"NexPress Solutions",NISS:"Nissei Sangyo America",NKON:"Nikon",NONE:"none",OCE:"Oce Technologies",OCEC:"OceColor",OKI:"Oki",OKID:"Okidata",OKIP:"Okidata",OLIV:"Olivetti",OLYM:"Olympus",ONYX:"Onyx Graphics",OPTI:"Optiquest",PACK:"Packard Bell",PANA:"Matsushita Electric Industrial",PANT:"Pantone",PBN:"Packard Bell",PFU:"PFU",PHIL:"Philips Consumer Electronics",PNTX:"HOYA",POne:"Phase One A/S",PREM:"Premier Computer Innovations",PRIN:"Princeton Graphic Systems",PRIP:"Princeton Publishing Labs",QLUX:"Hong Kong",QMS:"QMS",QPCD:"QPcard AB",QUAD:"QuadLaser",QUME:"Qume",RADI:"Radius",RDDx:"Integrated Color Solutions",RDG:"Roland DG",REDM:"REDMS Group",RELI:"Relisys",RGMS:"Rolf Gierling Multitools",RICO:"Ricoh",RNLD:"Edmund Ronald",ROYA:"Royal",RPC:"Ricoh Printing Systems",RTL:"Royal Information Electronics",SAMP:"Sampo",SAMS:"Samsung",SANT:"Jaime Santana Pomares",SCIT:"Scitex",SCRN:"Dainippon Screen",SDP:"Scitex",SEC:"Samsung",SEIK:"Seiko Instruments",SEIk:"Seikosha",SGUY:"ScanGuy.com",SHAR:"Sharp Laboratories",SICC:"International Color Consortium",SONY:"Sony",SPCL:"SpectraCal",STAR:"Star",STC:"Sampo Technology",Scit:"Scitex",Sdp:"Scitex",Sony:"Sony",TALO:"Talon Technology",TAND:"Tandy",TATU:"Tatung",TAXA:"TAXAN America",TDS:"Tokyo Denshi Sekei",TECO:"TECO Information Systems",TEGR:"Tegra",TEKT:"Tektronix",TI:"Texas Instruments",TMKR:"TypeMaker",TOSB:"Toshiba",TOSH:"Toshiba",TOTK:"TOTOKU ELECTRIC",TRIU:"Triumph",TSBT:"Toshiba",TTX:"TTX Computer Products",TVM:"TVM Professional Monitor",TW:"TW Casper",ULSX:"Ulead Systems",UNIS:"Unisys",UTZF:"Utz Fehlau & Sohn",VARI:"Varityper",VIEW:"Viewsonic",VISL:"Visual communication",VIVO:"Vivo Mobile Communication",WANG:"Wang",WLBR:"Wilbur Imaging",WTG2:"Ware To Go",WYSE:"WYSE Technology",XERX:"Xerox",XRIT:"X-Rite",ZRAN:"Zoran",Zebr:"Zebra Technologies",appl:"Apple Computer",bICC:"basICColor",berg:"bergdesign",ceyd:"Integrated Color Solutions",clsp:"MacDermid ColorSpan",ds:"Dainippon Screen",dupn:"DuPont",ffei:"FujiFilm Electronic Imaging",flux:"FluxData",iris:"Scitex",kart:"Scitex",lcms:"Little CMS",lino:"Linotronic",none:"none",ob4d:"Erdt Systems",obic:"Medigraph",quby:"Qubyx Sarl",scit:"Scitex",scrn:"Dainippon Screen",sdp:"Scitex",siwi:"SIWI GRAFIKA",yxym:"YxyMaster"},Pt={scnr:"Scanner",mntr:"Monitor",prtr:"Printer",link:"Device Link",abst:"Abstract",spac:"Color Space Conversion Profile",nmcl:"Named Color",cenc:"ColorEncodingSpace profile",mid:"MultiplexIdentification profile",mlnk:"MultiplexLink profile",mvis:"MultiplexVisualization profile",nkpf:"Nikon Input Device Profile (NON-STANDARD!)"};v(z,"icc",[[4,Se],[12,Pt],[40,Object.assign({},Se,Pt)],[48,Se],[80,Se],[64,{0:"Perceptual",1:"Relative Colorimetric",2:"Saturation",3:"Absolute Colorimetric"}],["tech",{amd:"Active Matrix Display",crt:"Cathode Ray Tube Display",kpcd:"Photo CD",pmd:"Passive Matrix Display",dcam:"Digital Camera",dcpj:"Digital Cinema Projector",dmpc:"Digital Motion Picture Camera",dsub:"Dye Sublimation Printer",epho:"Electrophotographic Printer",esta:"Electrostatic Printer",flex:"Flexography",fprn:"Film Writer",fscn:"Film Scanner",grav:"Gravure",ijet:"Ink Jet Printer",imgs:"Photo Image Setter",mpfr:"Motion Picture Film Recorder",mpfs:"Motion Picture Film Scanner",offs:"Offset Lithography",pjtv:"Projection Television",rpho:"Photographic Paper Printer",rscn:"Reflective Scanner",silk:"Silkscreen",twax:"Thermal Wax Printer",vidc:"Video Camera",vidm:"Video Monitor"}]]);class ke extends B{static canHandle(e,t,n){return e.getUint8(t+1)===237&&e.getString(t+4,9)==="Photoshop"&&this.containsIptc8bim(e,t,n)!==void 0}static headerLength(e,t,n){let i,s=this.containsIptc8bim(e,t,n);if(s!==void 0)return i=e.getUint8(t+s+7),i%2!=0&&(i+=1),i===0&&(i=4),s+8+i}static containsIptc8bim(e,t,n){for(let i=0;i<n;i++)if(this.isIptcSegmentHead(e,t+i))return i}static isIptcSegmentHead(e,t){return e.getUint8(t)===56&&e.getUint32(t)===943868237&&e.getUint16(t+4)===1028}parse(){let{raw:e}=this,t=this.chunk.byteLength-1,n=!1;for(let i=0;i<t;i++)if(this.chunk.getUint8(i)===28&&this.chunk.getUint8(i+1)===2){n=!0;let s=this.chunk.getUint16(i+3),a=this.chunk.getUint8(i+2),o=this.chunk.getLatin1String(i+5,s);e.set(a,this.pluralizeValue(e.get(a),o)),i+=4+s}else if(n)break;return this.translate(),this.output}pluralizeValue(e,t){return e!==void 0?e instanceof Array?(e.push(t),e):[e,t]:t}}g(ke,"type","iptc"),g(ke,"translateValues",!1),g(ke,"reviveValues",!1),P.set("iptc",ke),v(T,"iptc",[[0,"ApplicationRecordVersion"],[3,"ObjectTypeReference"],[4,"ObjectAttributeReference"],[5,"ObjectName"],[7,"EditStatus"],[8,"EditorialUpdate"],[10,"Urgency"],[12,"SubjectReference"],[15,"Category"],[20,"SupplementalCategories"],[22,"FixtureIdentifier"],[25,"Keywords"],[26,"ContentLocationCode"],[27,"ContentLocationName"],[30,"ReleaseDate"],[35,"ReleaseTime"],[37,"ExpirationDate"],[38,"ExpirationTime"],[40,"SpecialInstructions"],[42,"ActionAdvised"],[45,"ReferenceService"],[47,"ReferenceDate"],[50,"ReferenceNumber"],[55,"DateCreated"],[60,"TimeCreated"],[62,"DigitalCreationDate"],[63,"DigitalCreationTime"],[65,"OriginatingProgram"],[70,"ProgramVersion"],[75,"ObjectCycle"],[80,"Byline"],[85,"BylineTitle"],[90,"City"],[92,"Sublocation"],[95,"State"],[100,"CountryCode"],[101,"Country"],[103,"OriginalTransmissionReference"],[105,"Headline"],[110,"Credit"],[115,"Source"],[116,"CopyrightNotice"],[118,"Contact"],[120,"Caption"],[121,"LocalCaption"],[122,"Writer"],[125,"RasterizedCaption"],[130,"ImageType"],[131,"ImageOrientation"],[135,"LanguageIdentifier"],[150,"AudioType"],[151,"AudioSamplingRate"],[152,"AudioSamplingResolution"],[153,"AudioDuration"],[154,"AudioOutcue"],[184,"JobID"],[185,"MasterDocumentID"],[186,"ShortDocumentID"],[187,"UniqueDocumentID"],[188,"OwnerID"],[200,"ObjectPreviewFileFormat"],[201,"ObjectPreviewFileVersion"],[202,"ObjectPreviewData"],[221,"Prefs"],[225,"ClassifyState"],[228,"SimilarityIndex"],[230,"DocumentNotes"],[231,"DocumentHistory"],[232,"ExifCameraInfo"],[255,"CatalogSets"]]),v(z,"iptc",[[10,{0:"0 (reserved)",1:"1 (most urgent)",2:"2",3:"3",4:"4",5:"5 (normal urgency)",6:"6",7:"7",8:"8 (least urgent)",9:"9 (user-defined priority)"}],[75,{a:"Morning",b:"Both Morning and Evening",p:"Evening"}],[131,{L:"Landscape",P:"Portrait",S:"Square"}]]);let Xe=null;async function Ln(){if(!Xe)try{const n=(await import("./joraw-1Lq5hXK7.js")).default;if(typeof n!="function")throw new Error("JoRaw WASM import failed");const i=new URL("/assets/joraw-DraTMNgX.wasm",import.meta.url).href;Xe=n({locateFile:(s,a)=>s.endsWith("joraw.wasm")?i:a+s})}catch(t){throw console.error("Failed to load joraw.js:",t),t}const r=await Xe,e=r.LibRaw||r.JoRaw;if(!e)throw new Error("JoRaw class not found");return e}const En=async r=>{var n,i,s,a,o,c,u,d,l;const e=await Ln(),t=new e;try{if(await t.open(r,{}),typeof t.getRawImage!="function")throw new Error("WASM mismatch");const f=t.getRawImage();let h=new Uint16Array(f.data);const p=await t.metadata(!0);let y={...p};try{const m=await On.parse(r.buffer);m&&(y={...y,...m})}catch(m){console.warn("exifr parsing failed for RAW buffer",m)}const w=((n=p.idata)==null?void 0:n.filters)||0,b=((i=p.idata)==null?void 0:i.colors)||0,C=w===0&&b===3,R=w===9;let S=[0,0,0,0],k=!1;if(t.getBlackLevels)try{const m=t.getBlackLevels();m.dng_cblack&&m.dng_cblack.length===4&&Array.from(m.dng_cblack).some(O=>O>0)?(S=Array.from(m.dng_cblack).map(Number),k=!0):m.cblack&&m.cblack.length===4&&Array.from(m.cblack).some(O=>O>0)?(S=Array.from(m.cblack).map(Number),k=!0):typeof m.black=="number"&&m.black>0&&(S=[m.black,m.black,m.black,m.black],k=!0)}catch(m){console.warn("getBlackLevels binding failed",m)}if(!k){let m=[];if((s=p.color_data)!=null&&s.cblack_rawpy_style)m=p.color_data.cblack_rawpy_style;else if((o=(a=p.color_data)==null?void 0:a.dng_levels)!=null&&o.dng_cblack)m=p.color_data.dng_levels.dng_cblack;else if(((c=p.black_level_per_channel)==null?void 0:c.length)>=4)m=p.black_level_per_channel;else if(((u=p.cblack)==null?void 0:u.length)>=4)m=p.cblack;else if(((l=(d=p.color)==null?void 0:d.cblack)==null?void 0:l.length)>=4)m=p.color.cblack;else{const O=p.black_level||p.color_data&&p.color_data.black||0;m=[O,O,O,O]}S=[Number(m[0])||0,Number(m[1])||0,Number(m[2])||0,Number(m[3])||0]}return{data:h,width:f.width,height:f.height,bayerPattern:p.color_desc||"RGGB",blackLevels:S,whiteLevel:p.white_level||16383,metadata:y,isThreePlane:C,threePlaneTransfer:C?"linear":void 0,isXTrans:R}}finally{t.delete?t.delete():t.close()}};async function Un(r){if(on(r)){const n=await cn(r);if(!n)throw new Error("Sony cRAW HQ decoder did not return image data.");return n.rawImageData}const t=new Uint8Array(r);return En(t)}function Nn(r,e,t,n,i,s=!1,a=!1,o={rows:!1,rowGroup:0,cols:!1,colGroup:0},c=!0,u=1,d=!1){const l=Math.max(0,Math.floor(n.x)),f=Math.max(0,Math.floor(n.y)),h=Math.min(e,Math.ceil(n.x+n.w)),p=Math.min(t,Math.ceil(n.y+n.h)),y=h-l,w=p-f;if(y<=1||w<=1)return null;const b=d&&r.length>=e*t*3;if(b){if(i<0||i>2)return null}else if(i<0||i>3)return null;let C=0,R=0;if(b)C=h-l,R=p-f;else{for(let _=l;_<h;_++)_%2===(i&1)&&C++;for(let _=f;_<p;_++)_%2===i>>1&&R++}if(C<4||R<4)return null;const S=Math.max(1,u),k=Math.floor(2*C/(S+1)),m=Math.floor(2*R/(S+1));if(k<4||m<4)return null;const O=new ve(k),fe=new ve(m),V=new Float32Array(Math.floor(k/2)+1),oe=new Float32Array(Math.floor(m/2)+1);let le=0,N=0;const ye=new Float32Array(C),ce=new Float32Array(R);for(let _=f;_<p;_++){const E=_%2;if(!b&&i>>1!==E||o.rows&&(b?_&1:_>>1&1)!==o.rowGroup)continue;const Q=_*e;let D=0;for(let A=l;A<h;A++){if(!b&&A%2!==(i&1)||D>=C)continue;const j=b?r[(Q+A)*3+i]:r[Q+A],be=Number.isFinite(j)?j:0;ye[D++]=be}for(let A=0;A<S;A++){const j=Math.floor(A*k/2);if(j+k>C)break;O.calculateSpectrumWindow(ye,j,V,a),le++}if(!c)break}for(let _=l;_<h;_++){const E=_%2;if(!b&&(i&1)!==E||o.cols&&(b?_&1:_>>1&1)!==o.colGroup)continue;let Q=0;for(let D=f;D<p;D++){if(!b&&D%2!==i>>1||Q>=R)continue;const A=b?r[(D*e+_)*3+i]:r[D*e+_],j=Number.isFinite(A)?A:0;ce[Q++]=j}for(let D=0;D<S;D++){const A=Math.floor(D*m/2);if(A+m>R)break;fe.calculateSpectrumWindow(ce,A,oe,a),N++}if(!c)break}if(le>0)for(let _=0;_<V.length;_++)V[_]/=le;if(N>0)for(let _=0;_<oe.length;_++)oe[_]/=N;return{avgHoriz:Array.from(V),avgVert:Array.from(oe),sampleCount:{h:le,v:N}}}const It=192*1024*1024,Oe=new Map;let xe=0,At=Promise.resolve();const Bn=r=>{r.lastAccess=Date.now()},zn=(r=new Set)=>{if(xe<=It)return;const e=Array.from(Oe.entries()).filter(([t])=>!r.has(t)).sort((t,n)=>t[1].lastAccess-n[1].lastAccess);for(const[t,n]of e){if(xe<=It)break;xe-=n.sizeBytes,Oe.delete(t)}},Vn=r=>{var e;return r.data.byteLength+(((e=r.floatData)==null?void 0:e.byteLength)||0)},jt=async r=>{const e=Oe.get(r.key);if(e)return Bn(e),e.raw;if(!r.file)throw new Error(`Worker cache miss for ${r.key}`);const t=await r.file.arrayBuffer(),n=await Un(t),i={raw:n,sizeBytes:Vn(n),lastAccess:Date.now()};return Oe.set(r.key,i),xe+=i.sizeBytes,zn(new Set([r.key])),n},Gn=async r=>{for(const e of r.files)await jt(e);self.postMessage({id:r.id,type:"result",success:!0})},Hn=async r=>{const e=new Map;for(const i of r.files){const s=await jt(i);for(const a of r.channelIndices){const o=Nn(s.data,s.width,s.height,r.roi,a,!1,r.useHanning,r.interlace,r.doAverage,r.nSegments,!!s.isThreePlane);if(!o)continue;let c=e.get(a);if(!c)c={sumH:new Float32Array(o.avgHoriz.length),sumV:new Float32Array(o.avgVert.length),sampleCount:o.sampleCount,contributingFrames:0},e.set(a,c);else if(c.sumH.length!==o.avgHoriz.length||c.sumV.length!==o.avgVert.length)throw new Error("Stack FFT requires matching frame dimensions and channel layout.");for(let u=0;u<o.avgHoriz.length;u++)c.sumH[u]+=o.avgHoriz[u];for(let u=0;u<o.avgVert.length;u++)c.sumV[u]+=o.avgVert[u];c.sampleCount=o.sampleCount,c.contributingFrames+=1}self.postMessage({id:r.id,type:"progress",completed:1})}const t=Array.from(e.entries()).map(([i,s])=>({ch:i,sumH:s.sumH,sumV:s.sumV,sampleCount:s.sampleCount,contributingFrames:s.contributingFrames})),n=[];for(const i of t)n.push(i.sumH.buffer,i.sumV.buffer);self.postMessage({id:r.id,type:"result",success:!0,processedFrames:r.files.length,results:t},n)};self.onmessage=r=>{const e=r.data;At=At.then(async()=>{try{e.type==="preload"?await Gn(e):await Hn(e)}catch(t){self.postMessage({id:e.id,type:"result",success:!1,error:(t==null?void 0:t.message)||"Power spectrum worker failed"})}})};
